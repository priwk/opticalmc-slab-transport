#!/usr/bin/env python3
"""Build macroscopic ZnS scintillation source CSVs from Stage B alpha/Li steps."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


STEP_FIELDS = [
    "source_id",
    "source_event_uid",
    "eventID",
    "trackID",
    "stepID",
    "particle",
    "ratio_tag",
    "bn_wt",
    "zns_wt",
    "thickness_um",
    "placement_file",
    "surface_mode",
    "capture_x_um",
    "capture_y_um",
    "corr_x_um",
    "corr_y_um",
    "depth_um",
    "macro_anchor_x_um",
    "macro_anchor_y_um",
    "macro_anchor_z_um",
    "local_capture_x_um",
    "local_capture_y_um",
    "local_capture_z_um",
    "phase_pre",
    "phase_post",
    "src_x0_um",
    "src_y0_um",
    "src_z0_um",
    "src_x1_um",
    "src_y1_um",
    "src_z1_um",
    "src_mid_x_um",
    "src_mid_y_um",
    "src_mid_z_um",
    "step_len_um",
    "edep_keV",
    "visible_edep_keV",
    "n_photon_step",
    "wavelength_nm",
]


EVENT_FIELDS = [
    "source_event_uid",
    "eventID",
    "ratio_tag",
    "bn_wt",
    "zns_wt",
    "thickness_um",
    "placement_file",
    "surface_mode",
    "capture_x_um",
    "capture_y_um",
    "corr_x_um",
    "corr_y_um",
    "depth_um",
    "macro_anchor_x_um",
    "macro_anchor_y_um",
    "macro_anchor_z_um",
    "local_capture_x_um",
    "local_capture_y_um",
    "local_capture_z_um",
    "n_total_steps",
    "n_zns_steps",
    "total_edep_keV",
    "total_visible_edep_keV",
    "total_n_photon",
    "has_zns_edep",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Stage B alpha/Li step CSVs into macroscopic ZnS source "
            "steps and event-level light summaries."
        )
    )
    parser.add_argument("input_csv", type=Path, help="Raw *_alpha_li_steps.csv file")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. Defaults to inputs/generated_sources/<ratio_tag> when the "
            "input is under inputs/alpha_li_steps/<ratio_tag>, otherwise the input CSV directory."
        ),
    )
    parser.add_argument(
        "--xy-anchor-mode",
        choices=("capture", "corr"),
        default="capture",
        help=(
            "Use capture_x/y_um as the absolute capture-position anchor, or corr_x/y_um "
            "(capture_x/y_um - source_x/y_um) as the source-relative transverse anchor."
        ),
    )
    parser.add_argument(
        "--yield-zns-per-MeV",
        type=float,
        default=60000.0,
        help="ZnS(Ag) scintillation yield in photons/MeV used for n_photon_step.",
    )
    parser.add_argument(
        "--wavelength-nm",
        type=float,
        default=450.0,
        help="Source wavelength written to the source CSV.",
    )
    parser.add_argument(
        "--ratio-tag",
        default=None,
        help="Ratio tag to write; defaults to the parent directory name.",
    )
    parser.add_argument(
        "--thickness-um",
        type=float,
        default=None,
        help="Override thickness inferred from the input file name.",
    )
    parser.add_argument(
        "--z-tolerance-um",
        type=float,
        default=1.0e-6,
        help="Numerical tolerance for small z-boundary excursions.",
    )
    parser.add_argument(
        "--quench-model",
        choices=("none", "birks"),
        default="none",
        help="Visible-energy model. none keeps visible_edep_keV = edep_keV; birks applies a simple step Birks factor.",
    )
    parser.add_argument(
        "--birks-kb-um-per-keV",
        type=float,
        default=0.0,
        help="Birks kB in um/keV for --quench-model birks: visible = edep/(1+kB*dE/dx).",
    )
    return parser.parse_args()


def parse_float(row: Dict[str, str], name: str, default: float = 0.0) -> float:
    text = row.get(name, "")
    if text is None or text == "":
        return default
    return float(text)


def infer_thickness(path: Path) -> Optional[float]:
    match = re.match(r"^\s*([0-9]+(?:\.[0-9]+)?)", path.name)
    if not match:
        return None
    return float(match.group(1))


def format_number(value: object) -> str:
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.12g}"
    return "" if value is None else str(value)


def event_uid(ratio_tag: str, thickness_um: float, event_id: str) -> str:
    return f"{ratio_tag}|{format_number(thickness_um)}|{event_id}"


def default_output_dir(input_csv: Path, ratio_tag: str) -> Path:
    parts = list(input_csv.resolve().parts)
    lowered = [p.lower() for p in parts]
    if "inputs" in lowered and "alpha_li_steps" in lowered:
        inputs_index = lowered.index("inputs")
        return Path(*parts[: inputs_index + 1]) / "generated_sources" / ratio_tag
    return input_csv.parent


def visible_energy_keV(edep_keV: float, step_len_um: float, args: argparse.Namespace) -> float:
    if args.quench_model == "none":
        return edep_keV
    if args.quench_model == "birks":
        if args.birks_kb_um_per_keV < 0.0:
            raise ValueError("--birks-kb-um-per-keV must be non-negative")
        dedx = edep_keV / max(step_len_um, 1.0e-15)
        return edep_keV / (1.0 + args.birks_kb_um_per_keV * dedx)
    raise ValueError(f"Unsupported quench model: {args.quench_model}")


def first_event_record(
    row: Dict[str, str],
    ratio_tag: str,
    thickness_um: float,
    anchor_x: float,
    anchor_y: float,
) -> Dict[str, object]:
    return {
        "source_event_uid": event_uid(ratio_tag, thickness_um, row.get("eventID", "")),
        "eventID": row.get("eventID", ""),
        "ratio_tag": ratio_tag,
        "bn_wt": row.get("bn_wt", ""),
        "zns_wt": row.get("zns_wt", ""),
        "thickness_um": thickness_um,
        "placement_file": row.get("placement_file", ""),
        "surface_mode": row.get("surface_mode", ""),
        "capture_x_um": row.get("capture_x_um", ""),
        "capture_y_um": row.get("capture_y_um", ""),
        "corr_x_um": row.get("corr_x_um", ""),
        "corr_y_um": row.get("corr_y_um", ""),
        "depth_um": row.get("depth_um", ""),
        "macro_anchor_x_um": anchor_x,
        "macro_anchor_y_um": anchor_y,
        "macro_anchor_z_um": parse_float(row, "depth_um"),
        "local_capture_x_um": row.get("local_capture_x_um", ""),
        "local_capture_y_um": row.get("local_capture_y_um", ""),
        "local_capture_z_um": row.get("local_capture_z_um", ""),
        "n_total_steps": 0,
        "n_zns_steps": 0,
        "total_edep_keV": 0.0,
        "total_visible_edep_keV": 0.0,
        "total_n_photon": 0.0,
        "has_zns_edep": 0,
    }


def clip_segment_z(
    p0: Tuple[float, float, float],
    p1: Tuple[float, float, float],
    z_min: float,
    z_max: float,
    tol: float,
) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    x0, y0, z0 = p0
    x1, y1, z1 = p1

    if z0 < z_min - tol and z1 < z_min - tol:
        return None
    if z0 > z_max + tol and z1 > z_max + tol:
        return None

    if abs(z1 - z0) < 1.0e-15:
        z = min(max(z0, z_min), z_max)
        return (x0, y0, z), (x1, y1, z)

    t0 = 0.0
    t1 = 1.0
    dz = z1 - z0

    if z0 < z_min:
        t0 = max(t0, (z_min - z0) / dz)
    elif z0 > z_max:
        t0 = max(t0, (z_max - z0) / dz)

    if z1 < z_min:
        t1 = min(t1, (z_min - z0) / dz)
    elif z1 > z_max:
        t1 = min(t1, (z_max - z0) / dz)

    lo = min(t0, t1)
    hi = max(t0, t1)
    lo = min(max(lo, 0.0), 1.0)
    hi = min(max(hi, 0.0), 1.0)
    if hi < lo:
        return None

    def point_at(t: float) -> Tuple[float, float, float]:
        return (
            x0 + (x1 - x0) * t,
            y0 + (y1 - y0) * t,
            min(max(z0 + dz * t, z_min), z_max),
        )

    return point_at(lo), point_at(hi)


def segment_length(p0: Tuple[float, float, float], p1: Tuple[float, float, float]) -> float:
    return math.sqrt(
        (p1[0] - p0[0]) * (p1[0] - p0[0])
        + (p1[1] - p0[1]) * (p1[1] - p0[1])
        + (p1[2] - p0[2]) * (p1[2] - p0[2])
    )


def make_outputs(args: argparse.Namespace) -> Tuple[Path, Path, int, int]:
    input_csv = args.input_csv
    inferred_thickness = infer_thickness(input_csv)
    if args.thickness_um is not None:
        thickness_um = args.thickness_um
    elif inferred_thickness is not None:
        thickness_um = inferred_thickness
    else:
        raise ValueError("Could not infer thickness from file name; pass --thickness-um.")

    ratio_tag = args.ratio_tag or input_csv.parent.name
    output_dir = args.output_dir or default_output_dir(input_csv, ratio_tag)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem_thickness = format_number(thickness_um)
    step_out = output_dir / f"{stem_thickness}_macro_zns_step_sources.csv"
    event_out = output_dir / f"{stem_thickness}_event_light_sources.csv"

    events: "OrderedDict[str, Dict[str, object]]" = OrderedDict()
    source_count = 0
    row_count = 0
    column_thickness_values = set()

    with input_csv.open("r", newline="", encoding="utf-8-sig") as fin, step_out.open(
        "w", newline="", encoding="utf-8"
    ) as fout:
        reader = csv.DictReader(fin)
        missing = [
            c
            for c in (
                "eventID",
                "thickness_um",
                "phase_pre",
                "x_pre_um",
                "y_pre_um",
                "z_pre_um",
                "x_post_um",
                "y_post_um",
                "z_post_um",
                "edep_keV",
                "step_len_um",
                "depth_um",
                "local_capture_x_um",
                "local_capture_y_um",
                "local_capture_z_um",
            )
            if c not in (reader.fieldnames or [])
        ]
        if missing:
            raise ValueError(f"Input CSV is missing required columns: {', '.join(missing)}")

        writer = csv.DictWriter(fout, fieldnames=STEP_FIELDS, lineterminator="\n")
        writer.writeheader()

        for row in reader:
            row_count += 1
            try:
                column_thickness_values.add(round(parse_float(row, "thickness_um"), 9))
                macro_anchor_x = parse_float(
                    row, "capture_x_um" if args.xy_anchor_mode == "capture" else "corr_x_um"
                )
                macro_anchor_y = parse_float(
                    row, "capture_y_um" if args.xy_anchor_mode == "capture" else "corr_y_um"
                )
                uid = event_uid(ratio_tag, thickness_um, row.get("eventID", ""))
                if uid not in events:
                    events[uid] = first_event_record(
                        row, ratio_tag, thickness_um, macro_anchor_x, macro_anchor_y
                    )
                events[uid]["n_total_steps"] = int(events[uid]["n_total_steps"]) + 1

                edep_keV = parse_float(row, "edep_keV")
                input_step_len_um = parse_float(row, "step_len_um")
                if (
                    row.get("phase_pre", "") != "ZnS"
                    or edep_keV <= 0.0
                    or input_step_len_um <= 0.0
                ):
                    continue

                local_capture_x = parse_float(row, "local_capture_x_um")
                local_capture_y = parse_float(row, "local_capture_y_um")
                local_capture_z = parse_float(row, "local_capture_z_um")
                macro_anchor_z = parse_float(row, "depth_um")

                p0 = (
                    macro_anchor_x + (parse_float(row, "x_pre_um") - local_capture_x),
                    macro_anchor_y + (parse_float(row, "y_pre_um") - local_capture_y),
                    macro_anchor_z + (parse_float(row, "z_pre_um") - local_capture_z),
                )
                p1 = (
                    macro_anchor_x + (parse_float(row, "x_post_um") - local_capture_x),
                    macro_anchor_y + (parse_float(row, "y_post_um") - local_capture_y),
                    macro_anchor_z + (parse_float(row, "z_post_um") - local_capture_z),
                )
                clipped = clip_segment_z(
                    p0, p1, 0.0, thickness_um, max(0.0, args.z_tolerance_um)
                )
                if clipped is None:
                    continue
                p0c, p1c = clipped
                clipped_len = segment_length(p0c, p1c)
                if clipped_len <= 0.0:
                    continue

                visible_edep_keV = visible_energy_keV(edep_keV, clipped_len, args)
                n_photon_step = args.yield_zns_per_MeV * visible_edep_keV / 1000.0
                source_id = source_count
                source_count += 1

                events[uid]["n_zns_steps"] = int(events[uid]["n_zns_steps"]) + 1
                events[uid]["total_edep_keV"] = float(events[uid]["total_edep_keV"]) + edep_keV
                events[uid]["total_visible_edep_keV"] = (
                    float(events[uid]["total_visible_edep_keV"]) + visible_edep_keV
                )
                events[uid]["total_n_photon"] = (
                    float(events[uid]["total_n_photon"]) + n_photon_step
                )
                events[uid]["has_zns_edep"] = 1

                writer.writerow(
                    {
                        "source_id": source_id,
                        "source_event_uid": uid,
                        "eventID": row.get("eventID", ""),
                        "trackID": row.get("trackID", ""),
                        "stepID": row.get("stepID", ""),
                        "particle": row.get("particle", ""),
                        "ratio_tag": ratio_tag,
                        "bn_wt": row.get("bn_wt", ""),
                        "zns_wt": row.get("zns_wt", ""),
                        "thickness_um": thickness_um,
                        "placement_file": row.get("placement_file", ""),
                        "surface_mode": row.get("surface_mode", ""),
                        "capture_x_um": row.get("capture_x_um", ""),
                        "capture_y_um": row.get("capture_y_um", ""),
                        "corr_x_um": row.get("corr_x_um", ""),
                        "corr_y_um": row.get("corr_y_um", ""),
                        "depth_um": row.get("depth_um", ""),
                        "macro_anchor_x_um": macro_anchor_x,
                        "macro_anchor_y_um": macro_anchor_y,
                        "macro_anchor_z_um": macro_anchor_z,
                        "local_capture_x_um": local_capture_x,
                        "local_capture_y_um": local_capture_y,
                        "local_capture_z_um": local_capture_z,
                        "phase_pre": row.get("phase_pre", ""),
                        "phase_post": row.get("phase_post", ""),
                        "src_x0_um": p0c[0],
                        "src_y0_um": p0c[1],
                        "src_z0_um": p0c[2],
                        "src_x1_um": p1c[0],
                        "src_y1_um": p1c[1],
                        "src_z1_um": p1c[2],
                        "src_mid_x_um": 0.5 * (p0c[0] + p1c[0]),
                        "src_mid_y_um": 0.5 * (p0c[1] + p1c[1]),
                        "src_mid_z_um": 0.5 * (p0c[2] + p1c[2]),
                        "step_len_um": clipped_len,
                        "edep_keV": edep_keV,
                        "visible_edep_keV": visible_edep_keV,
                        "n_photon_step": n_photon_step,
                        "wavelength_nm": args.wavelength_nm,
                    }
                )
            except Exception as exc:
                raise ValueError(f"Failed while processing row {row_count}: {exc}") from exc

    if column_thickness_values:
        for column_thickness in sorted(column_thickness_values):
            if abs(column_thickness - thickness_um) > max(1.0e-6, args.z_tolerance_um):
                print(
                    "warning: filename/config thickness "
                    f"{format_number(thickness_um)} um differs from thickness_um column "
                    f"value {format_number(column_thickness)} um",
                    file=sys.stderr,
                )

    with event_out.open("w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=EVENT_FIELDS, lineterminator="\n")
        writer.writeheader()
        for record in events.values():
            writer.writerow({name: format_number(record.get(name, "")) for name in EVENT_FIELDS})

    return step_out, event_out, row_count, source_count


def main() -> int:
    args = parse_args()
    try:
        step_out, event_out, row_count, source_count = make_outputs(args)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"read_rows,{row_count}")
    print(f"zns_source_steps,{source_count}")
    print(f"wrote,{step_out}")
    print(f"wrote,{event_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
