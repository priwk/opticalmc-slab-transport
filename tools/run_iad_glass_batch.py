#!/usr/bin/env python3
"""Run IAD for manually entered samples mounted on a glass substrate.

The intended input is one CSV row per sample, for example:

    配比,厚度,总透射,总反射,g
    1-2,100,35.2,18.4,0.4

By default the sample thickness is interpreted as um, reflectance and
transmittance are auto-detected as either fractions or percentages, and the
sample is treated as mounted on one 1.1 mm bottom glass slide.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IAD_EXE_NAME = "iad.exe" if os.name == "nt" else "iad"
DEFAULT_IAD_EXE = (
    REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "src" / DEFAULT_IAD_EXE_NAME
)
DEFAULT_INPUT = REPO_ROOT / "inputs" / "iad_inputs" / "glass_1p1mm" / "samples.csv"
DEFAULT_OUTPUT = REPO_ROOT / "outputs" / "iad" / "glass_1p1mm" / "iad_results.csv"

NUMBER_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
RESULT_RE = re.compile(
    rf"^\s*"
    rf"(?P<wave>{NUMBER_RE})\s+"
    rf"(?P<mr>{NUMBER_RE})\s+"
    rf"(?P<mr_fit>{NUMBER_RE})\s+"
    rf"(?P<mt>{NUMBER_RE})\s+"
    rf"(?P<mt_fit>{NUMBER_RE})\s+"
    rf"(?P<mua>{NUMBER_RE})\s+"
    rf"(?P<musp>{NUMBER_RE})\s+"
    rf"(?P<g>{NUMBER_RE})\s+"
    rf"(?P<status>\S)"
)

ALIASES = {
    "sample_id": ["sample_id", "sample", "id", "编号", "样品", "样品编号"],
    "ratio": ["ratio", "ratio_tag", "formula", "mix", "配比", "比例"],
    "thickness": [
        "thickness",
        "thickness_um",
        "thickness_mm",
        "厚度",
        "厚度_um",
        "厚度_mm",
        "厚度um",
        "厚度mm",
        "样品厚度",
    ],
    "total_transmission": [
        "total_transmission",
        "total_transmission_percent",
        "total_transmittance",
        "total_transmittance_percent",
        "transmission",
        "transmission_percent",
        "transmittance",
        "transmittance_percent",
        "mt",
        "mt_percent",
        "t",
        "t_percent",
        "总透射",
        "总透射百分比",
        "透射",
        "透射百分比",
        "总透过率",
        "透过率",
    ],
    "total_reflection": [
        "total_reflection",
        "total_reflection_percent",
        "total_reflectance",
        "total_reflectance_percent",
        "reflection",
        "reflection_percent",
        "reflectance",
        "reflectance_percent",
        "mr",
        "mr_percent",
        "r",
        "r_percent",
        "总反射",
        "总反射百分比",
        "反射",
        "反射百分比",
        "总反射率",
        "反射率",
    ],
    "g": ["g", "anisotropy", "anisotropy_g", "各向异性", "各向异性因子"],
    "wavelength_nm": ["wavelength", "wavelength_nm", "lambda", "lambda_nm", "波长", "波长_nm"],
}

SLIDE_CONFIG = {
    "0": "0",
    "none": "0",
    "no": "0",
    "2": "2",
    "two": "2",
    "both": "2",
    "t": "t",
    "top": "t",
    "b": "b",
    "bottom": "b",
    "n": "n",
    "near": "n",
    "f": "f",
    "far": "f",
}


@dataclass
class Measurement:
    row_index: int
    sample_id: str
    ratio_tag: str
    sample_index: float
    sample_index_source: str
    thickness_input: str
    thickness_unit: str
    thickness_mm: float
    wavelength_nm: float
    mr_input: str
    mt_input: str
    mr: float
    mt: float
    g: float
    g_source: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Scott Prahl IAD on a simple sample table. Default geometry is one "
            "1.1 mm glass slide below the sample: -G b -N 1.52 -D 1.1."
        )
    )
    parser.add_argument(
        "input_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input CSV with ratio/thickness/T/R/g columns. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output CSV. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--effective-output-file",
        type=Path,
        default=None,
        help=(
            "Optional compact output CSV with effective-model column names "
            "(mua_eff_1_per_mm, musp_eff_1_per_mm, g_model)."
        ),
    )
    parser.add_argument(
        "--write-template",
        type=Path,
        default=None,
        help="Write an empty input template CSV and exit.",
    )
    parser.add_argument(
        "--iad-exe",
        type=Path,
        default=None,
        help="Path to iad.exe/iad. Default searches third_party/iad and PATH.",
    )
    parser.add_argument(
        "--thickness-unit",
        choices=["um", "mm"],
        default="um",
        help="Unit for an unqualified thickness column. Default: um.",
    )
    parser.add_argument(
        "--rt-unit",
        choices=["auto", "fraction", "percent"],
        default="auto",
        help="Unit for total reflection/transmission. Auto treats values > 1 as percent.",
    )
    parser.add_argument(
        "--wavelength-nm",
        type=float,
        default=450.0,
        help="Default wavelength if the CSV has no wavelength column. Default: 450 nm.",
    )
    parser.add_argument(
        "--g",
        type=float,
        default=None,
        help="Fallback anisotropy if the CSV has no g column.",
    )
    parser.add_argument(
        "--override-g",
        type=float,
        default=None,
        help=(
            "Force this anisotropy for every row, ignoring the CSV g column. "
            "Use --override-g 0 for the transport-equivalent IAD model."
        ),
    )
    parser.add_argument(
        "--sample-index",
        type=float,
        default=None,
        help=(
            "Override sample refractive index passed to IAD with -n. If omitted, "
            "n_eff is calculated from each row's ratio using the density/mixing model."
        ),
    )
    parser.add_argument(
        "--n-model",
        choices=("lorentz-lorenz", "linear"),
        default="lorentz-lorenz",
        help="Effective-index mixing rule used when --sample-index is omitted.",
    )
    parser.add_argument("--n-bn", type=float, default=1.80, help="BN refractive index.")
    parser.add_argument("--n-zns", type=float, default=2.36, help="ZnS refractive index.")
    parser.add_argument("--n-pmma", type=float, default=1.49, help="PMMA refractive index.")
    parser.add_argument("--n-air", type=float, default=1.000293, help="Air refractive index.")
    parser.add_argument("--rho-bn", type=float, default=2.10, help="BN density in g/cm^3.")
    parser.add_argument("--rho-zns", type=float, default=4.09, help="ZnS density in g/cm^3.")
    parser.add_argument(
        "--bn-zns-volume",
        type=float,
        default=64.0,
        help="Relative volume of BN/ZnS powder in the composite.",
    )
    parser.add_argument(
        "--pmma-volume",
        type=float,
        default=21.6,
        help="Relative volume of PMMA in the composite.",
    )
    parser.add_argument(
        "--air-volume",
        type=float,
        default=14.4,
        help="Relative volume of air/void in the composite.",
    )
    parser.add_argument(
        "--slide-configuration",
        default="bottom",
        choices=sorted(SLIDE_CONFIG.keys()),
        help=(
            "Slide geometry passed to IAD with -G. Default 'bottom' means light hits "
            "the sample first and the 1.1 mm glass is behind it."
        ),
    )
    parser.add_argument(
        "--slide-index",
        type=float,
        default=1.52,
        help="Glass refractive index passed to IAD with -N. Default: 1.52.",
    )
    parser.add_argument(
        "--slide-thickness-mm",
        type=float,
        default=1.1,
        help="Glass thickness passed to IAD with -D. Default: 1.1 mm.",
    )
    parser.add_argument(
        "--slide-od",
        type=float,
        default=None,
        help="Optional slide optical depth passed to IAD with -E.",
    )
    parser.add_argument(
        "--mc-iterations",
        type=int,
        default=0,
        help="IAD Monte Carlo correction iterations passed with -M. Default: 0.",
    )
    parser.add_argument(
        "--quadrature",
        type=int,
        default=None,
        help="Optional AD quadrature points passed to IAD with -q.",
    )
    parser.add_argument(
        "--incident-angle-deg",
        type=float,
        default=None,
        help="Optional incident angle passed to IAD with -i.",
    )
    parser.add_argument(
        "--extra-iad-args",
        default="",
        help='Extra raw IAD arguments, e.g. "-X -i 8".',
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write commands without running IAD.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop on the first invalid row or IAD error.",
    )
    return parser.parse_args()


def resolve_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def resolve_iad_exe(path: Optional[Path]) -> Path:
    candidates: List[Path] = []
    if path is not None:
        candidates.append(path)
    candidates.extend(
        [
            DEFAULT_IAD_EXE,
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / DEFAULT_IAD_EXE_NAME,
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "src" / DEFAULT_IAD_EXE_NAME,
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "iad.exe",
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "iad",
            REPO_ROOT / "third_party" / "iad" / "iad.exe",
            REPO_ROOT / "third_party" / "iad" / "iad",
        ]
    )
    for candidate in candidates:
        resolved = resolve_path(candidate)
        if resolved.exists():
            return resolved
    on_path = shutil.which("iad.exe") or shutil.which("iad")
    if on_path:
        return Path(on_path)
    checked = ", ".join(str(p) for p in candidates)
    raise FileNotFoundError(
        "Could not find IAD executable. Build it first or pass --iad-exe. "
        f"Checked: {checked}"
    )


def normalize_name(value: str) -> str:
    value = value.replace("\ufeff", "").replace("µ", "u").replace("μ", "u")
    value = value.strip().lower()
    return re.sub(r"[\s_\-./\\()[\]{}%％]+", "", value)


def canonical_ratio_tag(value: str) -> str:
    text = value.strip().replace("：", ":").replace("－", "-").replace("／", "/")
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*月\s*([0-9]+(?:\.[0-9]+)?)\s*日?", text)
    if match:
        return f"{match.group(1)}-{match.group(2)}"
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*/\s*([0-9]+(?:\.[0-9]+)?)", text)
    if match:
        return f"{match.group(1)}-{match.group(2)}"
    return text


def find_column(fieldnames: Sequence[str], aliases: Sequence[str]) -> Optional[str]:
    normalized = {normalize_name(name): name for name in fieldnames if name is not None}
    for alias in aliases:
        match = normalized.get(normalize_name(alias))
        if match is not None:
            return match
    return None


def parse_number(raw: str, label: str, row_index: int) -> Tuple[float, bool]:
    text = str(raw).strip().replace("％", "%")
    if not text:
        raise ValueError(f"row {row_index}: missing {label}")
    has_percent = text.endswith("%")
    if has_percent:
        text = text[:-1].strip()
    text = text.replace(",", "")
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(f"row {row_index}: invalid {label}: {raw!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"row {row_index}: non-finite {label}: {raw!r}")
    return value, has_percent


def as_fraction(raw: str, label: str, row_index: int, unit: str) -> float:
    value, has_percent = parse_number(raw, label, row_index)
    if has_percent or unit == "percent":
        value /= 100.0
    elif unit == "auto" and value > 1.0:
        value /= 100.0
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"row {row_index}: {label} must be between 0 and 1 after unit conversion")
    return value


def thickness_unit_from_column(column_name: str, default_unit: str) -> str:
    name = normalize_name(column_name)
    if "mm" in name or "毫米" in name:
        return "mm"
    if "um" in name or "微米" in name or "µm" in name or "μm" in name:
        return "um"
    return default_unit


def thickness_to_mm(raw: str, column_name: str, row_index: int, default_unit: str) -> Tuple[float, str]:
    value, _ = parse_number(raw, "thickness", row_index)
    unit = thickness_unit_from_column(column_name, default_unit)
    if value <= 0:
        raise ValueError(f"row {row_index}: thickness must be positive")
    return (value / 1000.0 if unit == "um" else value), unit


def validate_g(value: float, row_index: int) -> float:
    if not -1.0 <= value <= 1.0:
        raise ValueError(f"row {row_index}: g must be between -1 and 1")
    return value


def ratio_numbers(ratio: str, row_index: int) -> Tuple[float, float]:
    text = ratio.strip()
    if ":" in text:
        left, right = text.split(":", 1)
    elif "-" in text:
        left, right = text.split("-", 1)
    else:
        raise ValueError(
            f"row {row_index}: cannot parse ratio {ratio!r}; expected forms like 1-2 or 1:2.5"
        )
    try:
        bn_mass = float(left)
        zns_mass = float(right)
    except ValueError as exc:
        raise ValueError(
            f"row {row_index}: cannot parse ratio {ratio!r}; expected numeric BN:ZnS parts"
        ) from exc
    if bn_mass <= 0 or zns_mass <= 0:
        raise ValueError(f"row {row_index}: ratio parts must be positive")
    return bn_mass, zns_mass


def effective_index_from_ratio(ratio: str, row_index: int, args: argparse.Namespace) -> float:
    bn_mass, zns_mass = ratio_numbers(ratio, row_index)
    if args.rho_bn <= 0 or args.rho_zns <= 0:
        raise ValueError("material densities must be positive")
    total_volume = args.bn_zns_volume + args.pmma_volume + args.air_volume
    if total_volume <= 0:
        raise ValueError("composite component volumes must sum to a positive value")

    v_bn_raw = bn_mass / args.rho_bn
    v_zns_raw = zns_mass / args.rho_zns
    bn_fraction_inside = v_bn_raw / (v_bn_raw + v_zns_raw)
    zns_fraction_inside = 1.0 - bn_fraction_inside
    powder_fraction = args.bn_zns_volume / total_volume
    fractions = [
        (powder_fraction * bn_fraction_inside, args.n_bn),
        (powder_fraction * zns_fraction_inside, args.n_zns),
        (args.pmma_volume / total_volume, args.n_pmma),
        (args.air_volume / total_volume, args.n_air),
    ]
    if args.n_model == "linear":
        return sum(fraction * index for fraction, index in fractions)

    ll_sum = 0.0
    for fraction, index in fractions:
        if index <= 0:
            raise ValueError("refractive indices must be positive")
        eps = index * index
        ll_sum += fraction * (eps - 1.0) / (eps + 2.0)
    eps_eff = (1.0 + 2.0 * ll_sum) / max(1.0e-15, 1.0 - ll_sum)
    return math.sqrt(eps_eff)


def sample_index_for_row(ratio_tag: str, row_index: int, args: argparse.Namespace) -> Tuple[float, str]:
    if args.sample_index is not None:
        if args.sample_index <= 0:
            raise ValueError("--sample-index must be positive")
        return args.sample_index, "manual --sample-index"
    return effective_index_from_ratio(ratio_tag, row_index, args), (
        f"ratio_density_{args.n_model}"
    )


def read_measurements(path: Path, args: argparse.Namespace) -> List[Measurement]:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing input CSV: {path}\n"
            f"Create one with: python .\\tools\\run_iad_glass_batch.py --write-template {path}"
        )
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Input CSV has no header: {path}")
        fieldnames = list(reader.fieldnames)
        sample_id_col = find_column(fieldnames, ALIASES["sample_id"])
        ratio_col = find_column(fieldnames, ALIASES["ratio"])
        thickness_col = find_column(fieldnames, ALIASES["thickness"])
        mt_col = find_column(fieldnames, ALIASES["total_transmission"])
        mr_col = find_column(fieldnames, ALIASES["total_reflection"])
        g_col = find_column(fieldnames, ALIASES["g"])
        wavelength_col = find_column(fieldnames, ALIASES["wavelength_nm"])

        missing = []
        if ratio_col is None and args.sample_index is None:
            missing.append("ratio/配比")
        elif sample_id_col is None and ratio_col is None:
            missing.append("ratio/配比 or sample_id/样品")
        for label, column in [
            ("thickness/厚度", thickness_col),
            ("total_transmission/总透射", mt_col),
            ("total_reflection/总反射", mr_col),
        ]:
            if column is None:
                missing.append(label)
        if g_col is None and args.g is None and args.override_g is None:
            missing.append("g")
        if missing:
            raise ValueError(f"Input CSV is missing required column(s): {', '.join(missing)}")

        measurements: List[Measurement] = []
        for row_index, row in enumerate(reader, start=2):
            if not any((value or "").strip() for value in row.values()):
                continue
            assert thickness_col is not None
            assert mt_col is not None
            assert mr_col is not None

            ratio_tag = canonical_ratio_tag((row.get(ratio_col or "", "") or "").strip())
            thickness_input = (row.get(thickness_col, "") or "").strip()
            sample_id = (row.get(sample_id_col or "", "") or "").strip()
            if not sample_id:
                sample_id = f"{ratio_tag}-{thickness_input}"
            if not ratio_tag:
                ratio_tag = sample_id
            sample_index, sample_index_source = sample_index_for_row(ratio_tag, row_index, args)

            thickness_mm, thickness_unit = thickness_to_mm(
                thickness_input, thickness_col, row_index, args.thickness_unit
            )
            mt_input = (row.get(mt_col, "") or "").strip()
            mr_input = (row.get(mr_col, "") or "").strip()
            mt = as_fraction(mt_input, "total_transmission", row_index, args.rt_unit)
            mr = as_fraction(mr_input, "total_reflection", row_index, args.rt_unit)
            if mr + mt > 1.0 + 1.0e-12:
                raise ValueError(f"row {row_index}: total_reflection + total_transmission exceeds 1")

            if args.override_g is not None:
                g_value = args.override_g
                g_source = "override --override-g"
            elif g_col is not None and (row.get(g_col, "") or "").strip():
                g_value, _ = parse_number(row[g_col], "g", row_index)
                g_source = f"csv column {g_col}"
            else:
                assert args.g is not None
                g_value = args.g
                g_source = "fallback --g"
            g_value = validate_g(g_value, row_index)

            if wavelength_col is not None and (row.get(wavelength_col, "") or "").strip():
                wavelength_nm, _ = parse_number(row[wavelength_col], "wavelength_nm", row_index)
            else:
                wavelength_nm = args.wavelength_nm

            measurements.append(
                Measurement(
                    row_index=row_index,
                    sample_id=sample_id,
                    ratio_tag=ratio_tag,
                    sample_index=sample_index,
                    sample_index_source=sample_index_source,
                    thickness_input=thickness_input,
                    thickness_unit=thickness_unit,
                    thickness_mm=thickness_mm,
                    wavelength_nm=wavelength_nm,
                    mr_input=mr_input,
                    mt_input=mt_input,
                    mr=mr,
                    mt=mt,
                    g=g_value,
                    g_source=g_source,
                )
            )
    if not measurements:
        raise ValueError(f"Input CSV has no data rows: {path}")
    return measurements


def parse_iad_result(stdout: str) -> Dict[str, str]:
    result: Optional[Dict[str, str]] = None
    for line in stdout.splitlines():
        match = RESULT_RE.match(line)
        if match:
            result = match.groupdict()
    if result is None:
        tail = "\n".join(stdout.splitlines()[-12:])
        raise ValueError(f"Could not parse IAD result row from output:\n{tail}")
    return result


def build_command(
    iad_exe: Path,
    measurement: Measurement,
    args: argparse.Namespace,
    extra_args: Sequence[str],
) -> List[str]:
    cmd = [
        str(iad_exe),
        "-r",
        f"{measurement.mr:.12g}",
        "-t",
        f"{measurement.mt:.12g}",
        "-d",
        f"{measurement.thickness_mm:.12g}",
        "-g",
        f"{measurement.g:.12g}",
        "-L",
        f"{measurement.wavelength_nm:.12g}",
        "-M",
        str(args.mc_iterations),
    ]
    slide_code = SLIDE_CONFIG[args.slide_configuration]
    cmd.extend(["-G", slide_code])
    if slide_code != "0":
        cmd.extend(["-N", f"{args.slide_index:.12g}", "-D", f"{args.slide_thickness_mm:.12g}"])
        if args.slide_od is not None:
            cmd.extend(["-E", f"{args.slide_od:.12g}"])
    cmd.extend(["-n", f"{measurement.sample_index:.12g}"])
    if args.quadrature is not None:
        cmd.extend(["-q", str(args.quadrature)])
    if args.incident_angle_deg is not None:
        cmd.extend(["-i", f"{args.incident_angle_deg:.12g}"])
    cmd.extend(extra_args)
    return cmd


def run_iad(
    iad_exe: Path,
    measurement: Measurement,
    args: argparse.Namespace,
    extra_args: Sequence[str],
) -> Dict[str, str]:
    cmd = build_command(iad_exe, measurement, args, extra_args)
    command_text = subprocess.list2cmdline(cmd)
    if args.dry_run:
        return {
            "mua_1_per_mm": "",
            "musp_1_per_mm": "",
            "mus_1_per_mm": "",
            "mua_per_um": "",
            "musp_per_um": "",
            "mus_per_um": "",
            "MR_fit": "",
            "MT_fit": "",
            "iad_status": "dry-run",
            "iad_wave": "",
            "iad_command": command_text,
            "error": "",
        }

    completed = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"IAD failed with exit code {completed.returncode}\n"
            f"command: {command_text}\n"
            f"stderr:\n{completed.stderr}\nstdout:\n{completed.stdout}"
        )
    parsed = parse_iad_result(completed.stdout)
    mua = float(parsed["mua"])
    musp = float(parsed["musp"])
    mus = musp / (1.0 - measurement.g) if measurement.g < 1.0 else math.nan
    return {
        "mua_1_per_mm": f"{mua:.12g}",
        "musp_1_per_mm": f"{musp:.12g}",
        "mus_1_per_mm": "" if not math.isfinite(mus) else f"{mus:.12g}",
        "mua_per_um": f"{mua / 1000.0:.12g}",
        "musp_per_um": f"{musp / 1000.0:.12g}",
        "mus_per_um": "" if not math.isfinite(mus) else f"{mus / 1000.0:.12g}",
        "MR_fit": parsed["mr_fit"],
        "MT_fit": parsed["mt_fit"],
        "iad_status": parsed["status"],
        "iad_wave": parsed["wave"],
        "iad_command": command_text,
        "error": "",
    }


def base_output_row(measurement: Measurement, args: argparse.Namespace) -> Dict[str, str]:
    slide_code = SLIDE_CONFIG[args.slide_configuration]
    return {
        "row_index": str(measurement.row_index),
        "sample_id": measurement.sample_id,
        "ratio_tag": measurement.ratio_tag,
        "thickness_input": measurement.thickness_input,
        "thickness_unit": measurement.thickness_unit,
        "thickness_mm": f"{measurement.thickness_mm:.12g}",
        "wavelength_nm": f"{measurement.wavelength_nm:.12g}",
        "MR_input": measurement.mr_input,
        "MT_input": measurement.mt_input,
        "MR": f"{measurement.mr:.12g}",
        "MT": f"{measurement.mt:.12g}",
        "g_input": f"{measurement.g:.12g}",
        "g_source": measurement.g_source,
        "sample_index_input": f"{measurement.sample_index:.12g}",
        "sample_index_source": measurement.sample_index_source,
        "n_model": args.n_model,
        "slide_configuration": slide_code,
        "slide_index": "" if slide_code == "0" else f"{args.slide_index:.12g}",
        "slide_thickness_mm": "" if slide_code == "0" else f"{args.slide_thickness_mm:.12g}",
        "slide_od": "" if args.slide_od is None or slide_code == "0" else f"{args.slide_od:.12g}",
        "mc_iterations": str(args.mc_iterations),
    }


def write_results(path: Path, rows: Iterable[Dict[str, str]]) -> int:
    fieldnames = [
        "row_index",
        "sample_id",
        "ratio_tag",
        "thickness_input",
        "thickness_unit",
        "thickness_mm",
        "wavelength_nm",
        "MR_input",
        "MT_input",
        "MR",
        "MT",
        "g_input",
        "g_source",
        "sample_index_input",
        "sample_index_source",
        "n_model",
        "slide_configuration",
        "slide_index",
        "slide_thickness_mm",
        "slide_od",
        "mc_iterations",
        "mua_1_per_mm",
        "musp_1_per_mm",
        "mus_1_per_mm",
        "mua_per_um",
        "musp_per_um",
        "mus_per_um",
        "MR_fit",
        "MT_fit",
        "iad_status",
        "iad_wave",
        "iad_command",
        "error",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
            count += 1
    return count


def write_effective_results(path: Path, rows: Iterable[Dict[str, str]]) -> int:
    fieldnames = [
        "sample_id",
        "ratio_tag",
        "thickness_um",
        "thickness_mm",
        "wavelength_nm",
        "MR",
        "MT",
        "sample_index_input",
        "slide_configuration",
        "slide_index",
        "slide_thickness_mm",
        "mua_eff_1_per_mm",
        "musp_eff_1_per_mm",
        "g_model",
        "MR_fit",
        "MT_fit",
        "iad_status",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out = {
                "sample_id": row.get("sample_id", ""),
                "ratio_tag": row.get("ratio_tag", ""),
                "thickness_um": row.get("thickness_input", "")
                if row.get("thickness_unit") == "um"
                else "",
                "thickness_mm": row.get("thickness_mm", ""),
                "wavelength_nm": row.get("wavelength_nm", ""),
                "MR": row.get("MR", ""),
                "MT": row.get("MT", ""),
                "sample_index_input": row.get("sample_index_input", ""),
                "slide_configuration": row.get("slide_configuration", ""),
                "slide_index": row.get("slide_index", ""),
                "slide_thickness_mm": row.get("slide_thickness_mm", ""),
                "mua_eff_1_per_mm": row.get("mua_1_per_mm", ""),
                "musp_eff_1_per_mm": row.get("musp_1_per_mm", ""),
                "g_model": row.get("g_input", ""),
                "MR_fit": row.get("MR_fit", ""),
                "MT_fit": row.get("MT_fit", ""),
                "iad_status": row.get("iad_status", ""),
            }
            writer.writerow(out)
            count += 1
    return count


def write_template(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["配比", "厚度", "总透射", "总反射", "g"])


def main() -> int:
    args = parse_args()
    if args.write_template is not None:
        template_path = resolve_path(args.write_template)
        write_template(template_path)
        print(f"template_csv,{template_path}")
        return 0

    input_csv = resolve_path(args.input_csv)
    output_file = resolve_path(args.output_file)
    iad_exe = resolve_iad_exe(args.iad_exe)
    extra_args = shlex.split(args.extra_iad_args)
    measurements = read_measurements(input_csv, args)

    rows: List[Dict[str, str]] = []
    for measurement in measurements:
        row = base_output_row(measurement, args)
        try:
            row.update(run_iad(iad_exe, measurement, args, extra_args))
        except Exception as exc:
            if args.fail_fast:
                raise
            row.update(
                {
                    "mua_1_per_mm": "",
                    "musp_1_per_mm": "",
                    "mus_1_per_mm": "",
                    "mua_per_um": "",
                    "musp_per_um": "",
                    "mus_per_um": "",
                    "MR_fit": "",
                    "MT_fit": "",
                    "iad_status": "error",
                    "iad_wave": "",
                    "iad_command": subprocess.list2cmdline(
                        build_command(iad_exe, measurement, args, extra_args)
                    ),
                    "error": str(exc).replace("\r", " ").replace("\n", " "),
                }
            )
        rows.append(row)

    count = write_results(output_file, rows)
    if args.effective_output_file is not None:
        effective_output_file = resolve_path(args.effective_output_file)
        effective_count = write_effective_results(effective_output_file, rows)
        print(f"effective_output_csv,{effective_output_file}")
        print(f"effective_rows,{effective_count}")
    errors = sum(1 for row in rows if row.get("iad_status") == "error")
    print(f"iad_exe,{iad_exe}")
    print(f"input_csv,{input_csv}")
    print(f"output_csv,{output_file}")
    print(f"rows,{count}")
    print(f"errors,{errors}")
    return 1 if errors and args.fail_fast else 0


if __name__ == "__main__":
    raise SystemExit(main())
