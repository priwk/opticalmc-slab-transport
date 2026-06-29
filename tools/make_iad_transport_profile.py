#!/usr/bin/env python3
"""Build a reduced-isotropic OpticalMC profile from IAD results."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


OUTPUT_FIELDS = [
    "ratio_tag",
    "wavelength_nm",
    "mu_a_per_um",
    "mu_s_per_um",
    "g",
    "n_eff",
    "transport_scattering_mode",
    "mu_s_prime_per_um",
    "phase_function_csv",
    "mu_s_input_per_um",
    "g_input",
    "mu_s_prime_input_per_um",
    "profile_source",
    "profile_note",
]

EFFECTIVE_MUA_COLUMN = "mua_eff_1_per_mm"
EFFECTIVE_MUSP_COLUMN = "musp_eff_1_per_mm"
DIRECT_MUA_COLUMN = "mua_per_um"
DIRECT_MUSP_COLUMN = "musp_per_um"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate IAD results into an OpticalMC optical_properties.csv for "
            "reduced-isotropic transport: g=0 and mu_s=mu_s_prime."
        )
    )
    parser.add_argument("--iad-csv", type=Path, required=True, help="Input IAD results CSV.")
    parser.add_argument("--ratio", required=True, help="Ratio tag to keep, e.g. 1-1.")
    parser.add_argument(
        "--wavelength-nm",
        type=float,
        default=450.0,
        help="Output wavelength and input wavelength filter when present. Default: 450.",
    )
    parser.add_argument(
        "--aggregate",
        choices=("mean", "median"),
        default="mean",
        help="Aggregation method for matching samples. Default: mean.",
    )
    parser.add_argument(
        "--exclude-sample",
        action="append",
        default=[],
        help="Sample id to exclude. Repeatable.",
    )
    parser.add_argument(
        "--sample-index-column",
        default="sample_index_input",
        help="Input column to average for n_eff. Default: sample_index_input.",
    )
    parser.add_argument(
        "--n-eff",
        type=float,
        default=None,
        help="Fallback n_eff if the IAD table has no sample index column.",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output optical_properties.csv.")
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=None,
        help="Optional one-row diagnostic CSV for the aggregation.",
    )
    return parser.parse_args()


def fmt(value: float) -> str:
    return f"{value:.12g}"


def read_rows(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows: List[Dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            row["_line_number"] = str(line_number)
            rows.append(row)
    if not fieldnames:
        raise ValueError(f"Input CSV has no header: {path}")
    return rows, fieldnames


def cell(row: Dict[str, str], names: Sequence[str]) -> str:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    return ""


def parse_float_cell(row: Dict[str, str], names: Sequence[str], label: str) -> float:
    value = cell(row, names)
    if value == "":
        line = row.get("_line_number", "?")
        raise ValueError(f"Missing {label} at input line {line}; tried {', '.join(names)}")
    try:
        parsed = float(value)
    except ValueError as exc:
        line = row.get("_line_number", "?")
        raise ValueError(f"Invalid {label}={value!r} at input line {line}") from exc
    if not math.isfinite(parsed):
        line = row.get("_line_number", "?")
        raise ValueError(f"Invalid non-finite {label}={value!r} at input line {line}")
    return parsed


def input_mode(fieldnames: Sequence[str]) -> str:
    fields = set(fieldnames)
    if EFFECTIVE_MUA_COLUMN in fields and EFFECTIVE_MUSP_COLUMN in fields:
        return "effective_mm"
    if DIRECT_MUA_COLUMN in fields and DIRECT_MUSP_COLUMN in fields:
        return "direct_um"
    raise ValueError(
        "Input must contain either "
        f"{EFFECTIVE_MUA_COLUMN}/{EFFECTIVE_MUSP_COLUMN} or "
        f"{DIRECT_MUA_COLUMN}/{DIRECT_MUSP_COLUMN}."
    )


def sample_id(row: Dict[str, str]) -> str:
    return cell(row, ["sample_id", "sample", "id", "sample_name", "样品", "样品编号"])


def row_ratio(row: Dict[str, str]) -> str:
    return cell(row, ["ratio_tag", "ratio", "formula", "mix", "配比", "比例"])


def row_wavelength(row: Dict[str, str]) -> str:
    return cell(row, ["wavelength_nm", "wavelength", "lambda_nm", "lambda", "波长", "波长_nm"])


def row_matches(row: Dict[str, str], ratio: str, wavelength_nm: float) -> bool:
    ratio_value = row_ratio(row)
    if ratio_value:
        if ratio_value != ratio:
            return False
    else:
        sid = sample_id(row)
        if sid and not sid.startswith(f"{ratio}-"):
            return False

    wavelength_value = row_wavelength(row)
    if wavelength_value:
        try:
            if abs(float(wavelength_value) - wavelength_nm) > 1.0e-6:
                return False
        except ValueError as exc:
            line = row.get("_line_number", "?")
            raise ValueError(f"Invalid wavelength={wavelength_value!r} at input line {line}") from exc
    return True


def aggregate(values: Sequence[float], method: str) -> float:
    if method == "mean":
        return statistics.fmean(values)
    if method == "median":
        return statistics.median(values)
    raise ValueError(f"Unsupported aggregate method: {method}")


def collect_values(
    rows: Iterable[Dict[str, str]],
    mode: str,
    args: argparse.Namespace,
) -> Tuple[List[float], List[float], List[float], List[str]]:
    mua_values: List[float] = []
    musp_values: List[float] = []
    n_values: List[float] = []
    sample_ids: List[str] = []
    excluded = {str(value).strip() for value in args.exclude_sample}

    for row in rows:
        if not row_matches(row, args.ratio, args.wavelength_nm):
            continue
        sid = sample_id(row)
        if sid in excluded:
            continue

        if mode == "effective_mm":
            mua = parse_float_cell(row, [EFFECTIVE_MUA_COLUMN], "mua_eff_1_per_mm") / 1000.0
            musp = parse_float_cell(row, [EFFECTIVE_MUSP_COLUMN], "musp_eff_1_per_mm") / 1000.0
        else:
            mua = parse_float_cell(row, [DIRECT_MUA_COLUMN], "mua_per_um")
            musp = parse_float_cell(row, [DIRECT_MUSP_COLUMN], "musp_per_um")
        mua_values.append(mua)
        musp_values.append(musp)
        sample_ids.append(sid)

        n_value = cell(row, [args.sample_index_column, "sample_index_input", "n_eff", "sample_index"])
        if n_value:
            try:
                n_parsed = float(n_value)
            except ValueError as exc:
                line = row.get("_line_number", "?")
                raise ValueError(
                    f"Invalid {args.sample_index_column}={n_value!r} at input line {line}"
                ) from exc
            if math.isfinite(n_parsed):
                n_values.append(n_parsed)

    return mua_values, musp_values, n_values, sample_ids


def write_profile(
    path: Path,
    ratio: str,
    wavelength_nm: float,
    mu_a: float,
    mu_s_prime: float,
    n_eff: float,
    profile_source: str,
    profile_note: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "ratio_tag": ratio,
                "wavelength_nm": fmt(wavelength_nm),
                "mu_a_per_um": fmt(mu_a),
                "mu_s_per_um": fmt(mu_s_prime),
                "g": "0",
                "n_eff": fmt(n_eff),
                "transport_scattering_mode": "reduced-isotropic",
                "mu_s_prime_per_um": fmt(mu_s_prime),
                "phase_function_csv": "",
                "mu_s_input_per_um": fmt(mu_s_prime),
                "g_input": "",
                "mu_s_prime_input_per_um": fmt(mu_s_prime),
                "profile_source": profile_source,
                "profile_note": profile_note,
            }
        )


def write_summary(
    path: Path,
    args: argparse.Namespace,
    sample_ids: Sequence[str],
    mua_values: Sequence[float],
    musp_values: Sequence[float],
    mu_a: float,
    mu_s_prime: float,
    n_eff: float,
    profile_source: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "ratio_tag",
        "wavelength_nm",
        "aggregate",
        "input_rows",
        "sample_ids",
        "excluded_samples",
        "mu_a_per_um",
        "mu_a_min_per_um",
        "mu_a_max_per_um",
        "mu_s_prime_per_um",
        "mu_s_prime_min_per_um",
        "mu_s_prime_max_per_um",
        "n_eff",
        "profile_source",
        "source_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "ratio_tag": args.ratio,
                "wavelength_nm": fmt(args.wavelength_nm),
                "aggregate": args.aggregate,
                "input_rows": str(len(mua_values)),
                "sample_ids": ";".join(sid for sid in sample_ids if sid),
                "excluded_samples": ";".join(args.exclude_sample),
                "mu_a_per_um": fmt(mu_a),
                "mu_a_min_per_um": fmt(min(mua_values)),
                "mu_a_max_per_um": fmt(max(mua_values)),
                "mu_s_prime_per_um": fmt(mu_s_prime),
                "mu_s_prime_min_per_um": fmt(min(musp_values)),
                "mu_s_prime_max_per_um": fmt(max(musp_values)),
                "n_eff": fmt(n_eff),
                "profile_source": profile_source,
                "source_csv": str(args.iad_csv),
            }
        )


def main() -> int:
    args = parse_args()
    try:
        rows, fieldnames = read_rows(args.iad_csv)
        mode = input_mode(fieldnames)
        mua_values, musp_values, n_values, sample_ids = collect_values(rows, mode, args)
        if not mua_values:
            raise ValueError(
                f"No rows matched ratio={args.ratio!r} and wavelength={args.wavelength_nm:.12g}."
            )
        if n_values:
            n_eff = aggregate(n_values, args.aggregate)
        elif args.n_eff is not None:
            n_eff = args.n_eff
        else:
            raise ValueError(
                f"No n_eff values found in {args.sample_index_column!r}; pass --n-eff."
            )

        mu_a = aggregate(mua_values, args.aggregate)
        mu_s_prime = aggregate(musp_values, args.aggregate)
        if mode == "effective_mm":
            profile_source = "iad_g0_effective"
            profile_note = "transport-equivalent isotropic profile from IAD g=0 effective results"
        else:
            profile_source = args.iad_csv.stem
            profile_note = "transport-equivalent isotropic profile from IAD reduced scattering results"

        write_profile(
            args.output,
            args.ratio,
            args.wavelength_nm,
            mu_a,
            mu_s_prime,
            n_eff,
            profile_source,
            profile_note,
        )
        if args.summary_output is not None:
            write_summary(
                args.summary_output,
                args,
                sample_ids,
                mua_values,
                musp_values,
                mu_a,
                mu_s_prime,
                n_eff,
                profile_source,
            )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"profile_csv,{args.output}", flush=True)
    print(f"input_rows,{len(mua_values)}", flush=True)
    print(f"mu_a_per_um,{fmt(mu_a)}", flush=True)
    print(f"mu_s_prime_per_um,{fmt(mu_s_prime)}", flush=True)
    print(f"n_eff,{fmt(n_eff)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
