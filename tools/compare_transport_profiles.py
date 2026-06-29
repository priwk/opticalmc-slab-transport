#!/usr/bin/env python3
"""Compare StageD and IAD transport profiles using mu_a and mu_s_prime only."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


STAGED_MUA_COLUMNS = [
    "recommended_mu_a_per_um",
    "mu_a_expected_mean_per_um",
    "mu_a_expected_per_um",
    "mu_a_per_um",
    "mu_a_mean_per_um",
]
STAGED_MUSP_COLUMNS = [
    "recommended_mu_s_prime_per_um",
    "mu_s_prime_mean_per_um",
    "mu_s_prime_direct_per_um",
    "mu_s_prime_per_um",
]
IAD_MUA_COLUMNS = ["mu_a_per_um"]
IAD_MUSP_COLUMNS = ["mu_s_prime_per_um", "musp_per_um"]

OUTPUT_FIELDS = [
    "ratio_tag",
    "wavelength_nm",
    "stageD_mu_a_per_um",
    "iad_mu_a_per_um",
    "stageD_mu_s_prime_per_um",
    "iad_mu_s_prime_per_um",
    "mu_a_ratio_iad_over_stageD",
    "mu_s_prime_ratio_iad_over_stageD",
    "recommended_action",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare transport-relevant optical parameters between a StageD "
            "recommended-input CSV and an IAD reduced-isotropic profile."
        )
    )
    parser.add_argument("--stageD-csv", type=Path, required=True, help="StageD recommended CSV.")
    parser.add_argument("--iad-profile", type=Path, required=True, help="IAD optical_properties.csv.")
    parser.add_argument("--output", type=Path, required=True, help="Output comparison CSV.")
    return parser.parse_args()


def fmt(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def cell(row: Dict[str, str], names: Sequence[str]) -> str:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    return ""


def first_float(row: Dict[str, str], columns: Sequence[str], label: str, path: Path) -> float:
    value = cell(row, columns)
    if value == "":
        raise ValueError(f"Missing {label} in {path}; tried {', '.join(columns)}")
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid {label}={value!r} in {path}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"Invalid non-finite {label}={value!r} in {path}")
    return parsed


def ratio_of(numerator: float, denominator: float) -> Optional[float]:
    if denominator == 0.0:
        return None
    return numerator / denominator


def row_ratio(row: Dict[str, str]) -> str:
    return cell(row, ["ratio_tag", "ratio", "配比"])


def row_wavelength(row: Dict[str, str]) -> str:
    return cell(row, ["wavelength_nm", "wavelength", "lambda_nm"])


def wavelengths_match(a: str, b: str) -> bool:
    if not a or not b:
        return True
    try:
        return abs(float(a) - float(b)) <= 1.0e-6
    except ValueError:
        return a == b


def find_stage_row(
    stage_rows: Sequence[Dict[str, str]], iad_row: Dict[str, str], stage_path: Path
) -> Dict[str, str]:
    iad_ratio = row_ratio(iad_row)
    iad_wavelength = row_wavelength(iad_row)
    for row in stage_rows:
        if row_ratio(row) == iad_ratio and wavelengths_match(row_wavelength(row), iad_wavelength):
            return row
    for row in stage_rows:
        if row_ratio(row) == iad_ratio:
            return row
    if len(stage_rows) == 1:
        return stage_rows[0]
    raise ValueError(
        f"No matching StageD row for ratio={iad_ratio!r}, wavelength={iad_wavelength!r} in {stage_path}"
    )


def compare_rows(
    stage_rows: Sequence[Dict[str, str]],
    iad_rows: Sequence[Dict[str, str]],
    stage_path: Path,
    iad_path: Path,
) -> List[Dict[str, str]]:
    output_rows: List[Dict[str, str]] = []
    for iad_row in iad_rows:
        stage_row = find_stage_row(stage_rows, iad_row, stage_path)
        ratio = row_ratio(iad_row) or row_ratio(stage_row)
        wavelength = row_wavelength(iad_row) or row_wavelength(stage_row)

        stage_mu_a = first_float(stage_row, STAGED_MUA_COLUMNS, "StageD mu_a", stage_path)
        stage_musp = first_float(stage_row, STAGED_MUSP_COLUMNS, "StageD mu_s_prime", stage_path)
        iad_mu_a = first_float(iad_row, IAD_MUA_COLUMNS, "IAD mu_a", iad_path)
        iad_musp = first_float(iad_row, IAD_MUSP_COLUMNS, "IAD mu_s_prime", iad_path)

        output_rows.append(
            {
                "ratio_tag": ratio,
                "wavelength_nm": wavelength,
                "stageD_mu_a_per_um": fmt(stage_mu_a),
                "iad_mu_a_per_um": fmt(iad_mu_a),
                "stageD_mu_s_prime_per_um": fmt(stage_musp),
                "iad_mu_s_prime_per_um": fmt(iad_musp),
                "mu_a_ratio_iad_over_stageD": fmt(ratio_of(iad_mu_a, stage_mu_a)),
                "mu_s_prime_ratio_iad_over_stageD": fmt(ratio_of(iad_musp, stage_musp)),
                "recommended_action": "run A/B optical sensitivity: stageD vs iad_g0_transport",
            }
        )
    return output_rows


def main() -> int:
    args = parse_args()
    try:
        stage_rows = read_rows(args.stageD_csv)
        iad_rows = read_rows(args.iad_profile)
        if not stage_rows:
            raise ValueError(f"No rows in {args.stageD_csv}")
        if not iad_rows:
            raise ValueError(f"No rows in {args.iad_profile}")
        rows = compare_rows(stage_rows, iad_rows, args.stageD_csv, args.iad_profile)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"comparison_csv,{args.output}", flush=True)
    print(f"rows,{len(rows)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
