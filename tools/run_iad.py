#!/usr/bin/env python3
"""Independent wrapper for Scott Prahl's IAD command-line program.

This script is intentionally separate from the OpticalMC batch pipeline. It reads
cleaned reflectance/transmittance CSV files from inputs/iad_inputs and writes IAD
inversion results under outputs/iad.
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
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IAD_EXE_NAME = "iad.exe" if os.name == "nt" else "iad"
DEFAULT_IAD_EXE = REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / DEFAULT_IAD_EXE_NAME
RESULT_RE = re.compile(
    r"^\s*"
    r"(?P<wave>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<mr>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<mr_fit>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<mt>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<mt_fit>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<mua>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<musp>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<g>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<status>\S)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run IAD independently on cleaned diffuse reflectance/transmittance "
            "CSV files. This does not call or modify OpticalMC."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("inputs/iad_inputs/511"),
        help="Cleaned IAD input dataset directory. Default: inputs/iad_inputs/511",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/iad/511"),
        help="Output directory for IAD results. Default: outputs/iad/511",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        default=None,
        help="Optional explicit output CSV path. Overrides the automatic filename.",
    )
    parser.add_argument(
        "--iad-exe",
        type=Path,
        default=None,
        help="Path to iad.exe/iad. Default searches third_party/iad and PATH.",
    )
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        help="Sample id or glob to process, e.g. X-X-100. Repeatable. Default: all samples.",
    )
    parser.add_argument(
        "--wavelength-nm",
        type=float,
        default=450.0,
        help="Single wavelength to process. Default: 450.0 nm.",
    )
    parser.add_argument(
        "--all-wavelengths",
        action="store_true",
        help="Process every wavelength row. This can be slow because IAD is run once per row.",
    )
    parser.add_argument(
        "--g",
        type=float,
        default=0.4,
        help="Scattering anisotropy passed to IAD with -g. Default: 0.4.",
    )
    parser.add_argument(
        "--sample-index",
        type=float,
        default=None,
        help="Optional sample refractive index passed to IAD with -n.",
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
        "--extra-iad-args",
        default="",
        help="Extra raw IAD arguments, e.g. \"-X -i 8\".",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running IAD.",
    )
    return parser.parse_args()


def resolve_iad_exe(path: Optional[Path]) -> Path:
    candidates: List[Path] = []
    if path is not None:
        candidates.append(path)
    candidates.extend(
        [
            DEFAULT_IAD_EXE,
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "src" / DEFAULT_IAD_EXE_NAME,
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "iad.exe",
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "iad",
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "src" / "iad.exe",
            REPO_ROOT / "third_party" / "iad" / "iad-3.16.3" / "src" / "iad",
            REPO_ROOT / "third_party" / "iad" / "iad.exe",
            REPO_ROOT / "third_party" / "iad" / "iad",
        ]
    )
    for candidate in candidates:
        resolved = candidate if candidate.is_absolute() else REPO_ROOT / candidate
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


def read_sample_catalog(input_dir: Path) -> List[Dict[str, str]]:
    path = input_dir / "samples.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing sample catalog: {path}\n"
            "请先准备 IAD 输入数据。例如默认 511 日期批次需要原始文件在 inputs/511，"
            "然后运行: python .\\tools\\prepare_iad_inputs.py"
        )
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def selected_samples(catalog: Sequence[Dict[str, str]], patterns: Sequence[str]) -> List[Dict[str, str]]:
    if not patterns:
        return list(catalog)
    import fnmatch

    selected: List[Dict[str, str]] = []
    for row in catalog:
        sample_id = row["sample_id"]
        if any(fnmatch.fnmatchcase(sample_id, pattern) for pattern in patterns):
            selected.append(row)
    if not selected:
        raise ValueError(f"No samples matched: {', '.join(patterns)}")
    return selected


def read_measurements(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing per-sample IAD input CSV: {path}")
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def filter_measurements(
    rows: Sequence[Dict[str, str]], all_wavelengths: bool, wavelength_nm: float
) -> List[Dict[str, str]]:
    if all_wavelengths:
        return list(rows)
    exact = [row for row in rows if math.isclose(float(row["wavelength_nm"]), wavelength_nm, abs_tol=1.0e-9)]
    if exact:
        return exact
    closest = min(rows, key=lambda row: abs(float(row["wavelength_nm"]) - wavelength_nm))
    return [closest]


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


def run_iad(
    iad_exe: Path,
    measurement: Dict[str, str],
    g: float,
    sample_index: Optional[float],
    mc_iterations: int,
    quadrature: Optional[int],
    extra_args: Sequence[str],
    dry_run: bool,
) -> Dict[str, str]:
    wavelength_nm = float(measurement["wavelength_nm"])
    thickness_um = float(measurement["thickness_um"])
    thickness_mm = thickness_um / 1000.0
    mr = float(measurement["MR"])
    mt = float(measurement["MT"])

    cmd = [
        str(iad_exe),
        "-r",
        f"{mr:.12g}",
        "-t",
        f"{mt:.12g}",
        "-d",
        f"{thickness_mm:.12g}",
        "-g",
        f"{g:.12g}",
        "-L",
        f"{wavelength_nm:.12g}",
        "-M",
        str(mc_iterations),
    ]
    if sample_index is not None:
        cmd.extend(["-n", f"{sample_index:.12g}"])
    if quadrature is not None:
        cmd.extend(["-q", str(quadrature)])
    cmd.extend(extra_args)

    if dry_run:
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
            "iad_command": subprocess.list2cmdline(cmd),
        }

    completed = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"IAD failed with exit code {completed.returncode}\n"
            f"command: {subprocess.list2cmdline(cmd)}\n"
            f"stderr:\n{completed.stderr}\nstdout:\n{completed.stdout}"
        )

    parsed = parse_iad_result(completed.stdout)
    mua = float(parsed["mua"])
    musp = float(parsed["musp"])
    mus = musp / (1.0 - g) if g < 1.0 else math.nan
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
        "iad_command": subprocess.list2cmdline(cmd),
    }


def write_results(path: Path, rows: Iterable[Dict[str, str]]) -> int:
    fieldnames = [
        "sample_id",
        "ratio_tag",
        "thickness_um",
        "thickness_mm",
        "wavelength_nm",
        "MR",
        "MT",
        "MR_corrected_percent",
        "MT_corrected_percent",
        "MR_raw_percent",
        "MT_raw_percent",
        "glass_R_percent",
        "glass_T_percent",
        "valid_iad_input",
        "correction_method",
        "g_input",
        "sample_index_input",
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
        "source_R_file",
        "source_T_file",
    ]
    count = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
            count += 1
    return count


def safe_filename_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "selection"


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir if args.input_dir.is_absolute() else REPO_ROOT / args.input_dir
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    iad_exe = resolve_iad_exe(args.iad_exe)
    extra_args = shlex.split(args.extra_iad_args)

    catalog = read_sample_catalog(input_dir)
    samples = selected_samples(catalog, args.sample)

    results: List[Dict[str, str]] = []
    for sample in samples:
        per_sample = input_dir / "per_sample" / f"{sample['sample_id']}_iad_input.csv"
        measurements = filter_measurements(
            read_measurements(per_sample),
            all_wavelengths=args.all_wavelengths,
            wavelength_nm=args.wavelength_nm,
        )
        for measurement in measurements:
            iad_result = run_iad(
                iad_exe=iad_exe,
                measurement=measurement,
                g=args.g,
                sample_index=args.sample_index,
                mc_iterations=args.mc_iterations,
                quadrature=args.quadrature,
                extra_args=extra_args,
                dry_run=args.dry_run,
            )
            thickness_um = float(measurement["thickness_um"])
            results.append(
                {
                    "sample_id": measurement["sample_id"],
                    "ratio_tag": measurement["ratio_tag"],
                    "thickness_um": measurement["thickness_um"],
                    "thickness_mm": f"{thickness_um / 1000.0:.12g}",
                    "wavelength_nm": measurement["wavelength_nm"],
                    "MR": measurement["MR"],
                    "MT": measurement["MT"],
                    "MR_corrected_percent": measurement.get("MR_corrected_percent", ""),
                    "MT_corrected_percent": measurement.get("MT_corrected_percent", ""),
                    "MR_raw_percent": measurement.get("MR_raw_percent", ""),
                    "MT_raw_percent": measurement.get("MT_raw_percent", ""),
                    "glass_R_percent": measurement.get("glass_R_percent", ""),
                    "glass_T_percent": measurement.get("glass_T_percent", ""),
                    "valid_iad_input": measurement.get("valid_iad_input", ""),
                    "correction_method": measurement.get("correction_method", ""),
                    "g_input": f"{args.g:.12g}",
                    "sample_index_input": "" if args.sample_index is None else f"{args.sample_index:.12g}",
                    "source_R_file": measurement.get("source_R_file", ""),
                    "source_T_file": measurement.get("source_T_file", ""),
                    **iad_result,
                }
            )

    wavelength_suffix = "all_wavelengths" if args.all_wavelengths else f"{args.wavelength_nm:.12g}nm"
    if args.output_file is not None:
        output_path = args.output_file if args.output_file.is_absolute() else REPO_ROOT / args.output_file
    elif args.sample:
        sample_suffix = "_".join(safe_filename_token(pattern) for pattern in args.sample)
        output_path = output_dir / f"iad_results_{sample_suffix}_{wavelength_suffix}.csv"
    else:
        output_path = output_dir / f"iad_results_{wavelength_suffix}.csv"
    count = write_results(output_path, results)
    print(f"iad_exe,{iad_exe}")
    print(f"input_dir,{input_dir}")
    print(f"output_csv,{output_path}")
    print(f"rows,{count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
