#!/usr/bin/env python3
"""One-command analysis wrapper for OpticalMC outputs."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect tables and draw figures for OpticalMC runs.")
    parser.add_argument("ratio", nargs="*", help="Optional ratio tags. Empty means all ratios in outputs.")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--analysis-dir", type=Path, default=Path("analysis_results"))
    parser.add_argument("--thickness", nargs="*", default=None, help="Optional thickness list for PSF/LSF/MTF plots.")
    parser.add_argument("--include-events", action="store_true", help="Also concatenate large event summary tables.")
    parser.add_argument("--no-psf", action="store_true")
    parser.add_argument("--no-depth", action="store_true")
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    collect_cmd: List[str] = [
        sys.executable,
        str(script_dir / "collect_mc_results.py"),
        "--outputs-dir",
        str(args.outputs_dir),
        "--analysis-dir",
        str(args.analysis_dir),
    ]
    plot_cmd: List[str] = [
        sys.executable,
        str(script_dir / "plot_mc_results.py"),
        "--outputs-dir",
        str(args.outputs_dir),
        "--analysis-dir",
        str(args.analysis_dir),
        "--dpi",
        str(args.dpi),
    ]
    if args.ratio:
        collect_cmd.extend(["--ratio", *args.ratio])
        plot_cmd.extend(["--ratio", *args.ratio])
    if args.thickness:
        plot_cmd.extend(["--thickness", *args.thickness])
    if args.include_events:
        collect_cmd.append("--include-events")
    if args.no_psf:
        plot_cmd.append("--no-psf")
    if args.no_depth:
        plot_cmd.append("--no-depth")

    print(" ".join(collect_cmd), flush=True)
    subprocess.run(collect_cmd, check=True)
    print(" ".join(plot_cmd), flush=True)
    subprocess.run(plot_cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
