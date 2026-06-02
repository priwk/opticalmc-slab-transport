#!/usr/bin/env python3
"""Collect OpticalMC run CSVs into analysis_results tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import discover_runs, ensure_dir, read_summary, thickness_label


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect OpticalMC summary outputs.")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--analysis-dir", type=Path, default=Path("analysis_results"))
    parser.add_argument("--ratio", nargs="*", default=None, help="Optional ratio tags to collect.")
    parser.add_argument(
        "--include-events",
        action="store_true",
        help="Also concatenate event summaries. This can create large files.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    table_dir = ensure_dir(args.analysis_dir / "tables")
    runs = discover_runs(args.outputs_dir, args.ratio)
    if runs.empty:
        print(f"No completed OpticalMC runs found under {args.outputs_dir}")
        return 1

    summary_rows = []
    for run in runs.itertuples(index=False):
        row = read_summary(Path(run.summary_csv))
        row["ratio_tag"] = run.ratio_tag
        row["thickness_um"] = run.thickness_um
        row["run_dir"] = run.run_dir
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows).sort_values(["ratio_tag", "thickness_um"])
    summary_all = table_dir / "summary_all.csv"
    summary.to_csv(summary_all, index=False)
    print(f"wrote {summary_all}")

    for ratio, group in summary.groupby("ratio_tag"):
        out = table_dir / f"summary_{ratio}.csv"
        group.to_csv(out, index=False)
        print(f"wrote {out}")

    if args.include_events:
        for ratio, ratio_runs in runs.groupby("ratio_tag"):
            event_frames = []
            for run in ratio_runs.itertuples(index=False):
                path = Path(run.event_summary_csv)
                if not path.exists():
                    continue
                df = pd.read_csv(path)
                df.insert(0, "ratio_tag", run.ratio_tag)
                df.insert(1, "thickness_um", run.thickness_um)
                event_frames.append(df)
            if event_frames:
                out = table_dir / f"event_summary_{ratio}.csv"
                pd.concat(event_frames, ignore_index=True).to_csv(out, index=False)
                print(f"wrote {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
