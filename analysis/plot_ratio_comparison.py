#!/usr/bin/env python3
"""Compare thickness trends across OpticalMC ratio batches."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from common import RUN_TAG_SEPARATOR, discover_runs, read_summary, run_filter_matches, run_tag


NUMERIC_COLUMNS = [
    "thickness_um",
    "mean_light_per_capture",
    "mean_detected_light_per_capture",
    "mean_light_per_incident",
    "mean_detected_light_per_incident",
    "detection_efficiency",
    "capture_fraction",
    "spot_rms_x",
    "spot_rms_y",
    "spot_rms_r",
    "fwhm_x",
    "fwhm_y",
    "total_source_weight",
    "total_detected_weight",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create comparison plots for different OpticalMC ratio batches."
    )
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--analysis-dir", type=Path, default=Path("analysis_results"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for comparison figures and CSV. Defaults to analysis_results/ratio_comparison.",
    )
    parser.add_argument(
        "--ratio",
        nargs="*",
        default=None,
        help="Optional ratio tags to compare. Defaults to all ratios with completed outputs.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        help="Output figure formats, e.g. --formats png pdf svg.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow writing into an existing comparison output directory.",
    )
    return parser.parse_args()


def ratio_sort_key(value: str) -> tuple[float, float, str]:
    ratio = value.split(RUN_TAG_SEPARATOR, 1)[0].split("/", 1)[0]
    parts = ratio.split("-", 1)
    if len(parts) == 2:
        try:
            return (float(parts[0]), float(parts[1]), value)
        except ValueError:
            pass
    return (math.inf, math.inf, value)


def format_ratio_label(value: str) -> str:
    return value.replace(RUN_TAG_SEPARATOR, " ").replace("-", ":")


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def unique_output_dir(base: Path, overwrite: bool) -> Path:
    if overwrite or not base.exists():
        return base
    for index in range(1, 1000):
        candidate = base.with_name(f"{base.name}_{index}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"could not choose a unique output directory near {base}")


def numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def summary_path_for_tag(outputs_dir: Path, tag: str) -> Optional[Path]:
    direct = outputs_dir / tag / "thickness_light_summary.csv"
    if direct.exists():
        return direct

    for separator in (RUN_TAG_SEPARATOR, "/"):
        if separator not in tag:
            continue
        ratio, run_label = tag.split(separator, 1)
        nested = outputs_dir / ratio / run_label / "thickness_light_summary.csv"
        if nested.exists():
            return nested
    return None


def load_ratio_summary(outputs_dir: Path, ratio: str) -> pd.DataFrame:
    summary_path = summary_path_for_tag(outputs_dir, ratio)
    if summary_path is not None:
        df = pd.read_csv(summary_path)
        df["ratio_tag"] = ratio
    else:
        runs = discover_runs(outputs_dir, [ratio])
        rows = []
        for run in runs.itertuples(index=False):
            row = read_summary(Path(run.summary_csv))
            row["ratio_tag"] = run.ratio_tag
            row["thickness_um"] = run.thickness_um
            rows.append(row)
        df = pd.DataFrame(rows)

    if df.empty:
        return df
    if "ratio_tag" not in df.columns:
        df["ratio_tag"] = ratio
    numeric(df, NUMERIC_COLUMNS)
    if "fwhm_mean" not in df.columns and {"fwhm_x", "fwhm_y"}.issubset(df.columns):
        df["fwhm_mean"] = df[["fwhm_x", "fwhm_y"]].mean(axis=1)
    if "detection_efficiency_percent" not in df.columns and "detection_efficiency" in df.columns:
        df["detection_efficiency_percent"] = 100.0 * df["detection_efficiency"]
    if "capture_fraction_percent" not in df.columns and "capture_fraction" in df.columns:
        df["capture_fraction_percent"] = 100.0 * df["capture_fraction"]
    light_col = detected_light_column(df)
    if light_col and "detected_light_norm_by_ratio" not in df.columns:
        max_value = df[light_col].max()
        df["detected_light_norm_by_ratio"] = df[light_col] / max_value if max_value > 0 else 0.0
    return df.sort_values("thickness_um").reset_index(drop=True)


def discover_ratios(outputs_dir: Path, requested: Optional[List[str]]) -> List[str]:
    if not outputs_dir.exists():
        return []
    ratio_filter = set(requested or [])
    ratios = []
    seen = set()

    def add(tag: str) -> None:
        if tag not in seen:
            seen.add(tag)
            ratios.append(tag)

    for path in outputs_dir.iterdir():
        if not path.is_dir():
            continue
        ratio = path.name
        if ratio.startswith("_") and not ratio_filter:
            continue
        if run_filter_matches(ratio_filter, ratio, None) and (path / "thickness_light_summary.csv").exists():
            add(ratio)
        for run_dir in sorted(p for p in path.iterdir() if p.is_dir()):
            if not (run_dir / "thickness_light_summary.csv").exists():
                continue
            tag = run_tag(ratio, run_dir.name)
            if run_filter_matches(ratio_filter, ratio, run_dir.name):
                add(tag)

    runs = discover_runs(outputs_dir, requested)
    if not runs.empty:
        for tag in runs["ratio_tag"].drop_duplicates():
            add(str(tag))
    return sorted(ratios, key=ratio_sort_key)


def detected_light_column(df: pd.DataFrame) -> Optional[str]:
    if "mean_detected_light_per_incident" in df.columns:
        return "mean_detected_light_per_incident"
    if "mean_detected_light_per_capture" in df.columns:
        return "mean_detected_light_per_capture"
    return None


def save_all(fig: plt.Figure, out_dir: Path, stem: str, formats: List[str], dpi: int) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    for fmt in formats:
        path = out_dir / f"{stem}.{fmt.lstrip('.')}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def style_axes(ax: plt.Axes) -> None:
    ax.grid(True, color="#d9d9d9", linewidth=0.7, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_metric(
    df: pd.DataFrame,
    y_col: str,
    y_label: str,
    title: str,
    out_dir: Path,
    stem: str,
    formats: List[str],
    dpi: int,
) -> List[Path]:
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    for ratio, group in df.groupby("ratio_tag", sort=False, observed=True):
        group = group.sort_values("thickness_um")
        ax.plot(
            group["thickness_um"],
            group[y_col],
            marker="o",
            linewidth=1.8,
            label=format_ratio_label(str(ratio)),
        )
    ax.set_xlabel("Screen thickness (um)")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    style_axes(ax)
    ax.legend(title="Ratio", frameon=False)
    return save_all(fig, out_dir, stem, formats, dpi)


def plot_summary_grid(
    df: pd.DataFrame,
    light_col: str,
    out_dir: Path,
    formats: List[str],
    dpi: int,
) -> List[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.8))
    metrics = [
        (light_col, "Detected photons / incident neutron" if light_col.endswith("incident") else "Detected photons / capture", "(a) Detected light"),
        ("detection_efficiency_percent", "Readout efficiency (%)", "(b) Readout efficiency"),
        ("capture_fraction_percent", "Capture fraction (%)", "(c) Neutron capture"),
        ("fwhm_mean", "Mean FWHM (um)", "(d) Optical spread"),
    ]
    for ax, (column, ylabel, title) in zip(axes.ravel(), metrics):
        for ratio, group in df.groupby("ratio_tag", sort=False, observed=True):
            group = group.sort_values("thickness_um")
            if column not in group.columns:
                continue
            ax.plot(
                group["thickness_um"],
                group[column],
                marker="o",
                linewidth=1.6,
                label=format_ratio_label(str(ratio)),
            )
        ax.set_xlabel("Thickness (um)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        style_axes(ax)
    axes[0, 0].legend(title="Ratio", frameon=False, fontsize=8)
    fig.suptitle("OpticalMC ratio comparison", y=0.99)
    fig.tight_layout()
    return save_all(fig, out_dir, "ratio_comparison_summary", formats, dpi)


def write_plot_data(df: pd.DataFrame, out_dir: Path, light_col: str) -> Path:
    columns = [
        "ratio_tag",
        "thickness_um",
        light_col,
        "mean_light_per_incident",
        "mean_light_per_capture",
        "mean_detected_light_per_capture",
        "detected_light_norm_by_ratio",
        "capture_fraction",
        "capture_fraction_percent",
        "detection_efficiency",
        "detection_efficiency_percent",
        "spot_rms_r",
        "fwhm_x",
        "fwhm_y",
        "fwhm_mean",
    ]
    existing = []
    for column in columns:
        if column in df.columns and column not in existing:
            existing.append(column)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "ratio_comparison_plot_data.csv"
    df[existing].to_csv(out_path, index=False)
    return out_path


def main() -> int:
    args = parse_args()
    base_output_dir = args.output_dir or (args.analysis_dir / "ratio_comparison")
    inputs_dir = args.outputs_dir.parent / "inputs"
    if is_relative_to(base_output_dir, inputs_dir):
        raise SystemExit(f"error: refusing to write comparison outputs under raw input directory: {inputs_dir}")
    out_dir = unique_output_dir(base_output_dir, args.overwrite)

    ratios = discover_ratios(args.outputs_dir, args.ratio)
    if not ratios:
        print(f"No completed OpticalMC summaries found under {args.outputs_dir}")
        return 1

    frames = []
    skipped = []
    for ratio in ratios:
        df = load_ratio_summary(args.outputs_dir, ratio)
        if df.empty:
            skipped.append(ratio)
            continue
        frames.append(df)
    if not frames:
        print(f"No non-empty ratio summaries found under {args.outputs_dir}")
        return 1

    combined = pd.concat(frames, ignore_index=True)
    combined["ratio_tag"] = pd.Categorical(combined["ratio_tag"], categories=ratios, ordered=True)
    combined = combined.sort_values(["ratio_tag", "thickness_um"]).reset_index(drop=True)

    light_col = detected_light_column(combined)
    if light_col is None:
        raise SystemExit("error: no detected light column found in summaries")

    paths: List[Path] = [write_plot_data(combined, out_dir, light_col)]
    paths += plot_metric(
        combined,
        light_col,
        "Detected photons / incident neutron" if light_col.endswith("incident") else "Detected photons / capture",
        "Detected light by ratio",
        out_dir,
        "ratio_compare_detected_light",
        args.formats,
        args.dpi,
    )
    if "detection_efficiency_percent" in combined.columns:
        paths += plot_metric(
            combined,
            "detection_efficiency_percent",
            "Readout efficiency (%)",
            "Readout efficiency by ratio",
            out_dir,
            "ratio_compare_detection_efficiency",
            args.formats,
            args.dpi,
        )
    if "capture_fraction_percent" in combined.columns:
        paths += plot_metric(
            combined,
            "capture_fraction_percent",
            "Capture fraction (%)",
            "Neutron capture by ratio",
            out_dir,
            "ratio_compare_capture_fraction",
            args.formats,
            args.dpi,
        )
    if "fwhm_mean" in combined.columns:
        paths += plot_metric(
            combined,
            "fwhm_mean",
            "Mean FWHM (um)",
            "Optical spread by ratio",
            out_dir,
            "ratio_compare_fwhm",
            args.formats,
            args.dpi,
        )
    paths += plot_summary_grid(combined, light_col, out_dir, args.formats, args.dpi)

    for path in paths:
        print(f"wrote,{path}")
    if skipped:
        print(f"skipped empty ratios: {', '.join(skipped)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
