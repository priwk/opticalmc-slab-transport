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
import numpy as np
import pandas as pd

from common import (
    READOUT_LIGHT_SCALE_DEFAULT,
    RUN_TAG_SEPARATOR,
    add_plot_metrics,
    configure_plot_style,
    detected_light_column,
    display_run_tag,
    discover_runs,
    load_phase_function,
    phase_function_polar_curve,
    phase_function_with_theta,
    read_summary,
    resolve_existing_path,
    run_filter_matches,
    run_tag,
    select_panel_thicknesses,
    style_axes,
    style_figure_axes,
)


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
        default=["png"],
        help="Output figure formats, e.g. --formats png svg.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--light-scale",
        type=float,
        default=READOUT_LIGHT_SCALE_DEFAULT,
        help="Scale for readout light per incident neutron. Defaults to 1000000.",
    )
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
    return display_run_tag(value).replace(RUN_TAG_SEPARATOR, " ").replace("-", ":")


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


def table_summary_path(analysis_dir: Path, ratio: str) -> Optional[Path]:
    candidate = analysis_dir / "tables" / f"summary_{ratio}.csv"
    if candidate.exists():
        return candidate
    return None


def load_ratio_summary(outputs_dir: Path, analysis_dir: Path, ratio: str, light_scale: float) -> pd.DataFrame:
    summary_path = summary_path_for_tag(outputs_dir, ratio)
    if summary_path is not None:
        df = pd.read_csv(summary_path)
        df["ratio_tag"] = ratio
    elif table_summary_path(analysis_dir, ratio) is not None:
        df = pd.read_csv(table_summary_path(analysis_dir, ratio))
        if "ratio_tag" not in df.columns:
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
    df = add_plot_metrics(df, light_scale)
    if "detected_light_norm_by_ratio" not in df.columns and "readout_light_per_million_incident" in df.columns:
        max_value = df["readout_light_per_million_incident"].max()
        df["detected_light_norm_by_ratio"] = (
            df["readout_light_per_million_incident"] / max_value if max_value > 0 else 0.0
        )
    return df.sort_values("thickness_um").reset_index(drop=True)


def discover_ratios(outputs_dir: Path, analysis_dir: Path, requested: Optional[List[str]]) -> List[str]:
    ratio_filter = set(requested or [])
    ratios = []
    seen = set()

    def add(tag: str) -> None:
        if tag not in seen:
            seen.add(tag)
            ratios.append(tag)

    if outputs_dir.exists():
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

    table_dir = analysis_dir / "tables"
    if table_dir.exists():
        for path in sorted(table_dir.glob("summary_*.csv")):
            if path.name in {"summary_all.csv", "summary_for_plots.csv"}:
                continue
            ratio = path.stem.removeprefix("summary_")
            if run_filter_matches(ratio_filter, ratio, None):
                add(ratio)

    runs = discover_runs(outputs_dir, requested)
    if not runs.empty:
        for tag in runs["ratio_tag"].drop_duplicates():
            add(str(tag))
    return sorted(ratios, key=ratio_sort_key)


def save_all(fig: plt.Figure, out_dir: Path, stem: str, formats: List[str], dpi: int) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    style_figure_axes(fig)
    paths: List[Path] = []
    for fmt in formats:
        path = out_dir / f"{stem}.{fmt.lstrip('.')}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths

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
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
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


def plot_metric_by_thickness_panels(
    df: pd.DataFrame,
    y_col: str,
    y_label: str,
    out_dir: Path,
    stem: str,
    formats: List[str],
    dpi: int,
) -> List[Path]:
    if y_col not in df.columns:
        return []
    panel_thicknesses = select_panel_thicknesses(df["thickness_um"], 6)
    if not panel_thicknesses:
        return []
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), squeeze=False)
    if isinstance(df["ratio_tag"].dtype, pd.CategoricalDtype):
        categories = list(df["ratio_tag"].cat.categories)
    else:
        categories = list(df["ratio_tag"].drop_duplicates())
    for ax in axes.ravel():
        ax.set_visible(False)
    for panel_index, (ax, thickness) in enumerate(zip(axes.ravel(), panel_thicknesses)):
        sub = df[np.isclose(df["thickness_um"].astype(float), thickness)].copy()
        if sub.empty:
            continue
        sub["ratio_tag"] = pd.Categorical(sub["ratio_tag"], categories=categories, ordered=True)
        sub = sub.sort_values("ratio_tag")
        x = np.arange(len(sub))
        ax.set_visible(True)
        ax.plot(x, sub[y_col], marker="o", linewidth=1.4)
        ax.set_title(f"{thickness_label_for_axis(thickness)} um", fontsize=9)
        ax.set_xticks(x)
        if panel_index // 3 == 1:
            ax.set_xticklabels(
                [format_ratio_label(str(v)) for v in sub["ratio_tag"]],
                rotation=45,
                ha="right",
                fontsize=7,
            )
        else:
            ax.set_xticklabels([])
        if panel_index % 3 == 0:
            ax.set_ylabel(y_label, fontsize=8)
        else:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)
        ax.tick_params(labelsize=7)
        style_axes(ax, square=True)
    fig.suptitle(y_label, y=0.99)
    fig.tight_layout(rect=(0.02, 0.0, 1.0, 0.96))
    return save_all(fig, out_dir, stem, formats, dpi)


def thickness_label_for_axis(value: float) -> str:
    return f"{float(value):.12g}"


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
        "readout_light_per_incident",
        "readout_light_per_million_incident",
        "readout_light_relative",
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


def load_mtf_metrics(outputs_dir: Path, analysis_dir: Path, requested: Optional[List[str]]) -> pd.DataFrame:
    table_path = analysis_dir / "tables" / "mtf_metrics.csv"
    if table_path.exists():
        df = pd.read_csv(table_path)
        if not df.empty:
            numeric(df, ["thickness_um", "mtf50_lp_per_mm", "mtf10_lp_per_mm"])
            return df

    from common import crossing_frequency, load_lsf, normalized_mtf_from_lsf

    runs = discover_runs(outputs_dir, requested)
    rows = []
    for run in runs.itertuples(index=False):
        for axis, lsf_path in [("x", Path(run.lsf_x_csv)), ("y", Path(run.lsf_y_csv))]:
            if not lsf_path.exists():
                continue
            mtf = normalized_mtf_from_lsf(load_lsf(lsf_path))
            rows.append(
                {
                    "ratio_tag": run.ratio_tag,
                    "thickness_um": run.thickness_um,
                    "axis": axis,
                    "mtf50_lp_per_mm": crossing_frequency(mtf, 0.5),
                    "mtf10_lp_per_mm": crossing_frequency(mtf, 0.1),
                }
            )
    return pd.DataFrame(rows)


def plot_mtf_comparison(
    mtf_df: pd.DataFrame, out_dir: Path, formats: List[str], dpi: int
) -> List[Path]:
    if mtf_df.empty:
        return []
    grouped = (
        mtf_df.groupby(["ratio_tag", "thickness_um"], as_index=False)[
            ["mtf50_lp_per_mm", "mtf10_lp_per_mm"]
        ]
        .mean(numeric_only=True)
        .sort_values(["ratio_tag", "thickness_um"])
    )
    paths: List[Path] = []
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    for ratio, group in grouped.groupby("ratio_tag", sort=False):
        ax.plot(group["thickness_um"], group["mtf50_lp_per_mm"], marker="o", linewidth=1.6, label=f"{format_ratio_label(str(ratio))} MTF50")
        ax.plot(group["thickness_um"], group["mtf10_lp_per_mm"], marker="s", linestyle="--", linewidth=1.4, label=f"{format_ratio_label(str(ratio))} MTF10")
    ax.set_xlabel("Thickness (um)")
    ax.set_ylabel("Spatial frequency (lp/mm)")
    ax.set_title("MTF50 and MTF10")
    ax.legend(fontsize=7, ncols=2)
    style_axes(ax, square=True)
    paths += save_all(fig, out_dir, "ratio_compare_mtf50_mtf10", formats, dpi)
    paths += plot_metric_by_thickness_panels(
        grouped,
        "mtf50_lp_per_mm",
        "MTF50 (lp/mm)",
        out_dir,
        "ratio_compare_mtf50_by_thickness_panels",
        formats,
        dpi,
    )
    paths += plot_metric_by_thickness_panels(
        grouped,
        "mtf10_lp_per_mm",
        "MTF10 (lp/mm)",
        out_dir,
        "ratio_compare_mtf10_by_thickness_panels",
        formats,
        dpi,
    )
    return paths


def plot_phase_functions(
    df: pd.DataFrame, base_dir: Path, out_dir: Path, formats: List[str], dpi: int
) -> List[Path]:
    if "phase_function_csv" not in df.columns:
        return []
    rows = []
    for ratio, group in df.groupby("ratio_tag", sort=False, observed=True):
        paths = [p for p in group["phase_function_csv"].dropna().astype(str).unique() if p.strip()]
        if not paths:
            continue
        path = resolve_existing_path(paths[0], base_dir)
        if path is None:
            continue
        try:
            phase = load_phase_function(path)
        except ValueError:
            continue
        if phase.empty:
            continue
        rows.append((ratio, phase))
    if not rows:
        return []

    paths: List[Path] = []
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    for ratio, phase in rows:
        ax.plot(phase["cos_theta"], phase["phase_function"], linewidth=1.4, label=format_ratio_label(str(ratio)))
    ax.set_xlabel("cos(theta)")
    ax.set_ylabel("p(cos(theta))")
    ax.set_title("Tabulated p(cos(theta))")
    ax.legend(fontsize=8)
    style_axes(ax, square=True)
    paths += save_all(fig, out_dir, "ratio_compare_phase_functions_mu", formats, dpi)

    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    for ratio, phase in rows:
        phase_theta = phase_function_with_theta(phase)
        ax.plot(
            phase_theta["theta_deg"],
            phase_theta["phase_function_theta"],
            linewidth=1.4,
            label=format_ratio_label(str(ratio)),
        )
    ax.set_xlabel("Scattering angle theta (deg)")
    ax.set_ylabel("p(theta) (per rad)")
    ax.set_title("Tabulated p(theta)")
    ax.legend(fontsize=8)
    style_axes(ax, square=True)
    paths += save_all(fig, out_dir, "ratio_compare_phase_functions_theta", formats, dpi)

    fig = plt.figure(figsize=(5.2, 5.2))
    ax = fig.add_subplot(111, projection="polar")
    for ratio, phase in rows:
        plot_phase_polar(ax, phase, label=format_ratio_label(str(ratio)))
    ax.set_title("Polar phase functions", pad=16)
    ax.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.18, 1.12))
    paths += save_all(fig, out_dir, "ratio_compare_phase_functions_polar", formats, dpi)

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, (ratio, phase) in zip(axes.ravel(), rows[:6]):
        ax.set_visible(True)
        ax.plot(phase["cos_theta"], phase["phase_function"], linewidth=1.4)
        ax.set_title(format_ratio_label(str(ratio)), fontsize=9)
        ax.set_xlabel("cos(theta)")
        ax.set_ylabel("p(cos(theta))")
        style_axes(ax, square=True)
    fig.suptitle("Tabulated p(cos(theta))", y=0.99)
    fig.tight_layout()
    paths += save_all(fig, out_dir, "ratio_compare_phase_function_mu_panels", formats, dpi)

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, (ratio, phase) in zip(axes.ravel(), rows[:6]):
        phase_theta = phase_function_with_theta(phase)
        ax.set_visible(True)
        ax.plot(phase_theta["theta_deg"], phase_theta["phase_function_theta"], linewidth=1.4)
        ax.set_title(format_ratio_label(str(ratio)), fontsize=9)
        ax.set_xlabel("theta (deg)")
        ax.set_ylabel("p(theta)")
        style_axes(ax, square=True)
    fig.suptitle("Tabulated p(theta)", y=0.99)
    fig.tight_layout()
    paths += save_all(fig, out_dir, "ratio_compare_phase_function_theta_panels", formats, dpi)

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), subplot_kw={"projection": "polar"}, squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, (ratio, phase) in zip(axes.ravel(), rows[:6]):
        ax.set_visible(True)
        plot_phase_polar(ax, phase)
        ax.set_title(format_ratio_label(str(ratio)), fontsize=9, pad=10)
    fig.suptitle("Polar phase functions", y=0.99)
    fig.tight_layout()
    paths += save_all(fig, out_dir, "ratio_compare_phase_function_polar_panels", formats, dpi)
    return paths


def plot_phase_polar(ax: plt.Axes, phase: pd.DataFrame, label: Optional[str] = None) -> None:
    angles, radii = phase_function_polar_curve(phase)
    ax.plot(angles, radii, linewidth=1.3, label=label)
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_thetamin(-180)
    ax.set_thetamax(180)
    ax.grid(True, linewidth=0.6, alpha=0.7)
    ax.tick_params(direction="in", labelsize=7)


def main() -> int:
    args = parse_args()
    configure_plot_style()
    base_output_dir = args.output_dir or (args.analysis_dir / "ratio_comparison")
    inputs_dir = args.outputs_dir.parent / "inputs"
    if is_relative_to(base_output_dir, inputs_dir):
        raise SystemExit(f"error: refusing to write comparison outputs under raw input directory: {inputs_dir}")
    out_dir = unique_output_dir(base_output_dir, args.overwrite)

    ratios = discover_ratios(args.outputs_dir, args.analysis_dir, args.ratio)
    if not ratios:
        print(f"No completed OpticalMC summaries found under {args.outputs_dir}")
        return 1

    frames = []
    skipped = []
    for ratio in ratios:
        df = load_ratio_summary(args.outputs_dir, args.analysis_dir, ratio, args.light_scale)
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
    if "readout_light_per_million_incident" in combined.columns:
        paths += plot_metric(
            combined,
            "readout_light_per_million_incident",
            "Readout light / 1e6 incident neutrons",
            "Relative readout light by ratio",
            out_dir,
            "ratio_compare_readout_light_per_million_incident",
            args.formats,
            args.dpi,
        )
        paths += plot_metric_by_thickness_panels(
            combined,
            "readout_light_per_million_incident",
            "Readout light / 1e6 incident neutrons",
            out_dir,
            "ratio_compare_readout_light_by_thickness_panels",
            args.formats,
            args.dpi,
        )
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
        paths += plot_metric_by_thickness_panels(
            combined,
            "fwhm_mean",
            "Mean FWHM (um)",
            out_dir,
            "ratio_compare_fwhm_by_thickness_panels",
            args.formats,
            args.dpi,
        )
    paths += plot_summary_grid(combined, light_col, out_dir, args.formats, args.dpi)
    paths += plot_mtf_comparison(load_mtf_metrics(args.outputs_dir, args.analysis_dir, args.ratio), out_dir, args.formats, args.dpi)
    paths += plot_phase_functions(combined, args.outputs_dir.parent, out_dir, args.formats, args.dpi)

    for path in paths:
        print(f"wrote,{path}")
    if skipped:
        print(f"skipped empty ratios: {', '.join(skipped)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
