#!/usr/bin/env python3
"""Plot thickness-level light yield and imaging summaries for OpticalMC batches."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis.common import crossing_frequency, load_lsf, normalized_mtf_from_lsf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create thickness-light-yield and paper-style summary plots."
    )
    parser.add_argument(
        "ratio",
        nargs="?",
        default=None,
        help="Ratio tag, e.g. 1-2. Defaults to reading --summary directly.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Path to thickness_light_summary.csv. Defaults to outputs/<ratio>/thickness_light_summary.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for figures. Defaults to outputs/<ratio>/figures or the summary directory.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "pdf"],
        help="Output figure formats, e.g. --formats png pdf svg.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Raster image DPI.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional figure title prefix.",
    )
    return parser.parse_args()


def numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    numeric(
        df,
        [
            "thickness_um",
            "mean_light_per_capture",
            "mean_detected_light_per_capture",
            "detection_efficiency",
            "total_source_weight",
            "total_detected_weight",
            "front_escape_weight",
            "back_escape_weight",
            "absorbed_weight",
            "lost_weight",
            "spot_rms_x",
            "spot_rms_y",
            "spot_rms_r",
            "fwhm_x",
            "fwhm_y",
            "n_events",
            "incident_event_count",
            "capture_fraction",
            "n_source_steps",
        ],
    )
    df = df.sort_values("thickness_um").reset_index(drop=True)
    if "fwhm_mean" not in df.columns and {"fwhm_x", "fwhm_y"}.issubset(df.columns):
        df["fwhm_mean"] = df[["fwhm_x", "fwhm_y"]].mean(axis=1)
    if "detection_efficiency_percent" not in df.columns and "detection_efficiency" in df.columns:
        df["detection_efficiency_percent"] = 100.0 * df["detection_efficiency"]
    light_col = (
        "mean_detected_light_per_incident"
        if "mean_detected_light_per_incident" in df.columns
        else "mean_detected_light_per_capture"
    )
    if "detected_light_norm" not in df.columns and light_col in df.columns:
        max_value = df[light_col].max()
        df["detected_light_norm"] = (
            df[light_col] / max_value if max_value > 0 else 0.0
        )
    return df


def save_all(fig: plt.Figure, out_dir: Path, stem: str, formats: List[str], dpi: int) -> List[Path]:
    paths: List[Path] = []
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = out_dir / f"{stem}.{fmt.lstrip('.')}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def thickness_label(value: float) -> str:
    return f"{value:.12g}"


def style_axes(ax: plt.Axes) -> None:
    ax.grid(True, color="#d9d9d9", linewidth=0.7, alpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def collect_mtf_metrics(df: pd.DataFrame, summary_path: Path) -> pd.DataFrame:
    run_root = summary_path.parent
    rows = []
    for row in df.itertuples(index=False):
        thickness = float(getattr(row, "thickness_um"))
        thickness_dir = run_root / thickness_label(thickness)
        for axis in ("x", "y"):
            lsf_path = thickness_dir / f"lsf_{axis}.csv"
            if not lsf_path.exists():
                continue
            lsf = load_lsf(lsf_path)
            mtf = normalized_mtf_from_lsf(lsf)
            rows.append(
                {
                    "thickness_um": thickness,
                    "axis": axis,
                    "mtf50_lp_per_mm": crossing_frequency(mtf, 0.5),
                    "mtf20_lp_per_mm": crossing_frequency(mtf, 0.2),
                    "mtf10_lp_per_mm": crossing_frequency(mtf, 0.1),
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=[
                "thickness_um",
                "axis",
                "mtf50_lp_per_mm",
                "mtf20_lp_per_mm",
                "mtf10_lp_per_mm",
            ]
        )
    return pd.DataFrame(rows).sort_values(["thickness_um", "axis"]).reset_index(drop=True)


def plot_mtf_thresholds(
    mtf_metrics: pd.DataFrame, out_dir: Path, formats: List[str], dpi: int, title: str
) -> List[Path]:
    if mtf_metrics.empty:
        return []
    grouped = (
        mtf_metrics.groupby("thickness_um", as_index=False)[
            ["mtf10_lp_per_mm", "mtf20_lp_per_mm", "mtf50_lp_per_mm"]
        ]
        .mean(numeric_only=True)
        .sort_values("thickness_um")
    )
    if grouped.empty:
        return []
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for col, label, marker in [
        ("mtf50_lp_per_mm", "MTF50", "o"),
        ("mtf20_lp_per_mm", "MTF20", "s"),
        ("mtf10_lp_per_mm", "MTF10", "^"),
    ]:
        values = grouped[col].to_numpy(dtype=float)
        if np.isfinite(values).any():
            ax.plot(
                grouped["thickness_um"],
                values,
                marker=marker,
                linewidth=1.8,
                label=label,
            )
    ax.set_xlabel("Thickness (um)")
    ax.set_ylabel("Spatial frequency (lp/mm)")
    ax.set_title(f"{title}: MTF thresholds")
    ax.legend()
    style_axes(ax)
    return save_all(fig, out_dir, "thickness_mtf_thresholds", formats, dpi)


def plot_light_curve(df: pd.DataFrame, out_dir: Path, formats: List[str], dpi: int, title: str) -> List[Path]:
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    if "mean_light_per_incident" in df.columns:
        ax.plot(
            df["thickness_um"],
            df["mean_light_per_incident"],
            marker="o",
            linewidth=1.8,
            label="Generated photons per incident neutron",
        )
    if "mean_detected_light_per_incident" in df.columns:
        ax.plot(
            df["thickness_um"],
            df["mean_detected_light_per_incident"],
            marker="s",
            linewidth=1.8,
            label="Detected photons per incident neutron",
        )
    else:
        ax.plot(
            df["thickness_um"],
            df["mean_detected_light_per_capture"],
            marker="s",
            linewidth=1.8,
            label="Detected photons per capture",
        )
    ax.set_xlabel("Screen thickness (um)")
    ax.set_ylabel("Photons")
    ax.set_title(title)
    style_axes(ax)
    ax.legend(frameon=False)
    return save_all(fig, out_dir, "thickness_light_curve", formats, dpi)


def plot_paper_summary(df: pd.DataFrame, out_dir: Path, formats: List[str], dpi: int, title: str) -> List[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.4))
    ax = axes[0, 0]
    light_col = (
        "mean_detected_light_per_incident"
        if "mean_detected_light_per_incident" in df.columns
        else "mean_detected_light_per_capture"
    )
    light_ylabel = (
        "Detected photons / incident neutron"
        if light_col == "mean_detected_light_per_incident"
        else "Detected photons / capture"
    )
    ax.plot(df["thickness_um"], df[light_col], marker="o", linewidth=1.8)
    ax.set_xlabel("Thickness (um)")
    ax.set_ylabel(light_ylabel)
    ax.set_title("(a) Light output")
    style_axes(ax)

    ax = axes[0, 1]
    ax.plot(df["thickness_um"], df["detection_efficiency_percent"], marker="o", linewidth=1.8)
    ax.set_xlabel("Thickness (um)")
    ax.set_ylabel("Detection efficiency (%)")
    ax.set_title("(b) Readout efficiency")
    style_axes(ax)

    ax = axes[1, 0]
    if "fwhm_mean" in df.columns:
        ax.plot(df["thickness_um"], df["fwhm_mean"], marker="o", linewidth=1.8, label="Mean FWHM")
    if "spot_rms_r" in df.columns:
        ax.plot(df["thickness_um"], df["spot_rms_r"], marker="s", linewidth=1.8, label="RMS radius")
    ax.set_xlabel("Thickness (um)")
    ax.set_ylabel("Spot size (um)")
    ax.set_title("(c) Optical spread")
    style_axes(ax)
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(
        df["thickness_um"],
        df["detected_light_norm"],
        marker="o",
        linewidth=1.8,
        label="Detected light / max",
    )
    if "fwhm_mean" in df.columns and df["fwhm_mean"].max() > 0:
        ax.plot(
            df["thickness_um"],
            df["fwhm_mean"] / df["fwhm_mean"].max(),
            marker="s",
            linewidth=1.8,
            label="FWHM / max",
        )
    ax.set_xlabel("Thickness (um)")
    ax.set_ylabel("Normalized value")
    ax.set_title("(d) Light-resolution trade-off")
    style_axes(ax)
    ax.legend(frameon=False)

    fig.suptitle(title, y=0.99)
    fig.tight_layout()
    return save_all(fig, out_dir, "paper_thickness_summary", formats, dpi)


def plot_efficiency_curve(df: pd.DataFrame, out_dir: Path, formats: List[str], dpi: int, title: str) -> List[Path]:
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(df["thickness_um"], df["detection_efficiency_percent"], marker="o", linewidth=1.8)
    ax.set_xlabel("Screen thickness (um)")
    ax.set_ylabel("Detection efficiency (%)")
    ax.set_title(title)
    style_axes(ax)
    return save_all(fig, out_dir, "thickness_detection_efficiency", formats, dpi)


def write_analysis_csv(df: pd.DataFrame, out_dir: Path) -> Path:
    columns = [
        "thickness_um",
        "mean_light_per_capture",
        "mean_detected_light_per_capture",
        "mean_light_per_incident",
        "mean_detected_light_per_incident",
        "detected_light_norm",
        "incident_event_count",
        "capture_fraction",
        "detection_efficiency",
        "detection_efficiency_percent",
        "spot_rms_r",
        "fwhm_x",
        "fwhm_y",
        "fwhm_mean",
    ]
    existing = [col for col in columns if col in df.columns]
    out_path = out_dir / "thickness_plot_data.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    df[existing].to_csv(out_path, index=False)
    return out_path


def main() -> int:
    args = parse_args()
    if args.summary is None:
        if not args.ratio:
            raise SystemExit("error: pass a ratio or --summary")
        args.summary = Path("outputs") / args.ratio / "thickness_light_summary.csv"
    out_dir = args.output_dir
    if out_dir is None:
        out_dir = (Path("outputs") / args.ratio / "figures") if args.ratio else args.summary.parent / "figures"

    df = load_summary(args.summary)
    if df.empty:
        raise SystemExit(f"error: no rows in {args.summary}")
    ratio_label = args.ratio or str(df.get("ratio_tag", [""])[0])
    title = args.title or f"BN/ZnS(Ag) screen, ratio {ratio_label}"

    paths: List[Path] = []
    paths += plot_light_curve(df, out_dir, args.formats, args.dpi, title)
    paths += plot_efficiency_curve(df, out_dir, args.formats, args.dpi, title)
    paths += plot_paper_summary(df, out_dir, args.formats, args.dpi, title)
    mtf_metrics = collect_mtf_metrics(df, args.summary)
    if not mtf_metrics.empty:
        mtf_path = out_dir / "thickness_mtf_metrics.csv"
        out_dir.mkdir(parents=True, exist_ok=True)
        mtf_metrics.to_csv(mtf_path, index=False)
        paths.append(mtf_path)
        paths += plot_mtf_thresholds(mtf_metrics, out_dir, args.formats, args.dpi, title)
    paths.append(write_analysis_csv(df, out_dir))

    for path in paths:
        print(f"wrote,{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
