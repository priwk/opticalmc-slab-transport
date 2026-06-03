#!/usr/bin/env python3
"""Plot OpticalMC thickness trends, event-depth trends, LSF/MTF, and PSF heatmaps."""

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
    crossing_frequency,
    discover_runs,
    ensure_dir,
    load_lsf,
    normalized_mtf_from_lsf,
    read_summary,
    thickness_label,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create OpticalMC analysis figures.")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--analysis-dir", type=Path, default=Path("analysis_results"))
    parser.add_argument("--ratio", nargs="*", default=None, help="Optional ratio tags. Defaults to all.")
    parser.add_argument(
        "--thickness",
        nargs="*",
        default=None,
        help="Optional thickness list for per-thickness PSF/LSF/MTF plots.",
    )
    parser.add_argument("--no-psf", action="store_true", help="Skip PSF heatmaps.")
    parser.add_argument("--no-depth", action="store_true", help="Skip event depth plots.")
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def requested_thickness(values: Optional[List[str]]) -> Optional[set[float]]:
    if not values:
        return None
    out = set()
    for item in values:
        for token in item.split(","):
            token = token.strip()
            if token:
                out.add(float(token))
    return out


def savefig(path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    plt.close()
    print(f"wrote {path}")


def collect_summary(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for run in runs.itertuples(index=False):
        row = read_summary(Path(run.summary_csv))
        row["ratio_tag"] = run.ratio_tag
        row["thickness_um"] = run.thickness_um
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["ratio_tag", "thickness_um"])


def plot_thickness_trends(summary: pd.DataFrame, fig_dir: Path, dpi: int) -> None:
    for ratio, df in summary.groupby("ratio_tag"):
        df = df.sort_values("thickness_um")
        x = df["thickness_um"].to_numpy(dtype=float)

        plt.figure(figsize=(7.2, 4.6))
        plt.plot(x, df["detection_efficiency"], "o-", label="detected / source")
        plt.xlabel("Thickness (um)")
        plt.ylabel("Detection efficiency")
        plt.title(f"{ratio} detection efficiency")
        plt.grid(True, alpha=0.3)
        savefig(fig_dir / ratio / "thickness_detection_efficiency.png", dpi)

        plt.figure(figsize=(7.2, 4.6))
        plt.plot(x, df["mean_light_per_capture"], "o-", label="source")
        plt.plot(x, df["mean_detected_light_per_capture"], "s-", label="detected")
        plt.xlabel("Thickness (um)")
        plt.ylabel("Mean photons per capture")
        plt.title(f"{ratio} light yield per capture")
        plt.legend()
        plt.grid(True, alpha=0.3)
        savefig(fig_dir / ratio / "thickness_light_per_capture.png", dpi)

        plt.figure(figsize=(7.2, 4.6))
        plt.plot(x, df["spot_rms_r"], "o-", label="RMS r")
        if "fwhm_x" in df:
            plt.plot(x, df["fwhm_x"], "s-", label="FWHM x")
        if "fwhm_y" in df:
            plt.plot(x, df["fwhm_y"], "^-", label="FWHM y")
        plt.xlabel("Thickness (um)")
        plt.ylabel("Spot size (um)")
        plt.title(f"{ratio} spot spread")
        plt.legend()
        plt.grid(True, alpha=0.3)
        savefig(fig_dir / ratio / "thickness_spot_spread.png", dpi)

        total = df["total_source_weight"].replace(0, np.nan)
        fractions = pd.DataFrame(
            {
                "front_escape": df["front_escape_weight"] / total,
                "back_escape": df["back_escape_weight"] / total,
                "absorbed": df["absorbed_weight"] / total,
                "lost": df["lost_weight"] / total,
            }
        )
        plt.figure(figsize=(7.2, 4.6))
        bottom = np.zeros(len(df))
        for name, color in [
            ("front_escape", "#4C78A8"),
            ("back_escape", "#72B7B2"),
            ("absorbed", "#F58518"),
            ("lost", "#B279A2"),
        ]:
            vals = fractions[name].fillna(0).to_numpy(dtype=float)
            plt.bar(x, vals, bottom=bottom, width=np.maximum(2.0, np.diff(np.r_[x, x[-1] + 20]).min() * 0.65), label=name, color=color)
            bottom += vals
        plt.xlabel("Thickness (um)")
        plt.ylabel("Weight fraction")
        plt.title(f"{ratio} photon budget")
        plt.legend()
        plt.grid(True, axis="y", alpha=0.3)
        savefig(fig_dir / ratio / "thickness_photon_budget.png", dpi)


def depth_binned_event_plot(runs: pd.DataFrame, fig_dir: Path, dpi: int) -> None:
    for ratio, ratio_runs in runs.groupby("ratio_tag"):
        plt.figure(figsize=(7.6, 4.8))
        wrote_any = False
        for run in ratio_runs.sort_values("thickness_um").itertuples(index=False):
            path = Path(run.event_summary_csv)
            if not path.exists():
                continue
            df = pd.read_csv(path)
            df = df[df["total_n_photon"] > 0].copy()
            if df.empty:
                continue
            n_bins = min(20, max(6, int(math.sqrt(len(df)))))
            bins = np.linspace(0.0, float(run.thickness_um), n_bins + 1)
            df["depth_bin"] = pd.cut(df["depth_um"], bins=bins, include_lowest=True)
            grouped = df.groupby("depth_bin", observed=False).agg(
                depth_um=("depth_um", "mean"),
                detection_efficiency=("detection_efficiency", "mean"),
            )
            plt.plot(
                grouped["depth_um"],
                grouped["detection_efficiency"],
                marker="o",
                linewidth=1.3,
                label=f"{thickness_label(run.thickness_um)} um",
            )
            wrote_any = True
        if wrote_any:
            plt.xlabel("Capture depth (um)")
            plt.ylabel("Mean event detection efficiency")
            plt.title(f"{ratio} depth dependence")
            plt.legend(fontsize=8, ncols=2)
            plt.grid(True, alpha=0.3)
            savefig(fig_dir / ratio / "event_depth_detection_efficiency.png", dpi)
        else:
            plt.close()


def plot_lsf_mtf(runs: pd.DataFrame, fig_dir: Path, table_dir: Path, dpi: int, selected: Optional[set[float]]) -> None:
    mtf_rows = []
    for run in runs.itertuples(index=False):
        if selected is not None and float(run.thickness_um) not in selected:
            continue
        ratio = run.ratio_tag
        label = thickness_label(run.thickness_um)
        for axis, lsf_path in [("x", Path(run.lsf_x_csv)), ("y", Path(run.lsf_y_csv))]:
            if not lsf_path.exists():
                continue
            lsf = load_lsf(lsf_path)
            mtf = normalized_mtf_from_lsf(lsf)
            mtf50 = crossing_frequency(mtf, 0.5)
            mtf20 = crossing_frequency(mtf, 0.2)
            mtf10 = crossing_frequency(mtf, 0.1)
            mtf_rows.append(
                {
                    "ratio_tag": ratio,
                    "thickness_um": run.thickness_um,
                    "axis": axis,
                    "mtf50_lp_per_mm": mtf50,
                    "mtf20_lp_per_mm": mtf20,
                    "mtf10_lp_per_mm": mtf10,
                }
            )

            plt.figure(figsize=(7.2, 4.4))
            plt.plot(lsf["position_um"], lsf["weight"], "-")
            plt.xlabel(f"{axis} position (um)")
            plt.ylabel("Detected weight")
            plt.title(f"{ratio} {label} um LSF-{axis}")
            plt.grid(True, alpha=0.3)
            savefig(fig_dir / ratio / "lsf" / f"lsf_{axis}_{label}um.png", dpi)

            if not mtf.empty:
                plt.figure(figsize=(7.2, 4.4))
                plt.plot(mtf["frequency_lp_per_mm"], mtf["mtf"], "-")
                plt.axhline(0.5, color="0.4", linestyle="--", linewidth=1)
                plt.axhline(0.2, color="0.5", linestyle="-.", linewidth=1)
                plt.axhline(0.1, color="0.6", linestyle=":", linewidth=1)
                plt.xlim(left=0)
                plt.ylim(0, 1.05)
                plt.xlabel("Spatial frequency (lp/mm)")
                plt.ylabel("MTF")
                plt.title(f"{ratio} {label} um MTF-{axis}")
                plt.grid(True, alpha=0.3)
                savefig(fig_dir / ratio / "mtf" / f"mtf_{axis}_{label}um.png", dpi)

    if mtf_rows:
        out = table_dir / "mtf_metrics.csv"
        mtf_df = pd.DataFrame(mtf_rows).sort_values(["ratio_tag", "thickness_um", "axis"])
        mtf_df.to_csv(out, index=False)
        print(f"wrote {out}")
        plot_mtf_threshold_trends(mtf_df, fig_dir, dpi)


def plot_mtf_threshold_trends(mtf_df: pd.DataFrame, fig_dir: Path, dpi: int) -> None:
    for ratio, df in mtf_df.groupby("ratio_tag"):
        grouped = (
            df.groupby("thickness_um", as_index=False)[
                ["mtf10_lp_per_mm", "mtf20_lp_per_mm", "mtf50_lp_per_mm"]
            ]
            .mean(numeric_only=True)
            .sort_values("thickness_um")
        )
        if grouped.empty:
            continue
        plt.figure(figsize=(7.2, 4.6))
        for col, label, marker in [
            ("mtf50_lp_per_mm", "MTF50", "o"),
            ("mtf20_lp_per_mm", "MTF20", "s"),
            ("mtf10_lp_per_mm", "MTF10", "^"),
        ]:
            values = grouped[col].to_numpy(dtype=float)
            if np.isfinite(values).any():
                plt.plot(
                    grouped["thickness_um"],
                    values,
                    marker=marker,
                    linewidth=1.8,
                    label=label,
                )
        plt.xlabel("Thickness (um)")
        plt.ylabel("Spatial frequency (lp/mm)")
        plt.title(f"{ratio} MTF thresholds vs thickness")
        plt.legend()
        plt.grid(True, alpha=0.3)
        savefig(fig_dir / ratio / "thickness_mtf_thresholds.png", dpi)


def plot_psf(runs: pd.DataFrame, fig_dir: Path, dpi: int, selected: Optional[set[float]]) -> None:
    for run in runs.itertuples(index=False):
        if selected is not None and float(run.thickness_um) not in selected:
            continue
        path = Path(run.psf_2d_csv)
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if df.empty:
            continue
        xs = np.sort(df["x_bin_center_um"].unique())
        ys = np.sort(df["y_bin_center_um"].unique())
        image = df.pivot(index="y_bin_center_um", columns="x_bin_center_um", values="weight").reindex(index=ys, columns=xs).fillna(0)
        plt.figure(figsize=(6.2, 5.4))
        plt.imshow(
            image.to_numpy(dtype=float),
            origin="lower",
            extent=[xs.min(), xs.max(), ys.min(), ys.max()],
            aspect="equal",
            cmap="viridis",
        )
        plt.colorbar(label="Detected weight")
        plt.xlabel("x (um)")
        plt.ylabel("y (um)")
        plt.title(f"{run.ratio_tag} {thickness_label(run.thickness_um)} um PSF")
        savefig(fig_dir / run.ratio_tag / "psf" / f"psf_2d_{thickness_label(run.thickness_um)}um.png", dpi)


def main() -> int:
    args = parse_args()
    runs = discover_runs(args.outputs_dir, args.ratio)
    if runs.empty:
        print(f"No completed OpticalMC runs found under {args.outputs_dir}")
        return 1
    selected = requested_thickness(args.thickness)
    fig_dir = ensure_dir(args.analysis_dir / "figures")
    table_dir = ensure_dir(args.analysis_dir / "tables")

    summary = collect_summary(runs)
    summary.to_csv(table_dir / "summary_for_plots.csv", index=False)
    plot_thickness_trends(summary, fig_dir, args.dpi)
    plot_lsf_mtf(runs, fig_dir, table_dir, args.dpi, selected)
    if not args.no_depth:
        depth_binned_event_plot(runs, fig_dir, args.dpi)
    if not args.no_psf:
        plot_psf(runs, fig_dir, args.dpi, selected)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
