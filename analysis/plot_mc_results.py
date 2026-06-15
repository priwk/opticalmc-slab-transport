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
    READOUT_LIGHT_SCALE_DEFAULT,
    add_plot_metrics,
    configure_plot_style,
    crossing_frequency,
    display_run_tag,
    discover_runs,
    ensure_dir,
    load_lsf,
    load_phase_function,
    normalized_mtf_from_lsf,
    phase_function_polar_curve,
    phase_function_with_theta,
    read_summary,
    resolve_existing_path,
    select_panel_thicknesses,
    style_axes,
    style_figure_axes,
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
    parser.add_argument(
        "--light-scale",
        type=float,
        default=READOUT_LIGHT_SCALE_DEFAULT,
        help="Scale for readout light per incident neutron. Defaults to 1000000.",
    )
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
    style_figure_axes(plt.gcf())
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
        ratio_label = display_run_tag(str(ratio))
        df = df.sort_values("thickness_um")
        x = df["thickness_um"].to_numpy(dtype=float)

        if "readout_light_per_million_incident" in df.columns:
            plt.figure(figsize=(7.2, 4.8))
            plt.plot(x, df["readout_light_per_million_incident"], "o-", linewidth=1.8)
            plt.xlabel("Thickness (um)")
            plt.ylabel("Readout light / 1e6 incident neutrons")
            plt.title(f"{ratio_label} readout light")
            savefig(fig_dir / ratio / "thickness_readout_light_per_million_incident.png", dpi)

        plt.figure(figsize=(7.2, 4.8))
        plt.plot(x, df["detection_efficiency"], "o-", label="detected / source")
        plt.xlabel("Thickness (um)")
        plt.ylabel("Detection efficiency")
        plt.title(f"{ratio_label} detection efficiency")
        savefig(fig_dir / ratio / "thickness_detection_efficiency.png", dpi)

        plt.figure(figsize=(7.2, 4.8))
        plt.plot(x, df["mean_light_per_capture"], "o-", label="source")
        plt.plot(x, df["mean_detected_light_per_capture"], "s-", label="detected")
        plt.xlabel("Thickness (um)")
        plt.ylabel("Mean photons per capture")
        plt.title(f"{ratio_label} light yield per capture")
        plt.legend()
        savefig(fig_dir / ratio / "thickness_light_per_capture.png", dpi)

        plt.figure(figsize=(5.2, 5.2))
        if "fwhm_mean" in df:
            plt.plot(x, df["fwhm_mean"], "o-", linewidth=1.8, label="Mean FWHM")
        if "fwhm_x" in df:
            plt.plot(x, df["fwhm_x"], "s--", linewidth=1.2, label="FWHM x")
        if "fwhm_y" in df:
            plt.plot(x, df["fwhm_y"], "^--", linewidth=1.2, label="FWHM y")
        plt.xlabel("Thickness (um)")
        plt.ylabel("FWHM (um)")
        plt.title(f"{ratio_label} FWHM")
        plt.legend()
        ax = plt.gca()
        style_axes(ax, square=True)
        savefig(fig_dir / ratio / "thickness_fwhm.png", dpi)

        plt.figure(figsize=(5.2, 5.2))
        plt.plot(x, df["spot_rms_r"], "o-", label="RMS r")
        if "fwhm_x" in df:
            plt.plot(x, df["fwhm_x"], "s-", label="FWHM x")
        if "fwhm_y" in df:
            plt.plot(x, df["fwhm_y"], "^-", label="FWHM y")
        plt.xlabel("Thickness (um)")
        plt.ylabel("Spot size (um)")
        plt.title(f"{ratio_label} spot spread")
        plt.legend()
        ax = plt.gca()
        style_axes(ax, square=True)
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
        plt.figure(figsize=(7.2, 4.8))
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
        plt.title(f"{ratio_label} photon budget")
        plt.legend()
        savefig(fig_dir / ratio / "thickness_photon_budget.png", dpi)


def depth_binned_event_plot(runs: pd.DataFrame, fig_dir: Path, dpi: int) -> None:
    for ratio, ratio_runs in runs.groupby("ratio_tag"):
        ratio_label = display_run_tag(str(ratio))
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
            plt.title(f"{ratio_label} depth dependence")
            plt.legend(fontsize=8, ncols=2)
            savefig(fig_dir / ratio / "event_depth_detection_efficiency.png", dpi)
        else:
            plt.close()


def plot_lsf_mtf(runs: pd.DataFrame, fig_dir: Path, table_dir: Path, dpi: int, selected: Optional[set[float]]) -> None:
    mtf_rows = []
    for run in runs.itertuples(index=False):
        if selected is not None and float(run.thickness_um) not in selected:
            continue
        ratio = run.ratio_tag
        ratio_label = display_run_tag(str(ratio))
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
            plt.title(f"{ratio_label} {label} um LSF-{axis}")
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
                plt.title(f"{ratio_label} {label} um MTF-{axis}")
                savefig(fig_dir / ratio / "mtf" / f"mtf_{axis}_{label}um.png", dpi)

    if mtf_rows:
        out = table_dir / "mtf_metrics.csv"
        mtf_df = pd.DataFrame(mtf_rows).sort_values(["ratio_tag", "thickness_um", "axis"])
        mtf_df.to_csv(out, index=False)
        print(f"wrote {out}")
        plot_mtf_threshold_trends(mtf_df, fig_dir, dpi)
        plot_lsf_thickness_panels(runs, fig_dir, dpi, selected)


def plot_mtf_threshold_trends(mtf_df: pd.DataFrame, fig_dir: Path, dpi: int) -> None:
    for ratio, df in mtf_df.groupby("ratio_tag"):
        ratio_label = display_run_tag(str(ratio))
        grouped = (
            df.groupby("thickness_um", as_index=False)[
                ["mtf10_lp_per_mm", "mtf20_lp_per_mm", "mtf50_lp_per_mm"]
            ]
            .mean(numeric_only=True)
            .sort_values("thickness_um")
        )
        if grouped.empty:
            continue
        plt.figure(figsize=(5.2, 5.2))
        for col, label, marker in [
            ("mtf50_lp_per_mm", "MTF50", "o"),
            ("mtf10_lp_per_mm", "MTF10", "s"),
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
        plt.title(f"{ratio_label} MTF thresholds vs thickness")
        plt.legend()
        ax = plt.gca()
        style_axes(ax, square=True)
        savefig(fig_dir / ratio / "thickness_mtf_thresholds.png", dpi)


def plot_lsf_thickness_panels(
    runs: pd.DataFrame, fig_dir: Path, dpi: int, selected: Optional[set[float]]
) -> None:
    for ratio, ratio_runs in runs.groupby("ratio_tag"):
        ratio_label = display_run_tag(str(ratio))
        available = ratio_runs.sort_values("thickness_um")
        if selected is not None:
            available = available[available["thickness_um"].astype(float).isin(selected)]
        panel_thicknesses = select_panel_thicknesses(available["thickness_um"], 6)
        if not panel_thicknesses:
            continue
        for axis, column in [("x", "lsf_x_csv"), ("y", "lsf_y_csv")]:
            fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), squeeze=False)
            wrote_any = False
            for ax in axes.ravel():
                ax.set_visible(False)
            for ax, thickness in zip(axes.ravel(), panel_thicknesses):
                match = available[np.isclose(available["thickness_um"].astype(float), thickness)]
                if match.empty:
                    continue
                path = Path(str(match.iloc[0][column]))
                if not path.exists():
                    continue
                lsf = load_lsf(path)
                if lsf.empty:
                    continue
                y = lsf["weight"].to_numpy(dtype=float)
                max_y = np.nanmax(y) if len(y) else np.nan
                if np.isfinite(max_y) and max_y > 0:
                    y = y / max_y
                ax.set_visible(True)
                ax.plot(lsf["position_um"], y, linewidth=1.4)
                ax.set_title(f"{thickness_label(thickness)} um", fontsize=9)
                ax.set_xlabel(f"{axis} (um)")
                ax.set_ylabel("LSF / max")
                style_axes(ax, square=True)
                wrote_any = True
            if wrote_any:
                fig.suptitle(f"{ratio_label} LSF-{axis} vs thickness", y=0.99)
                fig.tight_layout()
                path = fig_dir / ratio / "lsf" / f"lsf_{axis}_thickness_panels.png"
                path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(path, dpi=dpi)
                print(f"wrote {path}")
            plt.close(fig)


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
        plt.title(f"{display_run_tag(str(run.ratio_tag))} {thickness_label(run.thickness_um)} um PSF")
        savefig(fig_dir / run.ratio_tag / "psf" / f"psf_2d_{thickness_label(run.thickness_um)}um.png", dpi)
    plot_psf_thickness_panels(runs, fig_dir, dpi, selected)


def plot_psf_thickness_panels(
    runs: pd.DataFrame, fig_dir: Path, dpi: int, selected: Optional[set[float]]
) -> None:
    for ratio, ratio_runs in runs.groupby("ratio_tag"):
        ratio_label = display_run_tag(str(ratio))
        available = ratio_runs.sort_values("thickness_um")
        if selected is not None:
            available = available[available["thickness_um"].astype(float).isin(selected)]
        panel_thicknesses = select_panel_thicknesses(available["thickness_um"], 6)
        if not panel_thicknesses:
            continue
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), squeeze=False)
        wrote_any = False
        for ax in axes.ravel():
            ax.set_visible(False)
        for ax, thickness in zip(axes.ravel(), panel_thicknesses):
            match = available[np.isclose(available["thickness_um"].astype(float), thickness)]
            if match.empty:
                continue
            path = Path(str(match.iloc[0]["psf_2d_csv"]))
            if not path.exists():
                continue
            df = pd.read_csv(path)
            if df.empty:
                continue
            xs = np.sort(df["x_bin_center_um"].unique())
            ys = np.sort(df["y_bin_center_um"].unique())
            image = (
                df.pivot(index="y_bin_center_um", columns="x_bin_center_um", values="weight")
                .reindex(index=ys, columns=xs)
                .fillna(0)
            )
            ax.set_visible(True)
            ax.imshow(
                image.to_numpy(dtype=float),
                origin="lower",
                extent=[xs.min(), xs.max(), ys.min(), ys.max()],
                aspect="equal",
                cmap="viridis",
            )
            ax.set_title(f"{thickness_label(thickness)} um", fontsize=9)
            ax.set_xlabel("x (um)")
            ax.set_ylabel("y (um)")
            style_axes(ax, square=True)
            wrote_any = True
        if wrote_any:
            fig.suptitle(f"{ratio_label} PSF vs thickness", y=0.99)
            fig.tight_layout()
            path = fig_dir / ratio / "psf" / "psf_2d_thickness_panels.png"
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=dpi)
            print(f"wrote {path}")
        plt.close(fig)


def plot_phase_functions(summary: pd.DataFrame, fig_dir: Path, dpi: int, base_dir: Path) -> None:
    phase_rows = []
    for ratio, df in summary.groupby("ratio_tag"):
        ratio_label = display_run_tag(str(ratio))
        if "phase_function_csv" not in df.columns:
            continue
        phase_paths = [p for p in df["phase_function_csv"].dropna().astype(str).unique() if p.strip()]
        if not phase_paths:
            continue
        path = resolve_existing_path(phase_paths[0], base_dir)
        if path is None:
            continue
        try:
            phase = load_phase_function(path)
        except ValueError:
            continue
        if phase.empty:
            continue
        phase_rows.append((ratio, phase))
        phase_theta = phase_function_with_theta(phase)

        fig, ax = plt.subplots(figsize=(5.2, 5.2))
        ax.plot(phase["cos_theta"], phase["phase_function"], linewidth=1.8)
        ax.set_xlabel("cos(theta)")
        ax.set_ylabel("p(cos(theta))")
        ax.set_title(f"{ratio_label} p(cos(theta))")
        style_axes(ax, square=True)
        out = fig_dir / ratio / "phase_function_mu.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(out, dpi=dpi)
        plt.close(fig)
        print(f"wrote {out}")

        fig, ax = plt.subplots(figsize=(5.2, 5.2))
        ax.plot(phase_theta["theta_deg"], phase_theta["phase_function_theta"], linewidth=1.8)
        ax.set_xlabel("Scattering angle theta (deg)")
        ax.set_ylabel("p(theta) (per rad)")
        ax.set_title(f"{ratio_label} p(theta)")
        style_axes(ax, square=True)
        out = fig_dir / ratio / "phase_function_theta.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(out, dpi=dpi)
        plt.close(fig)
        print(f"wrote {out}")

        fig = plt.figure(figsize=(5.2, 5.2))
        ax = fig.add_subplot(111, projection="polar")
        plot_phase_polar(ax, phase)
        ax.set_title(f"{ratio_label} polar phase function", pad=16)
        out = fig_dir / ratio / "phase_function_polar.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(out, dpi=dpi)
        plt.close(fig)
        print(f"wrote {out}")

    if not phase_rows:
        return
    panel_items = phase_rows[:6]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, (ratio, phase) in zip(axes.ravel(), panel_items):
        ax.set_visible(True)
        ax.plot(phase["cos_theta"], phase["phase_function"], linewidth=1.4)
        ax.set_title(display_run_tag(str(ratio)), fontsize=9)
        ax.set_xlabel("cos(theta)")
        ax.set_ylabel("p(cos(theta))")
        style_axes(ax, square=True)
    fig.suptitle("Tabulated p(cos(theta))", y=0.99)
    fig.tight_layout()
    out = fig_dir / "phase_functions_mu_panels.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"wrote {out}")

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, (ratio, phase) in zip(axes.ravel(), panel_items):
        phase_theta = phase_function_with_theta(phase)
        ax.set_visible(True)
        ax.plot(phase_theta["theta_deg"], phase_theta["phase_function_theta"], linewidth=1.4)
        ax.set_title(display_run_tag(str(ratio)), fontsize=9)
        ax.set_xlabel("theta (deg)")
        ax.set_ylabel("p(theta)")
        style_axes(ax, square=True)
    fig.suptitle("Tabulated p(theta)", y=0.99)
    fig.tight_layout()
    out = fig_dir / "phase_functions_theta_panels.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"wrote {out}")

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8), subplot_kw={"projection": "polar"}, squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, (ratio, phase) in zip(axes.ravel(), panel_items):
        ax.set_visible(True)
        plot_phase_polar(ax, phase)
        ax.set_title(display_run_tag(str(ratio)), fontsize=9, pad=10)
    fig.suptitle("Polar phase functions", y=0.99)
    fig.tight_layout()
    out = fig_dir / "phase_functions_polar_panels.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"wrote {out}")


def plot_phase_polar(ax: plt.Axes, phase: pd.DataFrame) -> None:
    angles, radii = phase_function_polar_curve(phase)
    ax.plot(angles, radii, linewidth=1.4)
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_thetamin(-180)
    ax.set_thetamax(180)
    ax.grid(True, linewidth=0.6, alpha=0.7)
    ax.tick_params(direction="in")


def main() -> int:
    args = parse_args()
    configure_plot_style()
    runs = discover_runs(args.outputs_dir, args.ratio)
    if runs.empty:
        print(f"No completed OpticalMC runs found under {args.outputs_dir}")
        return 1
    selected = requested_thickness(args.thickness)
    fig_dir = ensure_dir(args.analysis_dir / "figures")
    table_dir = ensure_dir(args.analysis_dir / "tables")

    summary = add_plot_metrics(collect_summary(runs), args.light_scale)
    summary.to_csv(table_dir / "summary_for_plots.csv", index=False)
    plot_thickness_trends(summary, fig_dir, args.dpi)
    plot_lsf_mtf(runs, fig_dir, table_dir, args.dpi, selected)
    if not args.no_depth:
        depth_binned_event_plot(runs, fig_dir, args.dpi)
    if not args.no_psf:
        plot_psf(runs, fig_dir, args.dpi, selected)
    plot_phase_functions(summary, fig_dir, args.dpi, args.outputs_dir.parent)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
