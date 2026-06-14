#!/usr/bin/env python3
"""Run OpticalMC through a source-position image mask without modifying core code."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import run_opticalmc_batch as batch  # noqa: E402


def default_opticalmc_path() -> Path:
    return PROJECT_ROOT / ("OpticalMC.exe" if os.name == "nt" else "OpticalMC")


DEFAULT_STRIP_THICKNESSES = ["50", "100", "200", "500"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply a black/white square source mask to StageB neutron source positions, "
            "run OpticalMC, and write only final photon-position imaging outputs."
        )
    )
    parser.add_argument(
        "--ratio",
        nargs="*",
        default=None,
        help="BN:ZnS mass ratio tag(s), e.g. 1-1. Omit or pass all to run every ratio under inputs/stageB.",
    )
    parser.add_argument(
        "--thickness",
        nargs="*",
        default=None,
        help="Thickness list in um. Defaults to 50 100 200 500.",
    )
    parser.add_argument("--mask-image", type=Path, required=True, help="Square black/white mask image.")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--inputs-dir", type=Path, default=PROJECT_ROOT / "inputs")
    parser.add_argument("--stageb-dirname", default="stageB")
    parser.add_argument("--optical-param-dirname", default="optical_params")
    parser.add_argument("--opticalmc", type=Path, default=default_opticalmc_path())
    parser.add_argument("--preprocessor", type=Path, default=PROJECT_ROOT / "make_macro_zns_sources.py")
    parser.add_argument("--optical-properties", type=Path, default=None)
    parser.add_argument("--phase-function-csv", type=Path, default=None)
    parser.add_argument("--scattering-model", choices=("auto", "tabulated", "hg"), default="auto")
    parser.add_argument("--optical-component", choices=("bulk", "total", "boundary"), default="bulk")
    parser.add_argument("--mu-a-scale", type=float, default=1.0)
    parser.add_argument("--mu-s-scale", type=float, default=1.0)
    parser.add_argument("--transparent-optics", action="store_true")
    parser.add_argument("--yield-zns-per-MeV", type=float, default=60000.0)
    parser.add_argument("--wavelength-nm", type=float, default=450.0)
    parser.add_argument("--quench-model", choices=("none", "birks"), default="none")
    parser.add_argument("--birks-kb-um-per-keV", type=float, default=0.0)
    parser.add_argument("--readout-surface", choices=("front", "back", "both"), default="back")
    parser.add_argument(
        "--front-reflection-model",
        choices=("effective", "aluminum_fresnel"),
        default="aluminum_fresnel",
    )
    parser.add_argument("--front-reflectance", type=float, default=0.0)
    parser.add_argument(
        "--front-reflection-mode",
        choices=("none", "specular", "diffuse"),
        default="specular",
    )
    parser.add_argument("--front-aluminum-n", type=float, default=0.65)
    parser.add_argument("--front-aluminum-k", type=float, default=5.3)
    parser.add_argument("--back-reflection-model", choices=("none", "air_fresnel"), default="air_fresnel")
    parser.add_argument("--back-air-n", type=float, default=1.000293)
    parser.add_argument("--samples-per-step", type=int, default=16)
    parser.add_argument("--num-threads", type=int, default=4)
    parser.add_argument("--random-seed", type=int, default=12345)
    parser.add_argument("--incident-event-count", type=float, default=100000.0)
    parser.add_argument("--max-steps", type=int, default=10000)
    parser.add_argument("--psf-bin-size-um", type=float, default=10.0)
    parser.add_argument("--psf-range-um", type=float, default=500.0)
    parser.add_argument("--lsf-range-um", type=float, default=5000.0)
    parser.add_argument("--n-model", choices=("lorentz-lorenz", "linear"), default="lorentz-lorenz")
    parser.add_argument("--n-bn", type=float, default=1.80)
    parser.add_argument("--n-zns", type=float, default=2.36)
    parser.add_argument("--n-pmma", type=float, default=1.49)
    parser.add_argument("--n-air", type=float, default=1.000293)
    parser.add_argument("--rho-bn", type=float, default=2.10)
    parser.add_argument("--rho-zns", type=float, default=4.09)
    parser.add_argument("--bn-zns-volume", type=float, default=64.0)
    parser.add_argument("--pmma-volume", type=float, default=21.6)
    parser.add_argument("--air-volume", type=float, default=14.4)
    parser.add_argument(
        "--source-range-um",
        type=float,
        default=5000.0,
        help="Mask image covers source_x/source_y from -range to +range um.",
    )
    parser.add_argument("--mask-threshold", type=float, default=0.5)
    parser.add_argument(
        "--block-mask-edge",
        action="store_true",
        help="Block all black pixels. By default, black pixels touching white/outside are treated as usable edges.",
    )
    parser.add_argument("--outside-mask", choices=("allow", "block"), default="allow")
    parser.add_argument(
        "--image-range-um",
        type=float,
        default=5000.0,
        help="Photon image covers readout_x/readout_y from -range to +range um.",
    )
    parser.add_argument("--raw-bin-size-um", type=float, default=25.0)
    parser.add_argument("--image-pixels", type=int, default=512)
    parser.add_argument("--blur-sigma-px", type=float, default=1.2)
    parser.add_argument("--gamma", type=float, default=0.8)
    parser.add_argument("--keep-intermediates", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def thickness_label(text: str) -> str:
    return batch.thickness_label(float(text))


def discover_ratios(inputs_dir: Path, stageb_dirname: str, optical_param_dirname: str) -> List[str]:
    stageb_root = inputs_dir / stageb_dirname
    if not stageb_root.exists():
        return []
    optical_root = inputs_dir / optical_param_dirname
    ratios = []
    for path in sorted(p for p in stageb_root.iterdir() if p.is_dir()):
        if has_optical_params(optical_root, path.name):
            ratios.append(path.name)
    return ratios


def has_optical_params(optical_root: Path, ratio: str) -> bool:
    base = optical_root / ratio
    return any(
        (base / name).exists()
        for name in (
            "monte_carlo_recommended_inputs.csv",
            "monte_carlo_recommended_inputs.json",
            "rve_raw_optical_params_by_ratio.csv",
            "optical_params.csv",
        )
    )


def selected_ratios(args: argparse.Namespace) -> List[str]:
    requested = args.ratio or ["all"]
    if any(str(item).lower() == "all" for item in requested):
        ratios = discover_ratios(args.inputs_dir, args.stageb_dirname, args.optical_param_dirname)
        if not ratios:
            raise FileNotFoundError(args.inputs_dir / args.stageb_dirname)
        return ratios
    return [str(item) for item in requested]


def selected_thicknesses(args: argparse.Namespace) -> List[str]:
    values = args.thickness or DEFAULT_STRIP_THICKNESSES
    out = []
    for item in values:
        for token in str(item).split(","):
            token = token.strip()
            if token:
                out.append(thickness_label(token))
    return out


def read_mask(path: Path, threshold: float, block_edge: bool) -> np.ndarray:
    image = mpimg.imread(path)
    if image.ndim == 3:
        rgb = image[..., :3].astype(float)
        luminance = 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]
    else:
        luminance = image.astype(float)
    if luminance.max(initial=0.0) > 1.0:
        luminance = luminance / 255.0
    if luminance.shape[0] != luminance.shape[1]:
        raise ValueError(f"Mask image must be square; got {luminance.shape[1]}x{luminance.shape[0]}.")
    black = luminance < threshold
    if block_edge:
        return black

    padded = np.pad(black, 1, mode="constant", constant_values=False)
    interior = black.copy()
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            interior &= padded[1 + dy : 1 + dy + black.shape[0], 1 + dx : 1 + dx + black.shape[1]]
    return interior


def source_pixel(x_um: float, y_um: float, half_range_um: float, n: int) -> Optional[Tuple[int, int]]:
    u = (x_um + half_range_um) / (2.0 * half_range_um)
    v = (half_range_um - y_um) / (2.0 * half_range_um)
    if u < 0.0 or u > 1.0 or v < 0.0 or v > 1.0:
        return None
    col = min(n - 1, max(0, int(np.floor(u * n))))
    row = min(n - 1, max(0, int(np.floor(v * n))))
    return row, col


def float_from_row(row: Dict[str, str], column: str) -> float:
    value = row.get(column, "")
    if value is None or str(value).strip() == "":
        raise ValueError(f"missing {column}")
    return float(value)


def source_columns(fieldnames: Sequence[str]) -> Tuple[str, str, bool]:
    names = set(fieldnames)
    if {"source_x_um", "source_y_um"}.issubset(names):
        return "source_x_um", "source_y_um", False
    if {"corr_x_um", "corr_y_um"}.issubset(names):
        return "corr_x_um", "corr_y_um", True
    raise ValueError(
        "Capture anchors must contain source_x_um/source_y_um. "
        "Legacy corr_x_um/corr_y_um are accepted only as a fallback."
    )


def allowed_by_mask(
    x_um: float,
    y_um: float,
    blocked: np.ndarray,
    half_range_um: float,
    outside_mask: str,
) -> bool:
    pixel = source_pixel(x_um, y_um, half_range_um, blocked.shape[0])
    if pixel is None:
        return outside_mask == "allow"
    row, col = pixel
    return not bool(blocked[row, col])


def filter_capture_anchors(
    anchor_in: Path,
    anchor_out: Path,
    blocked: np.ndarray,
    half_range_um: float,
    outside_mask: str,
) -> Tuple[Set[str], Set[str], Dict[str, int]]:
    with anchor_in.open("r", newline="", encoding="utf-8-sig") as fin:
        reader = csv.DictReader(fin)
        fieldnames = reader.fieldnames or []
        required = {"physical_event_uid", "source_event_uid"}
        missing = sorted(required - set(fieldnames))
        if missing:
            raise ValueError(f"Capture anchors CSV is missing required columns: {', '.join(missing)}")
        x_col, y_col, legacy = source_columns(fieldnames)
        rows_by_physical: "OrderedDict[str, List[Dict[str, str]]]" = OrderedDict()
        for row in reader:
            rows_by_physical.setdefault(row["physical_event_uid"], []).append(row)

    allowed_physical: Set[str] = set()
    allowed_sources: Set[str] = set()
    counts = {
        "total_physical_events": len(rows_by_physical),
        "allowed_physical_events": 0,
        "blocked_physical_events": 0,
        "outside_source_events": 0,
        "legacy_corr_columns_used": int(legacy),
    }

    anchor_out.parent.mkdir(parents=True, exist_ok=True)
    with anchor_out.open("w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for physical_uid, rows in rows_by_physical.items():
            first = rows[0]
            x_um = float_from_row(first, x_col)
            y_um = float_from_row(first, y_col)
            if source_pixel(x_um, y_um, half_range_um, blocked.shape[0]) is None:
                counts["outside_source_events"] += 1
            if not allowed_by_mask(x_um, y_um, blocked, half_range_um, outside_mask):
                counts["blocked_physical_events"] += 1
                continue
            allowed_physical.add(physical_uid)
            counts["allowed_physical_events"] += 1
            for row in rows:
                allowed_sources.add(row["source_event_uid"])
                writer.writerow(row)
    return allowed_physical, allowed_sources, counts


def filter_csv_by_ids(
    input_path: Path,
    output_path: Path,
    allowed_physical: Set[str],
    allowed_sources: Set[str],
) -> None:
    with input_path.open("r", newline="", encoding="utf-8-sig") as fin:
        reader = csv.DictReader(fin)
        fieldnames = reader.fieldnames or []
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="", encoding="utf-8") as fout:
            writer = csv.DictWriter(fout, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            for row in reader:
                physical_uid = row.get("physical_event_uid", "")
                source_uid = row.get("source_event_uid", "")
                if physical_uid in allowed_physical or source_uid in allowed_sources:
                    writer.writerow(row)


def run_command(cmd: Sequence[object], *, cwd: Path, dry_run: bool) -> None:
    argv = [str(item) for item in cmd]
    print(" ".join(argv), flush=True)
    if dry_run:
        return
    subprocess.run(argv, cwd=str(cwd), check=True)


def optical_args(args: argparse.Namespace) -> argparse.Namespace:
    ns = argparse.Namespace(**vars(args))
    ns.source_input_mode = "stageb"
    ns.alpha_li_dirname = "alpha_li_steps"
    ns.generated_dirname = "generated_sources"
    ns.thickness = [args.thickness]
    return ns


def make_gaussian_kernel(sigma: float) -> np.ndarray:
    if sigma <= 0.0:
        return np.array([1.0])
    radius = max(1, int(np.ceil(3.0 * sigma)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-(x * x) / (2.0 * sigma * sigma))
    return kernel / kernel.sum()


def gaussian_blur(image: np.ndarray, sigma: float) -> np.ndarray:
    kernel = make_gaussian_kernel(sigma)
    tmp = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), 1, image)
    return np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), 0, tmp)


def read_photons(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    xs: List[float] = []
    ys: List[float] = []
    weights: List[float] = []
    surfaces: List[str] = []
    with path.open("r", newline="", encoding="utf-8-sig") as fin:
        reader = csv.DictReader(fin)
        for row in reader:
            xs.append(float(row["readout_x_um"]))
            ys.append(float(row["readout_y_um"]))
            weights.append(float(row.get("weight", "1") or "1"))
            surfaces.append(row.get("readout_surface", ""))
    return np.asarray(xs), np.asarray(ys), np.asarray(weights), surfaces


def write_minimal_photons(
    path: Path,
    xs: np.ndarray,
    ys: np.ndarray,
    weights: np.ndarray,
    surfaces: Sequence[str],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout, lineterminator="\n")
        writer.writerow(["readout_surface", "readout_x_um", "readout_y_um", "weight"])
        for surface, x_um, y_um, weight in zip(surfaces, xs, ys, weights):
            writer.writerow([surface, f"{x_um:.12g}", f"{y_um:.12g}", f"{weight:.12g}"])


def photon_histogram(
    xs: np.ndarray,
    ys: np.ndarray,
    weights: np.ndarray,
    bins: int,
    half_range_um: float,
) -> np.ndarray:
    hist, _, _ = np.histogram2d(
        ys,
        xs,
        bins=bins,
        range=[[-half_range_um, half_range_um], [-half_range_um, half_range_um]],
        weights=weights,
    )
    return hist


def save_hit_map(hist: np.ndarray, output_path: Path, half_range_um: float) -> None:
    fig, ax = plt.subplots(figsize=(7, 6), dpi=180)
    positive = hist[hist > 0]
    norm = LogNorm(vmin=max(positive.min(), 1.0e-12), vmax=positive.max()) if positive.size else None
    im = ax.imshow(
        hist,
        origin="lower",
        extent=(-half_range_um, half_range_um, -half_range_um, half_range_um),
        cmap="magma",
        norm=norm,
        interpolation="nearest",
    )
    ax.set_xlabel("readout x (um)")
    ax.set_ylabel("readout y (um)")
    ax.set_title("Detected photon hit map")
    fig.colorbar(im, ax=ax, label="weighted photons")
    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_radiograph(
    hist: np.ndarray,
    output_path: Path,
    blur_sigma_px: float,
    gamma: float,
) -> None:
    image = gaussian_blur(hist.astype(float), blur_sigma_px)
    if np.any(image > 0.0):
        scale = np.percentile(image[image > 0.0], 99.5)
        image = np.clip(image / max(scale, 1.0e-12), 0.0, 1.0)
    if gamma > 0.0:
        image = np.power(image, gamma)
    plt.imsave(output_path, image, cmap="gray", vmin=0.0, vmax=1.0, origin="lower")


def make_outputs(args: argparse.Namespace, detected_photons: Path, out_dir: Path) -> None:
    xs, ys, weights, surfaces = read_photons(detected_photons)
    write_minimal_photons(out_dir / "detected_photon_positions.csv", xs, ys, weights, surfaces)

    raw_bins = max(1, int(np.ceil((2.0 * args.image_range_um) / args.raw_bin_size_um)))
    raw_hist = photon_histogram(xs, ys, weights, raw_bins, args.image_range_um)
    save_hit_map(raw_hist, out_dir / "photon_hit_map.png", args.image_range_um)

    image_hist = photon_histogram(xs, ys, weights, args.image_pixels, args.image_range_um)
    save_radiograph(
        image_hist,
        out_dir / "simulated_radiograph.png",
        args.blur_sigma_px,
        args.gamma,
    )


def save_ratio_strip(
    ratio: str,
    thicknesses: Sequence[str],
    ratio_dir: Path,
    output_path: Path,
) -> Optional[Path]:
    panels: List[Tuple[str, Optional[np.ndarray]]] = []
    for thickness in thicknesses:
        path = ratio_dir / thickness_label(thickness) / "simulated_radiograph.png"
        if path.exists():
            panels.append((thickness_label(thickness), mpimg.imread(path)))
        else:
            panels.append((thickness_label(thickness), None))

    if not any(image is not None for _, image in panels):
        return None

    fig, axes = plt.subplots(1, len(panels), figsize=(3.0 * len(panels), 3.0), dpi=220)
    if len(panels) == 1:
        axes = [axes]
    for ax, (label, image) in zip(axes, panels):
        ax.set_box_aspect(1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{label} um", fontsize=10)
        if image is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes)
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.8)
            continue
        ax.imshow(image, cmap="gray", origin="lower")
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle(f"Ratio {ratio}", y=0.98, fontsize=12)
    fig.tight_layout(pad=0.4, w_pad=0.2)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    print(f"wrote,{output_path}", flush=True)
    return output_path


def run_one_mask_imaging(args: argparse.Namespace) -> int:
    label = thickness_label(args.thickness)
    out_dir = args.output_dir or (PROJECT_ROOT / "outputs" / "mask_imaging" / args.ratio / label)
    out_dir.mkdir(parents=True, exist_ok=True)

    stageb_root = args.inputs_dir / args.stageb_dirname / args.ratio
    zns_steps = stageb_root / f"{label}_zns_track_steps.csv"
    anchors = stageb_root / f"{label}_capture_anchors.csv"
    unexpected = stageb_root / f"{label}_unexpected_boundary_exits.csv"
    for path in (zns_steps, anchors, unexpected):
        if not path.exists():
            raise FileNotFoundError(path)

    temp_parent = out_dir / "intermediates" if args.keep_intermediates else None
    if temp_parent is not None:
        temp_parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(dir=temp_parent))

    try:
        blocked = read_mask(args.mask_image, args.mask_threshold, args.block_mask_edge)
        filtered_stageb = temp_dir / "stageB"
        filtered_anchors = filtered_stageb / anchors.name
        filtered_unexpected = filtered_stageb / unexpected.name
        allowed_physical, allowed_sources, counts = filter_capture_anchors(
            anchors,
            filtered_anchors,
            blocked,
            args.source_range_um,
            args.outside_mask,
        )
        filter_csv_by_ids(unexpected, filtered_unexpected, allowed_physical, allowed_sources)
        if not allowed_physical:
            raise RuntimeError("Mask blocked all source events; no OpticalMC run was launched.")

        generated_dir = temp_dir / "generated_sources"
        run_command(
            [
                sys.executable,
                args.preprocessor,
                zns_steps,
                "--output-dir",
                generated_dir,
                "--ratio-tag",
                args.ratio,
                "--thickness-um",
                label,
                "--xy-anchor-mode",
                "capture",
                "--yield-zns-per-MeV",
                args.yield_zns_per_MeV,
                "--wavelength-nm",
                args.wavelength_nm,
                "--quench-model",
                args.quench_model,
                "--birks-kb-um-per-keV",
                args.birks_kb_um_per_keV,
                "--capture-anchors",
                filtered_anchors,
                "--unexpected-boundary-exits",
                filtered_unexpected,
            ],
            cwd=PROJECT_ROOT,
            dry_run=args.dry_run,
        )

        optical_properties = batch.build_optical_properties(optical_args(args), temp_dir / "optical")
        mc_dir = temp_dir / "mc"
        config = {
            "ratio_tag": args.ratio,
            "thickness_um": float(args.thickness),
            "readout_surface": args.readout_surface,
            "front_reflection_model": args.front_reflection_model,
            "front_reflectance": args.front_reflectance,
            "front_reflection_mode": args.front_reflection_mode,
            "front_aluminum_n": args.front_aluminum_n,
            "front_aluminum_k": args.front_aluminum_k,
            "back_reflection_model": args.back_reflection_model,
            "back_air_n": args.back_air_n,
            "samples_per_step": args.samples_per_step,
            "xy_boundary": "infinite",
            "random_seed": args.random_seed,
            "incident_event_count": args.incident_event_count,
            "num_threads": args.num_threads,
            "max_steps": args.max_steps,
            "roulette_threshold": 0.0,
            "roulette_survival_probability": 1.0,
            "output_detected_photons": True,
            "psf_bin_size_um": args.psf_bin_size_um,
            "psf_range_um": args.psf_range_um,
            "lsf_range_um": args.lsf_range_um,
            "optical_properties_csv": str(optical_properties),
            "source_steps_csv": str(generated_dir / f"{label}_macro_zns_step_sources.csv"),
            "event_sources_csv": str(generated_dir / f"{label}_event_light_sources.csv"),
            "output_dir": str(mc_dir),
            "wavelength_nm": args.wavelength_nm,
        }
        config_path = temp_dir / "run_config.mask_imaging.json"
        if not args.dry_run:
            config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        run_command([args.opticalmc, config_path], cwd=PROJECT_ROOT, dry_run=args.dry_run)

        if not args.dry_run:
            detected = mc_dir / "detected_photons.csv"
            if not detected.exists():
                raise FileNotFoundError(detected)
            make_outputs(args, detected, out_dir)
        return 0
    finally:
        if args.keep_intermediates:
            print(f"kept_intermediates,{temp_dir}", flush=True)
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)


def main() -> int:
    args = parse_args()
    ratios = selected_ratios(args)
    thicknesses = selected_thicknesses(args)
    failures: List[Tuple[str, str, str]] = []

    for ratio in ratios:
        for thickness in thicknesses:
            run_args = argparse.Namespace(**vars(args))
            run_args.ratio = ratio
            run_args.thickness = thickness
            run_args.output_dir = (
                args.output_dir / ratio / thickness_label(thickness)
                if args.output_dir is not None and len(ratios) * len(thicknesses) > 1
                else args.output_dir
            )
            try:
                print(f"mask_imaging_run,{ratio},{thickness_label(thickness)}", flush=True)
                run_one_mask_imaging(run_args)
            except Exception as exc:
                failures.append((ratio, thickness_label(thickness), str(exc)))
                print(f"failed,{ratio},{thickness_label(thickness)},{exc}", flush=True)

        if not args.dry_run:
            ratio_dir = (
                args.output_dir / ratio
                if args.output_dir is not None and len(ratios) * len(thicknesses) > 1
                else PROJECT_ROOT / "outputs" / "mask_imaging" / ratio
            )
            strip_name = "radiograph_strip_" + "_".join(thickness_label(t) for t in thicknesses) + "um.png"
            save_ratio_strip(ratio, thicknesses, ratio_dir, ratio_dir / strip_name)

    if failures:
        print("mask_imaging failures:", flush=True)
        for ratio, thickness, message in failures:
            print(f"  {ratio} {thickness} um: {message}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
