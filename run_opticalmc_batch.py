#!/usr/bin/env python3
"""Batch interface for preprocessing and OpticalMC runs by BN:ZnS ratio/thickness."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def default_opticalmc_path() -> Path:
    return Path("OpticalMC.exe" if os.name == "nt" else "OpticalMC")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run make_macro_zns_sources.py and OpticalMC for one ratio and selected thicknesses."
    )
    parser.add_argument("--ratio", required=True, help="BN:ZnS mass ratio tag, e.g. 2-1, 1-1, 1-2.5.")
    parser.add_argument(
        "--thickness",
        nargs="*",
        default=None,
        help="Thickness list in um, e.g. --thickness 30 40 50 or --thickness 30,40,50. Defaults to all files for the ratio.",
    )
    parser.add_argument("--inputs-dir", type=Path, default=Path("inputs"))
    parser.add_argument("--alpha-li-dirname", default="alpha_li_steps")
    parser.add_argument("--stageb-dirname", default="stageB")
    parser.add_argument(
        "--source-input-mode",
        choices=("auto", "stageb", "alpha-li"),
        default="auto",
        help="Source preprocessing input mode. auto prefers StageB *_zns_track_steps.csv files.",
    )
    parser.add_argument("--optical-param-dirname", default="optical_params")
    parser.add_argument("--generated-dirname", default="generated_sources")
    parser.add_argument("--outputs-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--opticalmc", type=Path, default=default_opticalmc_path())
    parser.add_argument("--preprocessor", type=Path, default=Path("make_macro_zns_sources.py"))
    parser.add_argument("--optical-properties", type=Path, default=None, help="Use an existing OpticalMC optical_properties.csv instead of generating one from merged RVE params.")
    parser.add_argument(
        "--phase-function-csv",
        type=Path,
        default=None,
        help="Optional tabulated phase_function.csv. When provided, OpticalMC samples cos(theta) from this table instead of HG.",
    )
    parser.add_argument(
        "--scattering-model",
        choices=("auto", "tabulated", "hg"),
        default="auto",
        help="Scattering model for OpticalMC. auto follows recommended inputs, tabulated requires/uses a phase function, hg forces HG(g).",
    )
    parser.add_argument(
        "--transport-scattering-mode",
        choices=("reduced-isotropic", "anisotropic"),
        default="reduced-isotropic",
        help=(
            "Transport mode. reduced-isotropic uses StageD mu_s_prime as the propagated "
            "mu_s with g=0; anisotropic uses StageD mu_s, g, and optional phase functions."
        ),
    )
    parser.add_argument(
        "--run-label",
        default=None,
        help="Optional subdirectory label under outputs/<ratio>/ for this run.",
    )
    parser.add_argument(
        "--optical-component",
        choices=("bulk", "total", "boundary"),
        default="bulk",
        help="Which merged optical parameter columns to use for mu_s/g. Default bulk uses mu_s_mean_per_um and g_mean.",
    )
    parser.add_argument("--mu-a-scale", type=float, default=1.0, help="Scale mu_a after reading merged optical params; useful for sensitivity checks.")
    parser.add_argument("--mu-s-scale", type=float, default=1.0, help="Scale selected mu_s after reading merged optical params; useful for sensitivity checks.")
    parser.add_argument("--transparent-optics", action="store_true", help="Debug mode: set mu_a=mu_s=g=0 in generated optical_properties.csv.")
    parser.add_argument(
        "--xy-anchor-mode",
        choices=("capture",),
        default="capture",
        help="X/Y anchor mode for generated sources. Only capture_x/y_um is supported.",
    )
    parser.add_argument("--yield-zns-per-MeV", type=float, default=60000.0)
    parser.add_argument("--wavelength-nm", type=float, default=450.0)
    parser.add_argument("--quench-model", choices=("none", "birks"), default="none")
    parser.add_argument("--birks-kb-um-per-keV", type=float, default=0.0)
    parser.add_argument("--readout-surface", choices=("front", "back", "both"), default="back")
    parser.add_argument(
        "--front-reflection-model",
        choices=("effective", "aluminum_fresnel"),
        default="aluminum_fresnel",
        help="Front reflection probability model: effective uses front-reflectance; aluminum_fresnel computes angle-dependent Fresnel reflection from n_eff to aluminum n+ik.",
    )
    parser.add_argument(
        "--front-reflectance",
        type=float,
        default=0.0,
        help="Effective reflectance of an aluminum front plate. Used with --front-reflection-model effective.",
    )
    parser.add_argument(
        "--front-reflection-mode",
        choices=("none", "specular", "diffuse"),
        default="specular",
        help="Front boundary behavior: none keeps the old escape boundary; specular/diffuse model a reflecting aluminum plate.",
    )
    parser.add_argument("--front-aluminum-n", type=float, default=0.65, help="Real part of aluminum refractive index for aluminum_fresnel.")
    parser.add_argument("--front-aluminum-k", type=float, default=5.3, help="Extinction coefficient of aluminum for aluminum_fresnel.")
    parser.add_argument(
        "--back-reflection-model",
        choices=("none", "air_fresnel"),
        default="air_fresnel",
        help="Back boundary behavior: air_fresnel computes n_eff-to-air Fresnel reflection before counting transmitted photons as detected.",
    )
    parser.add_argument("--back-air-n", type=float, default=1.000293, help="Refractive index of air outside the back readout surface.")
    parser.add_argument("--samples-per-step", type=int, default=16)
    parser.add_argument("--num-threads", type=int, default=4)
    parser.add_argument("--random-seed", type=int, default=12345)
    parser.add_argument("--incident-event-count", type=float, default=100000.0, help="Number of incident neutron histories used to normalize per-incident light output.")
    parser.add_argument("--max-steps", type=int, default=10000)
    parser.add_argument("--psf-bin-size-um", type=float, default=2.0)
    parser.add_argument("--psf-range-um", type=float, default=500.0)
    parser.add_argument(
        "--lsf-range-um",
        type=float,
        default=5000.0,
        help="Half-width of LSF histograms used for FWHM/MTF. Keep this much larger than psf_range_um to include long tails.",
    )
    parser.add_argument("--output-detected-photons", action="store_true")
    parser.add_argument("--overwrite-sources", action="store_true", help="Regenerate source/event CSVs even if they already exist.")
    parser.add_argument(
        "--keep-generated-sources",
        action="store_true",
        help=(
            "Keep StageB-generated source/event CSVs after a full MC run. By default, "
            "newly generated StageB intermediates are removed after OpticalMC finishes."
        ),
    )
    parser.add_argument("--preprocess-only", action="store_true")
    parser.add_argument("--mc-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--n-model", choices=("lorentz-lorenz", "linear"), default="lorentz-lorenz")
    parser.add_argument("--n-bn", type=float, default=1.80)
    parser.add_argument("--n-zns", type=float, default=2.36)
    parser.add_argument("--n-pmma", type=float, default=1.49)
    parser.add_argument("--n-air", type=float, default=1.000293)
    parser.add_argument("--rho-bn", type=float, default=2.10, help="BN density in g/cm^3.")
    parser.add_argument("--rho-zns", type=float, default=4.09, help="ZnS density in g/cm^3.")
    parser.add_argument("--bn-zns-volume", type=float, default=64.0)
    parser.add_argument("--pmma-volume", type=float, default=21.6)
    parser.add_argument("--air-volume", type=float, default=14.4)
    return parser.parse_args()


def thickness_from_path(path: Path) -> Optional[float]:
    match = re.match(
        r"^([0-9]+(?:\.[0-9]+)?)_(?:alpha_li_steps|zns_track_steps)\.csv$",
        path.name,
    )
    return float(match.group(1)) if match else None


def thickness_label(value: float) -> str:
    return f"{value:.12g}"


def parse_thickness_values(values: Optional[Sequence[str]]) -> Optional[List[float]]:
    if values is None:
        return None
    out: List[float] = []
    for item in values:
        for token in item.split(","):
            token = token.strip()
            if token:
                out.append(float(token))
    return out


def ratio_numbers(ratio: str) -> Tuple[float, float]:
    if ":" in ratio:
        a, b = ratio.split(":", 1)
    elif "-" in ratio:
        a, b = ratio.split("-", 1)
    else:
        raise ValueError(f"Cannot parse ratio {ratio!r}; expected forms like 2-1 or 1:2.5")
    return float(a), float(b)


def effective_index(args: argparse.Namespace) -> float:
    bn_mass, zns_mass = ratio_numbers(args.ratio)
    v_bn_raw = bn_mass / args.rho_bn
    v_zns_raw = zns_mass / args.rho_zns
    bn_fraction_inside = v_bn_raw / (v_bn_raw + v_zns_raw)
    zns_fraction_inside = 1.0 - bn_fraction_inside

    total = args.bn_zns_volume + args.pmma_volume + args.air_volume
    f_bn_zns = args.bn_zns_volume / total
    fractions = [
        (f_bn_zns * bn_fraction_inside, args.n_bn),
        (f_bn_zns * zns_fraction_inside, args.n_zns),
        (args.pmma_volume / total, args.n_pmma),
        (args.air_volume / total, args.n_air),
    ]
    if args.n_model == "linear":
        return sum(f * n for f, n in fractions)
    ll_sum = 0.0
    for f, n in fractions:
        eps = n * n
        ll_sum += f * (eps - 1.0) / (eps + 2.0)
    eps_eff = (1.0 + 2.0 * ll_sum) / max(1.0e-15, 1.0 - ll_sum)
    return math.sqrt(eps_eff)


def find_optical_params(args: argparse.Namespace) -> Path:
    base = args.inputs_dir / args.optical_param_dirname / args.ratio
    candidates = [
        base / "monte_carlo_recommended_inputs.csv",
        base / "monte_carlo_recommended_inputs.json",
        base / "rve_raw_optical_params_by_ratio.csv",
        base / "optical_params.csv",
        args.inputs_dir / args.optical_param_dirname / "monte_carlo_recommended_inputs.csv",
        args.inputs_dir / args.optical_param_dirname / "monte_carlo_recommended_inputs.json",
        args.inputs_dir / args.optical_param_dirname / "rve_raw_optical_params_by_ratio.csv",
        Path("Output") / args.optical_param_dirname / args.ratio / "monte_carlo_recommended_inputs.csv",
        Path("Output") / args.optical_param_dirname / args.ratio / "monte_carlo_recommended_inputs.json",
        Path("Output") / args.optical_param_dirname / args.ratio / "rve_raw_optical_params_by_ratio.csv",
        Path("Output") / args.optical_param_dirname / args.ratio / "optical_params.csv",
        Path("outputs") / args.optical_param_dirname / args.ratio / "monte_carlo_recommended_inputs.csv",
        Path("outputs") / args.optical_param_dirname / args.ratio / "monte_carlo_recommended_inputs.json",
        Path("outputs") / args.optical_param_dirname / args.ratio / "rve_raw_optical_params_by_ratio.csv",
        Path("outputs") / args.optical_param_dirname / args.ratio / "optical_params.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        f"No optical params found for ratio {args.ratio}. Checked: "
        + ", ".join(str(p) for p in candidates)
    )


def read_optical_param_rows(path: Path) -> List[Dict[str, str]]:
    if path.suffix.lower() != ".json":
        with path.open(newline="", encoding="utf-8-sig") as f:
            return list(csv.DictReader(f))

    with path.open(encoding="utf-8-sig") as f:
        data = json.load(f)
    recommended = data.get("recommended_inputs", {})
    fallback = data.get("fallback_inputs", {})
    return [
        {
            "ratio": str(data.get("ratio", "")),
            "wavelength_nm": str(data.get("wavelength_nm", "")),
            "param_model": str(data.get("param_model", "")),
            "recommended_mu_a_per_um": str(recommended.get("mu_a_per_um", "")),
            "recommended_mu_s_per_um": str(recommended.get("mu_s_per_um", "")),
            "recommended_g1": str(recommended.get("g1", fallback.get("g1", ""))),
            "recommended_g2": str(recommended.get("g2", "")),
            "recommended_mu_s_prime_per_um": str(recommended.get("mu_s_prime_per_um", "")),
            "phase_function_file": str(recommended.get("phase_function_file", "")),
            "recommended_scattering_model": str(recommended.get("scattering_model", "")),
        }
    ]


def require_float(row: Dict[str, str], column: str, source: Path) -> float:
    value = row.get(column)
    if value is None or str(value).strip() == "":
        raise ValueError(f"Missing required column {column!r} in {source}")
    return float(value)


def optional_float(row: Dict[str, str], column: str) -> Optional[float]:
    value = row.get(column)
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def first_float(
    row: Dict[str, str], columns: Sequence[str], source: Path, label: str
) -> Tuple[float, str]:
    for column in columns:
        value = optional_float(row, column)
        if value is not None:
            return value, column
    raise ValueError(
        f"Missing required {label} column in {source}; tried {', '.join(columns)}"
    )


def selected_transport_columns(component: str) -> Tuple[List[str], List[str], List[str]]:
    if component == "bulk":
        return (
            ["recommended_mu_s_per_um", "mu_s_mean_per_um", "mu_s_per_um"],
            ["recommended_g1", "g_mean", "g1", "g"],
            [
                "recommended_mu_s_prime_per_um",
                "mu_s_prime_mean_per_um",
                "mu_s_prime_direct_per_um",
                "mu_s_prime_per_um",
            ],
        )
    if component == "total":
        return (
            ["mu_s_total_mean_per_um", "mu_s_mean_per_um", "mu_s_per_um"],
            ["g_total_mean", "g_mean", "g1", "g"],
            ["mu_s_prime_total_mean_per_um", "mu_s_prime_mean_per_um", "mu_s_prime_direct_per_um"],
        )
    if component == "boundary":
        return (
            ["mu_s_boundary_mean_per_um"],
            ["g_boundary_mean"],
            ["mu_s_prime_boundary_mean_per_um"],
        )
    raise ValueError(f"Unsupported optical component: {component}")


def phase_function_from_row(
    args: argparse.Namespace, row: Dict[str, str], rve_path: Path
) -> str:
    if getattr(args, "transport_scattering_mode", "reduced-isotropic") == "reduced-isotropic":
        return ""
    if args.scattering_model == "hg":
        return ""
    candidate = args.phase_function_csv
    if candidate is None:
        recommended_model = (row.get("recommended_scattering_model") or "").strip()
        use_recommended_phase = args.scattering_model == "tabulated" or (
            args.scattering_model == "auto"
            and recommended_model in ("", "tabulated_phase_function", "tabulated")
        )
        if not use_recommended_phase:
            return ""
        for column in ("phase_function_csv", "phase_function_path", "phase_function_file"):
            value = row.get(column)
            if value and value.strip():
                candidate = Path(value.strip())
                break
        if candidate is None and (rve_path.parent / "phase_function_mean_by_ratio.csv").exists():
            candidate = Path("phase_function_mean_by_ratio.csv")
    if candidate is None:
        if args.scattering_model == "tabulated":
            raise ValueError(
                "Scattering model is tabulated but no phase function was found. "
                "Provide --phase-function-csv or phase_function_file in recommended inputs."
            )
        return ""
    if candidate.is_absolute():
        return str(candidate)
    if candidate.exists():
        return str(candidate)
    beside_rve = rve_path.parent / candidate
    if beside_rve.exists():
        return str(beside_rve)
    return str(candidate)


def run_root_for_args(args: argparse.Namespace) -> Path:
    if args.run_label:
        return args.outputs_dir / args.ratio / args.run_label
    if args.transport_scattering_mode == "anisotropic":
        return args.outputs_dir / args.ratio / "modeB_anisotropic"
    return args.outputs_dir / args.ratio


def build_optical_properties(args: argparse.Namespace, run_root: Path) -> Path:
    if args.optical_properties is not None:
        return args.optical_properties

    if not hasattr(args, "transport_scattering_mode"):
        args.transport_scattering_mode = "reduced-isotropic"
    if not hasattr(args, "run_label"):
        args.run_label = None

    rve_path = find_optical_params(args)
    n_eff = effective_index(args)
    out_path = run_root / "optical_properties.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mu_s_columns, g_columns, mu_s_prime_columns = selected_transport_columns(args.optical_component)
    mu_a_expected_columns = [
        "recommended_mu_a_per_um",
        "mu_a_expected_mean_per_um",
        "mu_a_expected_per_um",
    ]
    mu_a_fallback_columns = ["mu_a_mean_per_um", "mu_a_count_per_um", "mu_a_per_um"]

    rows_written = 0
    warned_mu_a_fallback = False
    with out_path.open("w", newline="", encoding="utf-8") as fout:
        reader = read_optical_param_rows(rve_path)
        writer = csv.DictWriter(
            fout,
            fieldnames=[
                "ratio_tag",
                "wavelength_nm",
                "mu_a_per_um",
                "mu_s_per_um",
                "g",
                "n_eff",
                "transport_scattering_mode",
                "mu_s_prime_per_um",
                "phase_function_csv",
                "mu_s_input_per_um",
                "g_input",
                "mu_s_prime_input_per_um",
                "mu_a_source_column",
                "mu_s_source_column",
                "g_source_column",
                "mu_s_prime_source_column",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            row_ratio = row.get("ratio") or row.get("ratio_tag") or args.ratio
            if row_ratio != args.ratio:
                continue
            if args.transparent_optics:
                mu_a = 0.0
                mu_s = 0.0
                g = 0.0
                mu_s_prime = 0.0
                mu_s_input = 0.0
                g_input = 0.0
                mu_s_prime_input = 0.0
                mu_a_source = "transparent_optics"
                mu_s_source = "transparent_optics"
                g_source = "transparent_optics"
                mu_s_prime_source = "transparent_optics"
            else:
                try:
                    mu_a_raw, mu_a_source = first_float(
                        row, mu_a_expected_columns, rve_path, "expected absorption"
                    )
                except ValueError:
                    mu_a_raw, mu_a_source = first_float(
                        row, mu_a_fallback_columns, rve_path, "absorption"
                    )
                    if not warned_mu_a_fallback:
                        print(
                            "warning: expected absorption column not found; "
                            f"using {mu_a_source} from {rve_path}. "
                            "Prefer recommended_mu_a_per_um, mu_a_expected_mean_per_um, or "
                            "mu_a_expected_per_um for new GO_RVE inputs.",
                            file=sys.stderr,
                        )
                        warned_mu_a_fallback = True
                mu_s_raw, mu_s_source = first_float(row, mu_s_columns, rve_path, "scattering")
                g, g_source = first_float(row, g_columns, rve_path, "anisotropy")
                mu_s_prime_raw, mu_s_prime_source = first_float(
                    row, mu_s_prime_columns, rve_path, "reduced scattering"
                )
                mu_s_input = mu_s_raw
                g_input = g
                mu_s_prime_input = mu_s_prime_raw
                mu_a = mu_a_raw * args.mu_a_scale
                if args.transport_scattering_mode == "reduced-isotropic":
                    mu_s = mu_s_prime_raw * args.mu_s_scale
                    g = 0.0
                    mu_s_prime = mu_s_prime_raw * args.mu_s_scale
                    mu_s_source = mu_s_prime_source
                    g_source = "reduced_isotropic"
                else:
                    mu_s = mu_s_raw * args.mu_s_scale
                    mu_s_prime = mu_s_prime_raw * args.mu_s_scale
                    mu_s_prime_from_g = mu_s * (1.0 - g)
                    denom = max(abs(mu_s_prime), abs(mu_s_prime_from_g), 1.0e-30)
                    if abs(mu_s_prime - mu_s_prime_from_g) / denom > 0.10:
                        print(
                            "warning: anisotropic transport has mu_s_prime_per_um "
                            f"({mu_s_prime:.12g}) differing from mu_s_per_um*(1-g) "
                            f"({mu_s_prime_from_g:.12g}) by more than 10%.",
                            file=sys.stderr,
                        )
            phase_function_csv = "" if args.transparent_optics else phase_function_from_row(args, row, rve_path)
            writer.writerow(
                {
                    "ratio_tag": args.ratio,
                    "wavelength_nm": row.get("wavelength_nm", args.wavelength_nm),
                    "mu_a_per_um": f"{mu_a:.12g}",
                    "mu_s_per_um": f"{mu_s:.12g}",
                    "g": f"{g:.12g}",
                    "n_eff": f"{n_eff:.12g}",
                    "transport_scattering_mode": args.transport_scattering_mode,
                    "mu_s_prime_per_um": f"{mu_s_prime:.12g}",
                    "phase_function_csv": phase_function_csv,
                    "mu_s_input_per_um": f"{mu_s_input:.12g}",
                    "g_input": f"{g_input:.12g}",
                    "mu_s_prime_input_per_um": f"{mu_s_prime_input:.12g}",
                    "mu_a_source_column": mu_a_source,
                    "mu_s_source_column": mu_s_source,
                    "g_source_column": g_source,
                    "mu_s_prime_source_column": mu_s_prime_source,
                }
            )
            rows_written += 1
    if rows_written == 0:
        raise ValueError(f"No row for ratio {args.ratio} in {rve_path}")
    return out_path


def discover_thicknesses(args: argparse.Namespace) -> List[float]:
    requested = parse_thickness_values(args.thickness)
    if requested is not None and requested:
        return requested

    values = []
    checked: List[Path] = []
    if args.source_input_mode in ("auto", "stageb"):
        ratio_dir = args.inputs_dir / args.stageb_dirname / args.ratio
        checked.append(ratio_dir)
        for path in ratio_dir.glob("*_zns_track_steps.csv"):
            t = thickness_from_path(path)
            if t is not None:
                values.append(t)
        if values or args.source_input_mode == "stageb":
            if not values:
                raise FileNotFoundError(f"No StageB ZnS track files found in {ratio_dir}")
            return sorted(values)

    if args.source_input_mode in ("auto", "alpha-li"):
        ratio_dir = args.inputs_dir / args.alpha_li_dirname / args.ratio
        checked.append(ratio_dir)
        for path in ratio_dir.glob("*_alpha_li_steps.csv"):
            t = thickness_from_path(path)
            if t is not None:
                values.append(t)
    if not values:
        raise FileNotFoundError(
            "No source input files found. Checked: " + ", ".join(str(path) for path in checked)
        )
    return sorted(values)


def source_input_paths(args: argparse.Namespace, label: str) -> Tuple[str, Path, Optional[Path], Optional[Path]]:
    if args.source_input_mode in ("auto", "stageb"):
        stageb_root = args.inputs_dir / args.stageb_dirname / args.ratio
        zns_path = stageb_root / f"{label}_zns_track_steps.csv"
        anchor_path = stageb_root / f"{label}_capture_anchors.csv"
        unexpected_path = stageb_root / f"{label}_unexpected_boundary_exits.csv"
        if zns_path.exists():
            if not anchor_path.exists():
                raise FileNotFoundError(f"Missing StageB capture anchors file: {anchor_path}")
            if not unexpected_path.exists():
                raise FileNotFoundError(
                    f"Missing StageB unexpected boundary exits file: {unexpected_path}"
                )
            return "stageb", zns_path, anchor_path, unexpected_path
        if args.source_input_mode == "stageb":
            raise FileNotFoundError(f"Missing StageB ZnS track file: {zns_path}")

    alpha_path = args.inputs_dir / args.alpha_li_dirname / args.ratio / f"{label}_alpha_li_steps.csv"
    if alpha_path.exists():
        return "alpha-li", alpha_path, None, None
    raise FileNotFoundError(f"Missing source input for thickness {label}: {alpha_path}")


def command_arg(value: object, *, executable: bool = False) -> str:
    if isinstance(value, Path):
        if executable and not value.is_absolute() and value.parent == Path("."):
            return str(value.resolve())
        return str(value)
    return str(value)


def run_command(cmd: Sequence[object], dry_run: bool) -> None:
    argv = [command_arg(c, executable=(i == 0)) for i, c in enumerate(cmd)]
    print(" ".join(argv), flush=True)
    if not dry_run:
        subprocess.run(argv, check=True)


def write_thickness_light_summary(run_root: Path, thicknesses: Sequence[float]) -> Optional[Path]:
    rows: List[Dict[str, str]] = []
    fieldnames: List[str] = []
    preferred = [
        "thickness_um",
        "ratio_tag",
        "wavelength_nm",
        "mu_a_per_um",
        "mu_s_per_um",
        "g",
        "n_eff",
        "transport_scattering_mode",
        "mu_s_input_per_um",
        "g_input",
        "mu_s_prime_input_per_um",
        "phase_function_mode",
        "phase_function_csv",
        "front_reflection_model",
        "front_reflectance",
        "front_reflection_mode",
        "front_aluminum_n",
        "front_aluminum_k",
        "back_reflection_model",
        "back_air_n",
        "mu_s_prime_per_um",
        "mu_s_prime_from_g_per_um",
        "mu_tr_per_um",
        "absorption_length_um",
        "scattering_length_um",
        "transport_mfp_um",
        "diffusion_length_um",
        "n_events",
        "incident_event_count",
        "capture_fraction",
        "n_source_steps",
        "samples_per_step",
        "total_source_weight",
        "total_detected_weight",
        "front_escape_weight",
        "back_escape_weight",
        "absorbed_weight",
        "lost_weight",
        "mean_light_per_capture",
        "mean_detected_light_per_capture",
        "mean_light_per_incident",
        "mean_detected_light_per_incident",
        "detection_efficiency",
        "spot_rms_x",
        "spot_rms_y",
        "spot_rms_r",
        "fwhm_x",
        "fwhm_y",
    ]

    for thickness in sorted(thicknesses):
        summary_path = run_root / thickness_label(thickness) / "optical_mc_summary.csv"
        if not summary_path.exists():
            continue
        with summary_path.open(newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            try:
                row = next(reader)
            except StopIteration:
                continue
        row["thickness_um"] = row.get("thickness_um") or thickness_label(thickness)
        rows.append(row)
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    if not rows:
        return None

    ordered = [k for k in preferred if k in fieldnames] + [
        k for k in fieldnames if k not in preferred
    ]
    out_path = run_root / "thickness_light_summary.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=ordered, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in ordered})
    return out_path


def main() -> int:
    args = parse_args()
    if args.transport_scattering_mode == "reduced-isotropic" and (
        args.phase_function_csv is not None or args.scattering_model == "tabulated"
    ):
        print(
            "warning: reduced-isotropic transport ignores --phase-function-csv and "
            "--scattering-model tabulated; using g=0 and no phase function.",
            file=sys.stderr,
        )
    run_root = run_root_for_args(args)
    source_root = args.inputs_dir / args.generated_dirname / args.ratio
    thicknesses = discover_thicknesses(args)

    try:
        optical_properties = build_optical_properties(args, run_root)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"ratio,{args.ratio}", flush=True)
    print(f"thicknesses,{','.join(thickness_label(t) for t in thicknesses)}", flush=True)
    print(f"n_eff,{effective_index(args):.12g}", flush=True)
    print(f"transport_scattering_mode,{args.transport_scattering_mode}", flush=True)
    print(f"run_root,{run_root}", flush=True)
    print(f"optical_properties,{optical_properties}", flush=True)

    for thickness in thicknesses:
        label = thickness_label(thickness)
        step_sources = source_root / f"{label}_macro_zns_step_sources.csv"
        event_sources = source_root / f"{label}_event_light_sources.csv"
        out_dir = run_root / label
        config_path = out_dir / "run_config.generated.json"

        try:
            input_mode, source_input, capture_anchors, unexpected_exits = source_input_paths(
                args, label
            )
        except Exception as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

        generated_sources_this_run = False
        if not args.mc_only and (args.overwrite_sources or not step_sources.exists() or not event_sources.exists()):
            cmd = [
                sys.executable,
                args.preprocessor,
                source_input,
                "--output-dir",
                source_root,
                "--ratio-tag",
                args.ratio,
                "--thickness-um",
                label,
                "--xy-anchor-mode",
                args.xy_anchor_mode,
                "--yield-zns-per-MeV",
                str(args.yield_zns_per_MeV),
                "--wavelength-nm",
                str(args.wavelength_nm),
                "--quench-model",
                args.quench_model,
                "--birks-kb-um-per-keV",
                str(args.birks_kb_um_per_keV),
            ]
            if input_mode == "stageb":
                cmd.extend(
                    [
                        "--capture-anchors",
                        capture_anchors,
                        "--unexpected-boundary-exits",
                        unexpected_exits,
                    ]
                )
            run_command(cmd, args.dry_run)
            generated_sources_this_run = input_mode == "stageb" and not args.dry_run

        if args.preprocess_only:
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        config = {
            "ratio_tag": args.ratio,
            "thickness_um": thickness,
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
            "output_detected_photons": args.output_detected_photons,
            "psf_bin_size_um": args.psf_bin_size_um,
            "psf_range_um": args.psf_range_um,
            "lsf_range_um": args.lsf_range_um,
            "optical_properties_csv": str(optical_properties),
            "source_steps_csv": str(step_sources),
            "event_sources_csv": str(event_sources),
            "output_dir": str(out_dir),
            "wavelength_nm": args.wavelength_nm,
        }
        if not args.dry_run:
            with config_path.open("w", encoding="utf-8") as f:
                json.dump(config, f, indent=2)
                f.write("\n")
        run_command([args.opticalmc, config_path], args.dry_run)
        if (
            generated_sources_this_run
            and not args.keep_generated_sources
            and not args.dry_run
        ):
            for generated_path in (step_sources, event_sources):
                try:
                    generated_path.unlink()
                    print(f"removed_generated_source,{generated_path}", flush=True)
                except FileNotFoundError:
                    pass

    if not args.preprocess_only and not args.dry_run:
        summary_path = write_thickness_light_summary(run_root, thicknesses)
        if summary_path is not None:
            print(f"thickness_light_summary,{summary_path}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
