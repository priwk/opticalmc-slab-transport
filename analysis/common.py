from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


READOUT_LIGHT_SCALE_DEFAULT = 1_000_000.0


def thickness_from_dir(path: Path) -> Optional[float]:
    try:
        return float(path.name)
    except ValueError:
        return None


def thickness_label(value: float) -> str:
    return f"{value:.12g}"


def configure_plot_style(font_dir: Optional[Path] = None) -> None:
    """Configure Matplotlib for boxed, inward-tick, local-font figures."""
    import matplotlib as mpl
    from matplotlib import font_manager

    font_dir = font_dir or Path(__file__).resolve().parent / "fonts"
    local_fonts = []
    for pattern in ("*.ttf", "*.TTF", "*.ttc", "*.TTC", "*.otf", "*.OTF"):
        local_fonts.extend(font_dir.glob(pattern))
    local_font_paths = {font_path.resolve() for font_path in local_fonts if font_path.exists()}
    local_font_names = set()
    for font_path in local_fonts:
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
            local_font_names.add(font_manager.FontProperties(fname=str(font_path)).get_name())

    if local_font_names:
        font_manager.fontManager.ttflist = [
            entry
            for entry in font_manager.fontManager.ttflist
            if entry.name not in local_font_names or Path(entry.fname).resolve() in local_font_paths
        ]
        if hasattr(font_manager, "_findfont_cached"):
            font_manager._findfont_cached.cache_clear()

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DengXian"],
            "axes.unicode_minus": False,
            "axes.linewidth": 1.0,
            "axes.spines.top": True,
            "axes.spines.right": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def style_axes(ax, square: bool = False, grid: bool = False) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)
    ax.tick_params(direction="in", top=True, right=True, which="both")
    if grid:
        ax.grid(True, color="#d9d9d9", linewidth=0.6, alpha=0.7)
    else:
        ax.grid(False)
    if square and hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect(1)


def style_figure_axes(fig, square: bool = False, grid: bool = False) -> None:
    for ax in fig.axes:
        style_axes(ax, square=square, grid=grid)


def add_plot_metrics(df: pd.DataFrame, light_scale: float = READOUT_LIGHT_SCALE_DEFAULT) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    numeric_columns = [
        "thickness_um",
        "mean_detected_light_per_incident",
        "mean_detected_light_per_capture",
        "mean_light_per_incident",
        "mean_light_per_capture",
        "incident_event_count",
        "detection_efficiency",
        "capture_fraction",
        "spot_rms_x",
        "spot_rms_y",
        "spot_rms_r",
        "fwhm_x",
        "fwhm_y",
    ]
    for column in numeric_columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")

    if "fwhm_mean" not in out.columns and {"fwhm_x", "fwhm_y"}.issubset(out.columns):
        out["fwhm_mean"] = out[["fwhm_x", "fwhm_y"]].mean(axis=1)
    if "detection_efficiency_percent" not in out.columns and "detection_efficiency" in out.columns:
        out["detection_efficiency_percent"] = 100.0 * out["detection_efficiency"]
    if "capture_fraction_percent" not in out.columns and "capture_fraction" in out.columns:
        out["capture_fraction_percent"] = 100.0 * out["capture_fraction"]

    light_col = detected_light_column(out)
    if light_col is not None:
        out["readout_light_per_incident"] = out[light_col]
        if light_col == "mean_detected_light_per_incident":
            out["readout_light_per_million_incident"] = out[light_col] * light_scale
        elif "incident_event_count" in out.columns:
            out["readout_light_per_million_incident"] = (
                out[light_col] * out.get("capture_fraction", 1.0) * light_scale
            )
        else:
            out["readout_light_per_million_incident"] = out[light_col] * light_scale
        max_value = out["readout_light_per_million_incident"].max(skipna=True)
        out["readout_light_relative"] = (
            out["readout_light_per_million_incident"] / max_value if max_value and max_value > 0 else np.nan
        )
    return out


def detected_light_column(df: pd.DataFrame) -> Optional[str]:
    if "mean_detected_light_per_incident" in df.columns:
        return "mean_detected_light_per_incident"
    if "mean_detected_light_per_capture" in df.columns:
        return "mean_detected_light_per_capture"
    return None


def select_panel_thicknesses(values: Iterable[float], max_count: int = 6) -> List[float]:
    unique = sorted({float(v) for v in values if pd.notna(v)})
    if len(unique) <= max_count:
        return unique
    indices = np.linspace(0, len(unique) - 1, max_count).round().astype(int)
    return [unique[int(i)] for i in sorted(set(indices))]


def resolve_existing_path(path_value: object, base_dir: Path) -> Optional[Path]:
    if path_value is None or (isinstance(path_value, float) and math.isnan(path_value)):
        return None
    text = str(path_value).strip()
    if not text:
        return None
    path = Path(text)
    candidates = [path]
    if not path.is_absolute():
        candidates.append(base_dir / path)
        candidates.append(Path.cwd() / path)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def load_phase_function(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if {"cos_theta_min", "cos_theta_max"}.issubset(df.columns):
        cos_theta_min = pd.to_numeric(df["cos_theta_min"], errors="coerce")
        cos_theta_max = pd.to_numeric(df["cos_theta_max"], errors="coerce")
        cos_theta = 0.5 * (cos_theta_min + cos_theta_max)
    elif "cos_theta" in df.columns:
        cos_theta_min = None
        cos_theta_max = None
        cos_theta = pd.to_numeric(df["cos_theta"], errors="coerce")
    elif "mu" in df.columns:
        cos_theta_min = None
        cos_theta_max = None
        cos_theta = pd.to_numeric(df["mu"], errors="coerce")
    else:
        raise ValueError(f"{path} does not contain cos(theta) columns")

    for column in ("probability_density_mean", "probability_density", "probability_mean", "probability", "weight"):
        if column in df.columns:
            values = pd.to_numeric(df[column], errors="coerce")
            break
    else:
        raise ValueError(f"{path} does not contain phase-function probability columns")

    out = pd.DataFrame({"cos_theta": cos_theta, "phase_function": values})
    if cos_theta_min is not None and cos_theta_max is not None:
        out["cos_theta_min"] = cos_theta_min
        out["cos_theta_max"] = cos_theta_max
    return out.dropna().sort_values("cos_theta").reset_index(drop=True)


def phase_function_with_theta(phase: pd.DataFrame) -> pd.DataFrame:
    out = phase.copy()
    mu = np.clip(out["cos_theta"].to_numpy(dtype=float), -1.0, 1.0)
    theta_rad = np.arccos(mu)
    out["theta_deg"] = np.degrees(theta_rad)
    out["phase_function_theta"] = out["phase_function"].to_numpy(dtype=float) * np.sin(theta_rad)
    return out.sort_values("theta_deg").reset_index(drop=True)


def phase_function_polar_curve(phase: pd.DataFrame, samples_per_bin: int = 12) -> tuple[np.ndarray, np.ndarray]:
    if {"cos_theta_min", "cos_theta_max"}.issubset(phase.columns):
        rows = []
        for row in phase.itertuples(index=False):
            mu_min = float(getattr(row, "cos_theta_min"))
            mu_max = float(getattr(row, "cos_theta_max"))
            radius = float(getattr(row, "phase_function"))
            theta_low = math.acos(float(np.clip(max(mu_min, mu_max), -1.0, 1.0)))
            theta_high = math.acos(float(np.clip(min(mu_min, mu_max), -1.0, 1.0)))
            rows.append((theta_low, theta_high, radius))
        rows.sort(key=lambda item: item[0])

        theta_pos: List[float] = []
        radius_pos: List[float] = []
        for theta_low, theta_high, radius in rows:
            segment = np.linspace(theta_low, theta_high, max(2, samples_per_bin))
            theta_pos.extend(segment.tolist())
            radius_pos.extend([radius] * len(segment))
            theta_pos.append(theta_high)
            radius_pos.append(radius)
        theta_arr = np.asarray(theta_pos, dtype=float)
        radius_arr = np.asarray(radius_pos, dtype=float)
    else:
        theta_arr = np.arccos(np.clip(phase["cos_theta"].to_numpy(dtype=float), -1.0, 1.0))
        radius_arr = phase["phase_function"].to_numpy(dtype=float)
        order = np.argsort(theta_arr)
        theta_arr = theta_arr[order]
        radius_arr = radius_arr[order]
        if len(theta_arr):
            theta_arr = np.r_[0.0, theta_arr, math.pi]
            radius_arr = np.r_[radius_arr[0], radius_arr, radius_arr[-1]]

    angles = np.concatenate((-theta_arr[:0:-1], theta_arr))
    radii = np.concatenate((radius_arr[:0:-1], radius_arr))
    return angles, radii


RUN_TAG_SEPARATOR = "__"


def run_tag(ratio: str, run_label: Optional[str] = None) -> str:
    if run_label:
        return f"{ratio}{RUN_TAG_SEPARATOR}{run_label}"
    return ratio


def display_run_tag(tag: str, hidden_run_labels: Optional[Iterable[str]] = None) -> str:
    hidden = set(hidden_run_labels or ["modeB_anisotropic"])
    for separator in (RUN_TAG_SEPARATOR, "/"):
        if separator not in tag:
            continue
        ratio, run_label = tag.split(separator, 1)
        if run_label in hidden:
            return ratio
    return tag


def run_filter_matches(ratio_filter: set[str], ratio: str, run_label: Optional[str]) -> bool:
    if not ratio_filter:
        return True
    candidates = {ratio, run_tag(ratio, run_label)}
    if run_label:
        candidates.add(run_label)
        candidates.add(f"{ratio}/{run_label}")
    return not candidates.isdisjoint(ratio_filter)


def add_runs_from_root(rows: List[Dict[str, object]], ratio: str, run_label: Optional[str], run_root: Path) -> None:
    tag = run_tag(ratio, run_label)
    for thickness_dir in sorted(p for p in run_root.iterdir() if p.is_dir()):
        thickness = thickness_from_dir(thickness_dir)
        if thickness is None:
            continue
        summary = thickness_dir / "optical_mc_summary.csv"
        if summary.exists():
            rows.append(
                {
                    "ratio_tag": tag,
                    "base_ratio_tag": ratio,
                    "run_label": run_label or "",
                    "thickness_um": thickness,
                    "run_dir": str(thickness_dir),
                    "summary_csv": str(summary),
                    "event_summary_csv": str(thickness_dir / "optical_mc_event_summary.csv"),
                    "step_summary_csv": str(thickness_dir / "optical_mc_source_step_summary.csv"),
                    "psf_2d_csv": str(thickness_dir / "psf_2d.csv"),
                    "lsf_x_csv": str(thickness_dir / "lsf_x.csv"),
                    "lsf_y_csv": str(thickness_dir / "lsf_y.csv"),
                }
            )


def discover_runs(outputs_dir: Path, ratios: Optional[Iterable[str]] = None) -> pd.DataFrame:
    ratio_filter = set(ratios or [])
    rows: List[Dict[str, object]] = []
    if not outputs_dir.exists():
        return pd.DataFrame()

    for ratio_dir in sorted(p for p in outputs_dir.iterdir() if p.is_dir()):
        ratio = ratio_dir.name
        if ratio.startswith("_") and not ratio_filter:
            continue
        if not run_filter_matches(ratio_filter, ratio, None):
            nested_requested = any(
                item.startswith(f"{ratio}{RUN_TAG_SEPARATOR}") or item.startswith(f"{ratio}/")
                for item in ratio_filter
            )
            if not nested_requested:
                continue

        if run_filter_matches(ratio_filter, ratio, None):
            add_runs_from_root(rows, ratio, None, ratio_dir)

        for run_dir in sorted(p for p in ratio_dir.iterdir() if p.is_dir()):
            if thickness_from_dir(run_dir) is not None:
                continue
            run_label = run_dir.name
            if run_filter_matches(ratio_filter, ratio, run_label):
                add_runs_from_root(rows, ratio, run_label, run_dir)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["ratio_tag", "thickness_um"]).reset_index(drop=True)


def read_summary(path: Path) -> Dict[str, object]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        row = next(csv.DictReader(f))
    out: Dict[str, object] = {}
    for key, value in row.items():
        if value is None or value == "":
            out[key] = value
            continue
        try:
            out[key] = float(value)
        except ValueError:
            out[key] = value
    return out


def load_lsf(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    center_col = df.columns[0]
    return df.rename(columns={center_col: "position_um", "weight": "weight"})


def normalized_mtf_from_lsf(lsf: pd.DataFrame) -> pd.DataFrame:
    if lsf.empty:
        return pd.DataFrame(columns=["frequency_lp_per_mm", "mtf"])
    pos = lsf["position_um"].to_numpy(dtype=float)
    weight = lsf["weight"].to_numpy(dtype=float)
    if len(pos) < 2 or np.all(weight == 0):
        return pd.DataFrame(columns=["frequency_lp_per_mm", "mtf"])

    spacing_um = float(np.median(np.diff(pos)))
    signal = np.clip(weight, 0.0, None)
    total = np.sum(signal)
    if total <= 0:
        return pd.DataFrame(columns=["frequency_lp_per_mm", "mtf"])

    mtf = np.abs(np.fft.rfft(signal))
    if mtf[0] > 0:
        mtf = mtf / mtf[0]
    freq_cycles_per_um = np.fft.rfftfreq(len(signal), d=spacing_um)
    return pd.DataFrame(
        {
            "frequency_lp_per_mm": freq_cycles_per_um * 1000.0,
            "mtf": mtf,
        }
    )


def crossing_frequency(mtf: pd.DataFrame, level: float) -> float:
    if mtf.empty:
        return math.nan
    freq = mtf["frequency_lp_per_mm"].to_numpy(dtype=float)
    values = mtf["mtf"].to_numpy(dtype=float)
    for i in range(1, len(values)):
        if values[i] <= level <= values[i - 1]:
            if abs(values[i] - values[i - 1]) < 1.0e-15:
                return float(freq[i])
            t = (level - values[i - 1]) / (values[i] - values[i - 1])
            return float(freq[i - 1] + t * (freq[i] - freq[i - 1]))
    return math.nan


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path
