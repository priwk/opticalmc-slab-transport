from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


def thickness_from_dir(path: Path) -> Optional[float]:
    try:
        return float(path.name)
    except ValueError:
        return None


def thickness_label(value: float) -> str:
    return f"{value:.12g}"


def discover_runs(outputs_dir: Path, ratios: Optional[Iterable[str]] = None) -> pd.DataFrame:
    ratio_filter = set(ratios or [])
    rows: List[Dict[str, object]] = []
    if not outputs_dir.exists():
        return pd.DataFrame()

    for ratio_dir in sorted(p for p in outputs_dir.iterdir() if p.is_dir()):
        ratio = ratio_dir.name
        if ratio_filter and ratio not in ratio_filter:
            continue
        for thickness_dir in sorted(p for p in ratio_dir.iterdir() if p.is_dir()):
            thickness = thickness_from_dir(thickness_dir)
            if thickness is None:
                continue
            summary = thickness_dir / "optical_mc_summary.csv"
            if summary.exists():
                rows.append(
                    {
                        "ratio_tag": ratio,
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
    signal = weight - np.min(weight)
    total = np.sum(signal)
    if total <= 0:
        signal = weight
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
