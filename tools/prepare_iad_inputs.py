#!/usr/bin/env python3
"""Prepare glass-corrected IAD input CSV files from raw spectra."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "从原始 R-D/T-D 光谱整理 IAD 输入。默认按同波长扣除 Glass-R-D/Glass-T-D，"
            "再把百分数转换为 0-1 分数。"
        )
    )
    parser.add_argument("--raw-dir", type=Path, default=Path("inputs/511"), help="原始光谱目录。")
    parser.add_argument("--output-dir", type=Path, default=Path("inputs/iad_inputs/511"), help="整理后的 IAD 输入目录。")
    return parser.parse_args()


def parse_measurement(path: Path) -> List[Tuple[float, float]]:
    rows: List[Tuple[float, float]] = []
    text = path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        text_line = line.strip().replace("\ufeff", "")
        match = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*[,\t ]\s*([-+0-9.]+)", text_line)
        if match:
            rows.append((float(match.group(1)), float(match.group(2))))
    if not rows:
        raise ValueError(f"没有从文件中解析到数值数据: {path}")
    return rows


def by_wavelength(rows: List[Tuple[float, float]]) -> Dict[float, float]:
    return {wavelength: value for wavelength, value in rows}


def value_at(rows: List[Tuple[float, float]], wavelength: float) -> float:
    values = by_wavelength(rows)
    return values[float(wavelength)]


def fmt(value: float) -> str:
    return f"{value:.10g}"


def fmt_wl(value: float) -> str:
    return f"{value:.1f}"


def resolve_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def write_readme(out_dir: Path) -> None:
    readme = """# IAD 输入数据：511 日期批次

本目录保存从原始 `inputs/511` 文本文件整理出的漫反射/漫透射光谱数据，
供后续 IAD 反演使用。这里的 `511` 是数据集/日期批次标签，不是默认反演波长。

样品数据已经做玻璃参考扣除：同一波长下用样品的漫反射/漫透射百分数分别减去
`Glass-R-D.txt` 和 `Glass-T-D.txt` 的百分数，然后再除以 100 得到 IAD 输入用的
`MR` 和 `MT`。

## 重新生成

如果本目录被清空，先确认原始文件还在 `inputs/511`，然后运行：

```powershell
python .\\tools\\prepare_iad_inputs.py
```

## 文件说明

- `samples.csv`：每个样品/厚度一行，包含扣玻璃后的 450 nm 和 511 nm 快速检查值。
- `measurements.csv`：所有样品光谱合并后的长表。
- `per_sample/*_iad_input.csv`：每个样品一个清洗后的光谱文件，包含
  `wavelength_nm`、`MR`、`MT` 等列，供 IAD 包装脚本读取。
- `references/glass_diffuse_reference.csv`：玻璃参考片的漫反射/漫透射光谱。
- `manifest.json`：机器可读的数据集摘要。
"""
    (out_dir / "README.md").write_text(readme, encoding="utf-8")


def main() -> int:
    args = parse_args()
    raw_dir = resolve_path(args.raw_dir)
    out_dir = resolve_path(args.output_dir)
    per_sample_dir = out_dir / "per_sample"
    ref_dir = out_dir / "references"
    per_sample_dir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)

    glass_r_path = raw_dir / "Glass-R-D.txt"
    glass_t_path = raw_dir / "Glass-T-D.txt"
    if not glass_r_path.exists() or not glass_t_path.exists():
        raise FileNotFoundError(
            f"缺少玻璃参考文件。需要同时存在: {glass_r_path} 和 {glass_t_path}"
        )

    glass_r_rows = parse_measurement(glass_r_path)
    glass_t_rows = parse_measurement(glass_t_path)
    if [wl for wl, _ in glass_r_rows] != [wl for wl, _ in glass_t_rows]:
        raise ValueError("Glass-R-D 和 Glass-T-D 的波长网格不一致。")

    glass_r = by_wavelength(glass_r_rows)
    glass_t = by_wavelength(glass_t_rows)
    with (ref_dir / "glass_diffuse_reference.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "wavelength_nm",
                "glass_diffuse_reflectance",
                "glass_diffuse_transmittance",
                "glass_reflectance_percent",
                "glass_transmittance_percent",
                "source_R_file",
                "source_T_file",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for wavelength, r_percent in glass_r_rows:
            t_percent = glass_t[wavelength]
            writer.writerow(
                {
                    "wavelength_nm": fmt_wl(wavelength),
                    "glass_diffuse_reflectance": fmt(r_percent / 100),
                    "glass_diffuse_transmittance": fmt(t_percent / 100),
                    "glass_reflectance_percent": fmt(r_percent),
                    "glass_transmittance_percent": fmt(t_percent),
                    "source_R_file": str(glass_r_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                    "source_T_file": str(glass_t_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                }
            )

    sample_pattern = re.compile(r"^(?P<ratio>.+)-(?P<thick>[0-9]+(?:\.[0-9]+)?)-(?P<kind>[RT])-D\.txt$", re.I)
    samples_by_id: Dict[str, Dict[str, object]] = {}
    for path in raw_dir.glob("*.txt"):
        if path.name.startswith("Glass"):
            continue
        match = sample_pattern.match(path.name)
        if not match:
            continue
        sample_id = f"{match.group('ratio')}-{match.group('thick')}"
        samples_by_id.setdefault(
            sample_id,
            {"ratio_tag": match.group("ratio"), "thickness_um": float(match.group("thick"))},
        )[match.group("kind").upper()] = path

    if not samples_by_id:
        raise FileNotFoundError(f"没有在 {raw_dir} 找到形如 X-X-100-R-D.txt/T-D.txt 的样品文件。")

    samples: List[Dict[str, object]] = []
    all_rows: List[Dict[str, str]] = []
    invalid_total = 0
    for sample_id, meta in sorted(samples_by_id.items(), key=lambda item: (float(item[1]["thickness_um"]), item[0])):
        if "R" not in meta or "T" not in meta:
            raise ValueError(f"样品缺少 R/T 成对文件: {sample_id}")

        r_path = meta["R"]
        t_path = meta["T"]
        assert isinstance(r_path, Path)
        assert isinstance(t_path, Path)
        r_rows = parse_measurement(r_path)
        t_rows = parse_measurement(t_path)
        wavelengths = [wavelength for wavelength, _ in r_rows]
        if wavelengths != [wavelength for wavelength, _ in t_rows]:
            raise ValueError(f"R/T 波长网格不一致: {sample_id}")
        if any(wavelength not in glass_r or wavelength not in glass_t for wavelength in wavelengths):
            raise ValueError(f"玻璃参考缺少样品波长点: {sample_id}")

        per_sample_path = per_sample_dir / f"{sample_id}_iad_input.csv"
        corrected_r_rows: List[Tuple[float, float]] = []
        corrected_t_rows: List[Tuple[float, float]] = []
        invalid_count = 0
        with per_sample_path.open("w", newline="", encoding="utf-8") as handle:
            fieldnames = [
                "wavelength_nm",
                "MR",
                "MT",
                "MR_corrected_percent",
                "MT_corrected_percent",
                "MR_raw_percent",
                "MT_raw_percent",
                "glass_R_percent",
                "glass_T_percent",
                "valid_iad_input",
                "ratio_tag",
                "thickness_um",
                "sample_id",
                "correction_method",
                "source_R_file",
                "source_T_file",
                "glass_R_file",
                "glass_T_file",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            for (wavelength, r_raw_percent), (_, t_raw_percent) in zip(r_rows, t_rows):
                r_corrected_percent = r_raw_percent - glass_r[wavelength]
                t_corrected_percent = t_raw_percent - glass_t[wavelength]
                is_valid = (
                    0.0 <= r_corrected_percent <= 100.0
                    and 0.0 <= t_corrected_percent <= 100.0
                    and r_corrected_percent + t_corrected_percent <= 100.0
                )
                if not is_valid:
                    invalid_count += 1
                corrected_r_rows.append((wavelength, r_corrected_percent))
                corrected_t_rows.append((wavelength, t_corrected_percent))
                row = {
                    "wavelength_nm": fmt_wl(wavelength),
                    "MR": fmt(r_corrected_percent / 100),
                    "MT": fmt(t_corrected_percent / 100),
                    "MR_corrected_percent": fmt(r_corrected_percent),
                    "MT_corrected_percent": fmt(t_corrected_percent),
                    "MR_raw_percent": fmt(r_raw_percent),
                    "MT_raw_percent": fmt(t_raw_percent),
                    "glass_R_percent": fmt(glass_r[wavelength]),
                    "glass_T_percent": fmt(glass_t[wavelength]),
                    "valid_iad_input": "1" if is_valid else "0",
                    "ratio_tag": str(meta["ratio_tag"]),
                    "thickness_um": f"{float(meta['thickness_um']):.12g}",
                    "sample_id": sample_id,
                    "correction_method": "sample_minus_glass_same_wavelength",
                    "source_R_file": str(r_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                    "source_T_file": str(t_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                    "glass_R_file": str(glass_r_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                    "glass_T_file": str(glass_t_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                }
                writer.writerow(row)
                all_rows.append(row)

        invalid_total += invalid_count
        samples.append(
            {
                "sample_id": sample_id,
                "ratio_tag": str(meta["ratio_tag"]),
                "thickness_um": f"{float(meta['thickness_um']):.12g}",
                "wavelength_min_nm": fmt_wl(wavelengths[0]),
                "wavelength_max_nm": fmt_wl(wavelengths[-1]),
                "wavelength_step_nm": fmt_wl(wavelengths[1] - wavelengths[0]) if len(wavelengths) > 1 else "",
                "n_wavelengths": len(wavelengths),
                "invalid_iad_rows": invalid_count,
                "MR_at_450_nm": fmt(value_at(corrected_r_rows, 450.0) / 100),
                "MT_at_450_nm": fmt(value_at(corrected_t_rows, 450.0) / 100),
                "MR_at_511_nm": fmt(value_at(corrected_r_rows, 511.0) / 100),
                "MT_at_511_nm": fmt(value_at(corrected_t_rows, 511.0) / 100),
                "correction_method": "sample_minus_glass_same_wavelength",
                "iad_input_csv": str(per_sample_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "source_R_file": str(r_path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "source_T_file": str(t_path.relative_to(REPO_ROOT)).replace("\\", "/"),
            }
        )

    with (out_dir / "samples.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "sample_id",
            "ratio_tag",
            "thickness_um",
            "wavelength_min_nm",
            "wavelength_max_nm",
            "wavelength_step_nm",
            "n_wavelengths",
            "invalid_iad_rows",
            "MR_at_450_nm",
            "MT_at_450_nm",
            "MR_at_511_nm",
            "MT_at_511_nm",
            "correction_method",
            "iad_input_csv",
            "source_R_file",
            "source_T_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(samples)

    with (out_dir / "measurements.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "sample_id",
            "ratio_tag",
            "thickness_um",
            "wavelength_nm",
            "MR",
            "MT",
            "MR_corrected_percent",
            "MT_corrected_percent",
            "MR_raw_percent",
            "MT_raw_percent",
            "glass_R_percent",
            "glass_T_percent",
            "valid_iad_input",
            "correction_method",
            "source_R_file",
            "source_T_file",
            "glass_R_file",
            "glass_T_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in all_rows:
            writer.writerow({key: row[key] for key in fieldnames})

    manifest = {
        "dataset_id": out_dir.name,
        "description": "Glass-corrected diffuse reflectance/transmittance spectra prepared as IAD input data.",
        "raw_dir": str(raw_dir.relative_to(REPO_ROOT)).replace("\\", "/"),
        "output_dir": str(out_dir.relative_to(REPO_ROOT)).replace("\\", "/"),
        "sample_count": len(samples),
        "measurement_rows": len(all_rows),
        "reference_rows": len(glass_r_rows),
        "invalid_iad_rows_total": invalid_total,
        "correction_method": "MR_percent = sample_R_percent - Glass_R_percent; MT_percent = sample_T_percent - Glass_T_percent; MR/MT are corrected_percent / 100.",
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    write_readme(out_dir)

    print(f"output_dir,{out_dir}")
    print(f"samples,{len(samples)}")
    print(f"measurement_rows,{len(all_rows)}")
    print(f"invalid_iad_rows,{invalid_total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
