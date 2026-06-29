#!/usr/bin/env python3
"""Simple front door for common OpticalMC batch runs."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List


REPO_ROOT = Path(__file__).resolve().parent

PRESETS = {
    "quick": {
        "samples_per_step": 4,
        "psf_bin_size_um": 2.0,
        "psf_range_um": 300.0,
        "lsf_range_um": 2000.0,
        "description": "快速检查流程和参数，统计噪声较大。",
    },
    "normal": {
        "samples_per_step": 32,
        "psf_bin_size_um": 2.0,
        "psf_range_um": 500.0,
        "lsf_range_um": 5000.0,
        "description": "日常批量结果，速度和统计质量较均衡。",
    },
    "fine": {
        "samples_per_step": 128,
        "psf_bin_size_um": 2.0,
        "psf_range_um": 500.0,
        "lsf_range_um": 5000.0,
        "description": "更高统计质量，适合最终图表，耗时更长。",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "傻瓜式批量入口。例：python .\\run.py 1-2 30 40 50；"
            "不写厚度则运行该配比的全部厚度。"
        )
    )
    parser.add_argument("ratio", help="BN:ZnS 质量配比目录名，例如 2-1、1-1、1-2.5、1-3。")
    parser.add_argument(
        "thickness",
        nargs="*",
        help="可选厚度列表，单位 um。例如 30 40 50；不填则跑该配比全部厚度。",
    )
    parser.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default="normal",
        help="运行预设：quick 快速检查，normal 日常批量，fine 高统计。",
    )
    parser.add_argument(
        "--coords",
        choices=("global",),
        default="global",
        help="坐标模式。仅支持 global，使用 capture_x/y 保留屏幕全局坐标。",
    )
    parser.add_argument(
        "--readout",
        choices=("back", "front", "both"),
        default="back",
        help="读出表面。",
    )
    parser.add_argument(
        "--front-reflection-model",
        choices=("effective", "aluminum_fresnel"),
        default="aluminum_fresnel",
        help="前表面反射概率模型：effective 使用常数反射率；aluminum_fresnel 用 n_eff 到铝 n+ik 的角度相关 Fresnel 反射率。",
    )
    parser.add_argument(
        "--front-reflectance",
        type=float,
        default=0.0,
        help="前表面铝板有效反射率；用于 --front-reflection-model effective。",
    )
    parser.add_argument(
        "--front-reflection-mode",
        choices=("none", "specular", "diffuse"),
        default="specular",
        help="前表面边界：none 为旧的逃逸边界；specular/diffuse 用于模拟铝板反射。",
    )
    parser.add_argument("--front-aluminum-n", type=float, default=0.65, help="aluminum_fresnel 使用的铝复折射率实部。")
    parser.add_argument("--front-aluminum-k", type=float, default=5.3, help="aluminum_fresnel 使用的铝消光系数。")
    parser.add_argument(
        "--back-reflection-model",
        choices=("none", "air_fresnel"),
        default="air_fresnel",
        help="后表面出屏边界：air_fresnel 用 n_eff 到空气的 Fresnel 反射率；none 为旧的直接出界。",
    )
    parser.add_argument("--back-air-n", type=float, default=1.000293, help="后表面外侧空气折射率。")
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 1),
        help="线程数，默认使用 CPU 核心数减一。",
    )
    parser.add_argument(
        "--yield-zns-per-MeV",
        type=float,
        default=60000.0,
        help="ZnS 光产额，单位 photons/MeV。",
    )
    parser.add_argument(
        "--incident-events",
        type=float,
        default=None,
        help="入射 neutron history 数，用于 per-incident 出光量归一。默认 1000000。",
    )
    parser.add_argument(
        "--quench-kb",
        type=float,
        default=None,
        help="启用 Birks 淬灭并设置 kB，单位 um/keV。不填则不淬灭。",
    )
    parser.add_argument(
        "--mu-a-scale",
        type=float,
        default=1.0,
        help="光学吸收系数缩放，用于敏感性检查。",
    )
    parser.add_argument(
        "--optical-component",
        choices=("bulk", "total", "boundary"),
        default="bulk",
        help="选择 merged optical params 中的 mu_s_prime/g 来源。默认 bulk 使用体散射主列；total/boundary 用于对照。",
    )
    parser.add_argument(
        "--mu-s-scale",
        type=float,
        default=1.0,
        help="光学散射强度缩放；各向异性模式下作用于 mu_s_prime，并同步派生 mu_s。",
    )
    parser.add_argument(
        "--phase-function-csv",
        type=Path,
        default=None,
        help="可选 phase_function.csv；需同时指定 --scattering-model tabulated 才会使用。",
    )
    parser.add_argument(
        "--scattering-model",
        choices=("auto", "tabulated", "hg"),
        default="hg",
        help="散射模型：默认 hg 使用 HG(g)；tabulated 显式使用表格相函数；auto 等同 HG(g)。",
    )
    parser.add_argument(
        "--transport-scattering-mode",
        choices=("reduced-isotropic", "anisotropic"),
        default="anisotropic",
        help="输运散射模式：默认 anisotropic，使用 StageD 的 mu_s_prime 和 g 计算 mu_s。",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=100000,
        help="每个光子最大相互作用步数。默认 100000，避免强散射参数下过早计入 lost_weight。",
    )
    parser.add_argument(
        "--run-label",
        default=None,
        help="可选输出子目录名，写到 outputs/<ratio>/<run-label>/。",
    )
    parser.add_argument(
        "--transparent-optics",
        action="store_true",
        help="调试模式：把 mu_a=mu_s=g=0，用来检查厚度趋势是否来自光学参数。",
    )
    parser.add_argument(
        "--detected-photons",
        action="store_true",
        help="输出 detected_photons.csv。文件可能很大。",
    )
    parser.add_argument(
        "--overwrite-sources",
        action="store_true",
        help="即使 source/event 次生文件已存在，也重新预处理生成。",
    )
    parser.add_argument(
        "--keep-generated-sources",
        action="store_true",
        help="StageB 全流程运行后保留新生成的 source/event 中间 CSV。",
    )
    parser.add_argument(
        "--preprocess-only",
        action="store_true",
        help="只生成 source/event 次生文件，不运行 Monte Carlo。",
    )
    parser.add_argument(
        "--mc-only",
        action="store_true",
        help="跳过预处理，直接使用已有 source/event 次生文件运行 Monte Carlo。",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="只打印将要执行的命令，不实际运行。",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="批量运行后不自动生成厚度-出光量汇总图。",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="如果 OpticalMC.exe 不存在，不尝试自动编译。",
    )
    return parser.parse_args()


def opticalmc_executable_name() -> str:
    return "OpticalMC.exe" if os.name == "nt" else "OpticalMC"


def opticalmc_executable_candidates() -> List[Path]:
    name = opticalmc_executable_name()
    if os.name == "nt":
        return [
            REPO_ROOT / name,
            REPO_ROOT / "build" / "Release" / name,
            REPO_ROOT / "build" / name,
        ]
    return [
        REPO_ROOT / name,
        REPO_ROOT / "build" / name,
    ]


def ensure_executable(dry_run: bool, no_build: bool) -> Path:
    candidates = opticalmc_executable_candidates()
    for candidate in candidates:
        if candidate.exists():
            return candidate
    exe = candidates[0]
    if no_build:
        raise FileNotFoundError(f"{exe.name} 不存在；请先编译，或去掉 --no-build 让脚本自动编译。")
    cmd = [
        "g++",
        "-std=c++17",
        "-O3",
        "-pthread",
        str(REPO_ROOT / "src" / "OpticalMC.cpp"),
        "-o",
        str(exe),
    ]
    print(f"{exe.name} 不存在，先自动编译：", flush=True)
    print(" ".join(cmd), flush=True)
    if not dry_run:
        subprocess.run(cmd, check=True)
    return exe


def main() -> int:
    args = parse_args()
    preset = PRESETS[args.preset]

    try:
        opticalmc = ensure_executable(args.dry_run, args.no_build)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    xy_anchor_mode = "capture"
    cmd: List[str] = [
        sys.executable,
        str(REPO_ROOT / "run_opticalmc_batch.py"),
        "--ratio",
        args.ratio,
        "--opticalmc",
        str(opticalmc),
        "--samples-per-step",
        str(preset["samples_per_step"]),
        "--psf-bin-size-um",
        str(preset["psf_bin_size_um"]),
        "--psf-range-um",
        str(preset["psf_range_um"]),
        "--lsf-range-um",
        str(preset["lsf_range_um"]),
        "--num-threads",
        str(args.threads),
        "--xy-anchor-mode",
        xy_anchor_mode,
        "--readout-surface",
        args.readout,
        "--front-reflection-model",
        args.front_reflection_model,
        "--front-reflectance",
        str(args.front_reflectance),
        "--front-reflection-mode",
        args.front_reflection_mode,
        "--front-aluminum-n",
        str(args.front_aluminum_n),
        "--front-aluminum-k",
        str(args.front_aluminum_k),
        "--back-reflection-model",
        args.back_reflection_model,
        "--back-air-n",
        str(args.back_air_n),
        "--yield-zns-per-MeV",
        str(args.yield_zns_per_MeV),
        "--optical-component",
        args.optical_component,
        "--mu-a-scale",
        str(args.mu_a_scale),
        "--mu-s-scale",
        str(args.mu_s_scale),
        "--scattering-model",
        args.scattering_model,
        "--transport-scattering-mode",
        args.transport_scattering_mode,
        "--max-steps",
        str(args.max_steps),
    ]
    if args.incident_events is not None:
        cmd.extend(["--incident-event-count", str(args.incident_events)])
    if args.phase_function_csv is not None:
        cmd.extend(["--phase-function-csv", str(args.phase_function_csv)])
    if args.run_label:
        cmd.extend(["--run-label", args.run_label])

    if args.thickness:
        cmd.extend(["--thickness", *args.thickness])
    if args.quench_kb is not None:
        cmd.extend(["--quench-model", "birks", "--birks-kb-um-per-keV", str(args.quench_kb)])
    if args.detected_photons:
        cmd.append("--output-detected-photons")
    if args.transparent_optics:
        cmd.append("--transparent-optics")
    if args.overwrite_sources:
        cmd.append("--overwrite-sources")
    if args.keep_generated_sources:
        cmd.append("--keep-generated-sources")
    if args.preprocess_only:
        cmd.append("--preprocess-only")
    if args.mc_only:
        cmd.append("--mc-only")
    if args.dry_run:
        cmd.append("--dry-run")

    print(f"配比: {args.ratio}", flush=True)
    print(
        "厚度: " + (", ".join(args.thickness) if args.thickness else "全部可用厚度"),
        flush=True,
    )
    print(f"预设: {args.preset} ({preset['description']})", flush=True)
    print("坐标: capture/global", flush=True)
    print(f"读出面: {args.readout}", flush=True)
    print(f"光学参数: {args.optical_component}", flush=True)
    print(f"输运散射模式: {args.transport_scattering_mode}", flush=True)
    print(
        "入射 neutron history: "
        + (f"{args.incident_events:.12g}" if args.incident_events is not None else "1000000 (默认)"),
        flush=True,
    )
    print(f"最大步数: {args.max_steps}", flush=True)
    print("执行命令:", flush=True)
    print(" ".join(cmd), flush=True)

    if not args.dry_run:
        subprocess.run(cmd, check=True)
        if not args.no_plot and not args.preprocess_only:
            if args.run_label is None:
                plot_cmd = [sys.executable, str(REPO_ROOT / "plot_thickness_light.py"), args.ratio]
            else:
                run_root = Path("outputs") / args.ratio / args.run_label
                plot_cmd = [
                    sys.executable,
                    str(REPO_ROOT / "plot_thickness_light.py"),
                    "--summary",
                    str(run_root / "thickness_light_summary.csv"),
                    "--output-dir",
                    str(run_root / "figures"),
                ]
            print("生成厚度-出光量汇总图:", flush=True)
            print(" ".join(plot_cmd), flush=True)
            subprocess.run(plot_cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
