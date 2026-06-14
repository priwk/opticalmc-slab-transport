# OpticalMC 后处理脚本

这些脚本只读取 `outputs/`，默认把表格和图片写到 `analysis_results/`。

最常用：

```powershell
python .\analysis\run_analysis.py 1-2
```

处理所有配比：

```powershell
python .\analysis\run_analysis.py
```

只为指定厚度绘制 PSF/LSF/MTF：

```powershell
python .\analysis\run_analysis.py 1-2 --thickness 30 50 100
```

对比不同配比的厚度趋势：

```powershell
python .\analysis\plot_ratio_comparison.py
```

只对比指定配比：

```powershell
python .\analysis\plot_ratio_comparison.py --ratio 1-1 1-2 1-3
```

默认把单位入射中子的读出光量按 `1000000` 个入射中子缩放：

```powershell
python .\analysis\plot_ratio_comparison.py --light-scale 1000000
```

这个脚本只读取 `outputs/<ratio>/thickness_light_summary.csv` 或已完成厚度目录，
默认输出到新的 `analysis_results/ratio_comparison` 目录；如果目录已存在，会自动使用
`ratio_comparison_1`、`ratio_comparison_2` 等新目录，避免覆盖已有对比结果。

输出目录：

```text
analysis_results/
  tables/
    summary_all.csv
    summary_<ratio>.csv
    mtf_metrics.csv
  figures/
    <ratio>/
      thickness_readout_light_per_million_incident.png
      thickness_detection_efficiency.png
      thickness_light_per_capture.png
      thickness_fwhm.png
      thickness_spot_spread.png
      thickness_photon_budget.png
      thickness_mtf_thresholds.png
      event_depth_detection_efficiency.png
      phase_function_mu.png
      phase_function_theta.png
      psf/
        psf_2d_thickness_panels.png
      lsf/
        lsf_x_thickness_panels.png
        lsf_y_thickness_panels.png
      mtf/
    phase_functions_mu_panels.png
    phase_functions_theta_panels.png
  ratio_comparison/
    ratio_comparison_plot_data.csv
    ratio_compare_readout_light_per_million_incident.png
    ratio_compare_readout_light_by_thickness_panels.png
    ratio_compare_detected_light.png
    ratio_compare_detection_efficiency.png
    ratio_compare_capture_fraction.png
    ratio_compare_fwhm.png
    ratio_compare_fwhm_by_thickness_panels.png
    ratio_compare_mtf50_mtf10.png
    ratio_compare_mtf50_by_thickness_panels.png
    ratio_compare_mtf10_by_thickness_panels.png
    ratio_compare_phase_functions_mu.png
    ratio_compare_phase_functions_theta.png
    ratio_compare_phase_functions_polar.png
    ratio_compare_phase_function_mu_panels.png
    ratio_compare_phase_function_theta_panels.png
    ratio_compare_phase_function_polar_panels.png
    ratio_comparison_summary.png
```

绘图会优先使用本机 `analysis/fonts/` 下的本地字体文件，例如 `arial.ttf` 和
`Deng.ttf`。这些字体文件不提交到 GitHub；需要时可从本机系统字体目录复制。坐标轴为
封闭框、刻度向内；六宫格图按最多 6 个代表厚度自动抽取，单个小图保持方形。
