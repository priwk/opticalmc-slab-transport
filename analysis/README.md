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
      thickness_detection_efficiency.png
      thickness_light_per_capture.png
      thickness_spot_spread.png
      thickness_photon_budget.png
      thickness_mtf_thresholds.png
      event_depth_detection_efficiency.png
      psf/
      lsf/
      mtf/
  ratio_comparison/
    ratio_comparison_plot_data.csv
    ratio_compare_detected_light.png
    ratio_compare_detection_efficiency.png
    ratio_compare_capture_fraction.png
    ratio_compare_fwhm.png
    ratio_comparison_summary.png
```
