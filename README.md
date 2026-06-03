# OpticalMC 使用说明

本项目分两步：

1. 用 Python 预处理 `alpha_li_steps`，生成宏观坐标下的 ZnS 闪烁光源。
2. 用 C++17 程序 `OpticalMC` 在均质 slab 中做光子 Monte Carlo 输运。

当前推荐目录结构：

```text
inputs/
  alpha_li_steps/
    1-2/
      30_alpha_li_steps.csv
      40_alpha_li_steps.csv
      ...
  optical_params/
    1-2/
      rve_raw_optical_params_by_ratio.csv
  generated_sources/
    1-2/
      30_macro_zns_step_sources.csv
      30_event_light_sources.csv
outputs/
  1-2/
    optical_properties.csv
    30/
      run_config.generated.json
      optical_mc_summary.csv
      optical_mc_event_summary.csv
      optical_mc_source_step_summary.csv
      psf_2d.csv
      lsf_x.csv
      lsf_y.csv
```

`inputs/alpha_li_steps` 只放原始数据。预处理生成的次生文件默认放到 `inputs/generated_sources`，不会污染原始数据目录。

## 编译 OpticalMC

如果已有 `OpticalMC.exe`，可以跳过这步。

用 g++ 编译：

```powershell
g++ -std=c++17 -O3 -pthread .\src\OpticalMC.cpp -o .\OpticalMC.exe
```

如果本机有 CMake，也可以：

```powershell
cmake -S . -B build
cmake --build build --config Release
```

## 最简单的批量运行方式

推荐用傻瓜入口 `run.py`。

跑某个配比的所有厚度：

```powershell
python .\run.py 1-2
```

这会自动运行 `inputs/alpha_li_steps/1-2` 下的所有厚度。

只跑指定厚度：

```powershell
python .\run.py 1-2 30 40 50 100
```

日常推荐：

```powershell
python .\run.py 1-2 30 40 50 70 100 150 200 300 --preset normal
```

先快速检查流程：

```powershell
python .\run.py 1-2 30 --preset quick
```

最终图表或更高统计：

```powershell
python .\run.py 1-2 30 40 50 70 100 150 200 300 --preset fine
```

`run.py` 会自动做这些事：

1. 读取 `inputs/alpha_li_steps/<ratio>/<thickness>_alpha_li_steps.csv`。
2. 读取 `inputs/optical_params/<ratio>/rve_raw_optical_params_by_ratio.csv` 中的 merged optical params。
3. 默认使用 bulk 主列估算宏观 OpticalMC 参数，并生成 `outputs/<ratio>/optical_properties.csv`。
4. 生成 ZnS source/event 文件到 `inputs/generated_sources/<ratio>`。
5. 为每个厚度生成 `outputs/<ratio>/<thickness>/run_config.generated.json`。
6. 调用 `OpticalMC.exe`。
7. 输出 MC 结果到 `outputs/<ratio>/<thickness>`。
8. 汇总所有已完成厚度到 `outputs/<ratio>/thickness_light_summary.csv`。

如果 `OpticalMC.exe` 不存在，`run.py` 会尝试用 g++ 自动编译。

三个预设：

```text
quick   samples_per_step = 4    适合确认流程
normal  samples_per_step = 32   适合日常批量
fine    samples_per_step = 128  适合最终结果，耗时更长
```

默认坐标模式是：

```text
--coords source-relative
```

这会使用 `corr_x_um/corr_y_um`，适合 PSF、LSF、MTF。如果要保留探测屏全局坐标，用：

```powershell
python .\run.py 1-2 30 --coords global
```

默认读出后表面：

```text
--readout back
```

可以改为：

```powershell
python .\run.py 1-2 30 --readout front
python .\run.py 1-2 30 --readout both
```

## 表面 Fresnel 反射

默认开启两个表面模型：

```text
front: n_eff -> Al(n + ik) 的 Fresnel 反射，镜面反射回闪烁体
back : n_eff -> air 的 Fresnel 反射，只有透射出后表面的光子才计入 detected/back_escape
```

默认前表面铝参数是 `front_aluminum_n = 0.65`、`front_aluminum_k = 5.3`，可按发光波长查表后覆盖：

```powershell
python .\run.py 1-2 30 --front-aluminum-n 0.65 --front-aluminum-k 5.3
```

后表面默认空气折射率是 `back_air_n = 1.000293`，可覆盖：

```powershell
python .\run.py 1-2 30 --back-air-n 1.000293
```

有效前表面模型保留常数反射率，用于灵敏度检查或独立实验标定：

```powershell
python .\run.py 1-2 30 40 50 --readout back --front-reflection-model effective --front-reflectance 0.85 --front-reflection-mode specular
```

`front_reflection_mode` 可选：

```text
none      保持旧的 front escape 边界
specular  镜面反射，反射时翻转 z 方向
diffuse   漫反射，反射时按 Lambertian 半球重新抽样方向
```

启用 `specular` 或 `diffuse` 后，光子打到前表面会按所选模型给出的概率反射回闪烁体；未反射的部分按铝板吸收计入 `absorbed_weight`，不再计入 `front_escape_weight`。

如果要回到旧边界行为，可显式关闭：

```powershell
python .\run.py 1-2 30 --front-reflection-mode none --back-reflection-model none
```

默认入射 neutron history 数是：

```text
100000
```

如果某批数据不是 100000 个入射 history，可以显式指定：

```powershell
python .\run.py 1-2 30 40 50 --incident-events 200000
```

## 高级批量入口

`run.py` 只是把常用参数包装得更简单。需要更细控制时，可以直接用底层批处理脚本：

```powershell
python .\run_opticalmc_batch.py --ratio 1-2 --thickness 30,40,50,100 --samples-per-step 32 --num-threads 8 --xy-anchor-mode corr --readout-surface back
```

## 配比和厚度

配比用目录名，例如：

```text
2-1
1-1
1-1.5
1-2
1-2.5
1-3
```

跑 `1-3` 所有厚度：

```powershell
python .\run.py 1-3
```

跑 `1-3` 的一批指定厚度：

```powershell
python .\run.py 1-3 30 40 50 60 70 80 90 100
```

注意：当前脚本需要对应配比下存在 RVE 光学参数文件：

```text
inputs/optical_params/<ratio>/rve_raw_optical_params_by_ratio.csv
```

新版 `inputs/optical_params/<ratio>/` 可以直接放 StageD 的四件套交付文件：

```text
monte_carlo_recommended_inputs.csv
monte_carlo_recommended_inputs.json
phase_function_mean_by_ratio.csv
rve_raw_optical_params_by_ratio.csv
```

脚本会优先读取 `monte_carlo_recommended_inputs.csv`，如果没有再读 JSON 或旧的
`rve_raw_optical_params_by_ratio.csv`。推荐输入的映射是：

```text
recommended_mu_a_per_um       -> mu_a_per_um
recommended_mu_s_per_um       -> mu_s_per_um
recommended_g1                -> g
recommended_mu_s_prime_per_um -> mu_s_prime_per_um
phase_function_file           -> phase_function_csv
```

如果回退到旧 merged 表，则读取：

```text
mu_a_expected_mean_per_um -> mu_a_per_um
mu_s_mean_per_um          -> mu_s_per_um
g_mean                    -> g
mu_s_prime_mean_per_um    -> mu_s_prime_per_um
```

逐散射 OpticalMC 使用 `mu_a + mu_s` 抽样自由程；`mu_s_prime` 只作为输运/扩散诊断输出，不会作为额外散射过程。默认 `--scattering-model auto` 会按推荐输入使用 `phase_function_mean_by_ratio.csv`；若要强制 HG(g)，加：

```powershell
python .\run.py 1-2 30 40 50 --scattering-model hg
```

`total_*` 和 `boundary_*` 列默认不会进入宏观 OpticalMC。若需要对照，可以显式指定：

```powershell
python .\run.py 1-2 30 40 50 --optical-component total
python .\run.py 1-2 30 40 50 --optical-component boundary
```

如果某个配比还没有这个文件，可以先用已有的 OpticalMC 格式参数表：

```powershell
python .\run_opticalmc_batch.py --ratio 1-2 --optical-properties .\some_optical_properties.csv
```

如果手动提供 tabulated phase function，MC 会按 `cos_theta_min/cos_theta_max` 和
`probability` 或 `probability_mean` 采样散射角余弦：

```powershell
python .\run.py 1-2 30 40 50 --phase-function-csv .\path\to\phase_function.csv
```

## 控制光子数

有两个不同层面的“光子数”。

物理光产额由 ZnS 光产额控制：

```powershell
--yield-zns-per-MeV 60000
```

计算方式：

```text
n_photon_step = yield_zns_per_MeV * visible_edep_keV / 1000
```

Monte Carlo 代表光子数由每个 step 的采样数控制。普通使用中直接选预设即可：

```powershell
python .\run.py 1-2 30 --preset quick
python .\run.py 1-2 30 --preset normal
python .\run.py 1-2 30 --preset fine
```

`samples_per_step` 越大，统计噪声越低，运行越慢。每个代表光子的权重是：

```text
n_photon_step / samples_per_step
```

## x/y 锚点模式

预处理时可选：

```powershell
--xy-anchor-mode corr
--xy-anchor-mode capture
```

推荐：

```powershell
--xy-anchor-mode corr
```

`corr` 会用 `corr_x_um/corr_y_um` 作为宏观 x/y 锚点。这里的 `corr` 不是任意校正后的绝对坐标，而是原 Geant4 项目中定义的“俘获点相对源点的横向偏移”：

```text
corr_x_um = capture_x_um - source_x_um
corr_y_um = capture_y_um - source_y_um
```

因此 `corr` 模式得到的是以入射源点为参考原点的横向坐标，适合 PSF、LSF、MTF 分析：它保留“中子从源点出发后，在屏内哪里俘获并发光”的相对位移。

`capture` 会用 `capture_x_um/capture_y_um` 作为宏观 x/y 锚点，适合保留探测屏全局坐标中的真实俘获位置分布。

由于当前 slab 的 x/y 边界默认无限，`corr` 和 `capture` 对总出光效率影响通常不大，主要影响光斑、PSF、LSF、MTF 这类空间分布结果。

## 淬灭接口

默认不淬灭：

```powershell
--quench-model none
```

启用简单 Birks step 级模型：

```powershell
python .\run.py 1-2 30 --quench-kb 0.001
```

当前公式：

```text
visible_edep_keV = edep_keV / (1 + kB * edep_keV / step_len_um)
```

第一版中 alpha 和 Li7 使用同一个 `kB`。后续可以扩展成按粒子类型或能量区间的淬灭表。

## n_eff

批处理脚本会临时用公式估算 `n_eff`。

总体体积分数默认：

```text
BN/ZnS 粉体 : PMMA : air = 64 : 21.6 : 14.4
```

BN/ZnS 粉体内部再按指定 BN:ZnS 质量比拆分，并用密度转成体积分数。默认参数：

```text
n_BN   = 1.80
n_ZnS  = 2.36
n_PMMA = 1.49
n_air  = 1.000293

rho_BN  = 2.10 g/cm3
rho_ZnS = 4.09 g/cm3
```

默认混合规则是 Lorentz-Lorenz：

```powershell
--n-model lorentz-lorenz
```

也可以用线性混合：

```powershell
--n-model linear
```

覆盖材料参数示例：

```powershell
python .\run_opticalmc_batch.py --ratio 1-2 --n-bn 1.9 --n-zns 2.4 --rho-bn 2.1 --rho-zns 4.09
```

说明：当前 OpticalMC 内部随机游走主要受 `mu_a`、`mu_s` 和散射相函数控制。默认散射相函数是 HG(g)；如果 `optical_properties.csv` 含有 `phase_function_csv`，或运行时传入 `--phase-function-csv`，则按表格中的 `mu = cos(theta)` 概率质量采样。程序默认用 `front_reflection_model`/`front_reflection_mode` 模拟前表面铝板反射，并用 `back_reflection_model = air_fresnel` 模拟后表面 `n_eff -> air` 的 Fresnel 反射；`detected` 表示实际透射出所选读出表面的光子。

## 光学参数敏感性检查

如果厚度-读出光量趋势明显不符合预期，建议先做两个对照。

透明介质对照：

```powershell
python .\run.py 1-2 30 40 50 70 100 --preset quick --transparent-optics
```

这个模式设置：

```text
mu_a = 0
mu_s = 0
g = 0
```

如果透明介质下 `mean_detected_light_per_incident` 随厚度增加，而真实参数下下降，说明 source 生成和厚度归一大体没反，下降主要来自光学吸收/散射或边界假设。

缩放 merged optical params 中选定的光学参数：

```powershell
python .\run.py 1-2 30 40 50 70 100 --preset quick --mu-a-scale 0.1 --mu-s-scale 0.1
```

这会把选定 optical component 的 `mu_a` 和 `mu_s` 同时乘以 0.1，用于判断趋势对光学参数强弱的敏感性。默认 component 是 `bulk`。

## 单独运行预处理器

一般不需要手动跑，但可以这样用：

```powershell
python .\make_macro_zns_sources.py .\inputs\alpha_li_steps\1-2\30_alpha_li_steps.csv --xy-anchor-mode corr
```

默认输出到：

```text
inputs/generated_sources/1-2/
```

也可以指定输出目录：

```powershell
python .\make_macro_zns_sources.py .\inputs\alpha_li_steps\1-2\30_alpha_li_steps.csv --output-dir .\inputs\generated_sources\1-2
```

## 单独运行 OpticalMC

批处理脚本会自动生成 `run_config.generated.json`。如果要手动运行：

```powershell
.\OpticalMC.exe .\outputs\1-2\30\run_config.generated.json
```

也可以用命令行覆盖路径：

```powershell
.\OpticalMC.exe .\run_config.example.json .\inputs\generated_sources\1-2\30_macro_zns_step_sources.csv .\inputs\generated_sources\1-2\30_event_light_sources.csv .\outputs\1-2\optical_properties.csv .\outputs\1-2\30
```

## 输出文件

每个厚度的输出目录中包括：

```text
optical_mc_summary.csv
optical_mc_event_summary.csv
optical_mc_source_step_summary.csv
psf_2d.csv
lsf_x.csv
lsf_y.csv
```

如果启用：

```powershell
--output-detected-photons
```

还会输出：

```text
detected_photons.csv
```

`optical_mc_summary.csv` 是厚度级总结果，重点看：

```text
mean_light_per_incident
mean_detected_light_per_incident
mean_light_per_capture
mean_detected_light_per_capture
detection_efficiency
capture_fraction
spot_rms_x
spot_rms_y
spot_rms_r
fwhm_x
fwhm_y
```

批量运行结束后，还会生成跨厚度汇总表：

```text
outputs/<ratio>/thickness_light_summary.csv
```

这张表就是“厚度-出光量/效率/光斑”的主表。重点列包括：

```text
thickness_um
mean_light_per_incident
mean_detected_light_per_incident
mean_light_per_capture
mean_detected_light_per_capture
detection_efficiency
capture_fraction
spot_rms_r
fwhm_x
fwhm_y
```

注意区分两个归一方式：

```text
per incident neutron：按入射中子数归一，适合画“厚度-总出光量/读出光量”。
per capture：只在已发生俘获的 event 上平均，适合看“每次俘获后的平均光学读出效率”。
```

如果研究屏厚优化，主图通常优先看：

```text
mean_detected_light_per_incident
```

而不是：

```text
mean_detected_light_per_capture
```

后者随厚度增加可能单调下降，因为厚屏中光子平均要走更长路径才到读出面。

## 厚度-出光量绘图

默认用 `run.py` 跑完 Monte Carlo 后，会自动调用：

```powershell
python .\plot_thickness_light.py <ratio>
```

例如：

```powershell
python .\plot_thickness_light.py 1-2
```

输入表：

```text
outputs/1-2/thickness_light_summary.csv
```

输出图：

```text
outputs/1-2/figures/thickness_light_curve.png
outputs/1-2/figures/thickness_light_curve.pdf
outputs/1-2/figures/thickness_detection_efficiency.png
outputs/1-2/figures/thickness_detection_efficiency.pdf
outputs/1-2/figures/thickness_mtf_thresholds.png
outputs/1-2/figures/thickness_mtf_thresholds.pdf
outputs/1-2/figures/paper_thickness_summary.png
outputs/1-2/figures/paper_thickness_summary.pdf
outputs/1-2/figures/thickness_plot_data.csv
outputs/1-2/figures/thickness_mtf_metrics.csv
```

其中最适合论文初稿的是：

```text
paper_thickness_summary.pdf
```

它是一个 2x2 多面板图：

```text
(a) detected photons per incident neutron vs thickness
(b) detection efficiency vs thickness
(c) FWHM/RMS spot size vs thickness
(d) normalized light-resolution trade-off
```

如果只想生成 PNG：

```powershell
python .\plot_thickness_light.py 1-2 --formats png
```

如果不想让 `run.py` 自动画图：

```powershell
python .\run.py 1-2 30 40 50 --no-plot
```

## 论文里通常怎么汇总

论文里一般不会“每个厚度一张独立图”作为主要结果，因为读者很难横向比较。更常见的是：

1. 主文放一张跨厚度曲线图：横轴厚度，纵轴出光量、探测效率或空间分辨率。
2. 用多面板图把几个核心指标放一起，例如出光量、效率、FWHM、RMS。
3. PSF/LSF 只选 3 到 5 个代表厚度展示，例如薄、中、厚、最佳厚度。
4. 完整每厚度 PSF/LSF 可以放补充材料，或者做成 heatmap/stacked curves。
5. 如果存在最佳厚度，主文通常强调 trade-off：出光量随厚度增加，但光斑扩散也可能变差。

因此建议主文图表结构：

```text
Figure 1: thickness_light_summary 的多面板汇总图
Figure 2: representative PSF/LSF for selected thicknesses
Figure 3: MTF curves or MTF50/MTF10 vs thickness
Supplementary: all thickness PSF/LSF maps
```

`optical_mc_event_summary.csv` 保留所有 event，包括 `total_n_photon = 0` 的零发光事件，因此平均光产额不会被偏高估计。

## 后处理和绘图

后处理脚本放在：

```text
analysis/
```

绘图和汇总结果默认写到：

```text
analysis_results/
```

不会写入 `inputs`，也不会把图片混到 `outputs` 里。

对某个配比画图：

```powershell
python .\analysis\run_analysis.py 1-2
```

对所有已经跑完的配比画图：

```powershell
python .\analysis\run_analysis.py
```

只给指定厚度生成 PSF/LSF/MTF 图：

```powershell
python .\analysis\run_analysis.py 1-2 --thickness 30 50 100
```

跨配比对比厚度趋势：

```powershell
python .\analysis\plot_ratio_comparison.py
```

只对比指定配比：

```powershell
python .\analysis\plot_ratio_comparison.py --ratio 1-1 1-2 1-3
```

跨配比脚本只读取项目里已有的 `outputs/<ratio>/thickness_light_summary.csv`
或已完成厚度目录，默认输出到 `analysis_results/ratio_comparison`。如果该目录已存在，
会自动改用 `ratio_comparison_1`、`ratio_comparison_2` 等新目录，避免覆盖已有结果。

主要输出：

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

## 独立 IAD 反演

IAD 反演工具与当前 OpticalMC 运行流程保持独立，不会被 `run.py` 或
`run_opticalmc_batch.py` 自动调用。

当前 511 日期批次的漫反射/漫透射数据已整理到，样品数据已扣除玻璃参考片：

```text
inputs/iad_inputs/511/
```

默认跑 450 nm 的所有厚度：

```powershell
python .\tools\run_iad.py
```

输出写到：

```text
outputs/iad/511/iad_results_450nm.csv
```

更多独立 IAD 用法见：

```text
tools/README_IAD.md
```
