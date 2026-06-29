# 独立 IAD 工作流

IAD 反演与当前 OpticalMC Monte Carlo 流程保持独立。
`tools/run_iad.py` 只读取整理后的漫反射/漫透射输入数据，并把 IAD 反演结果写到
`outputs/iad` 下；它不会调用 `run.py`、`run_opticalmc_batch.py` 或 `OpticalMC.exe`。

## IAD 可执行文件

Scott Prahl 的 IAD v3.16.3 源码放在：

```text
third_party/iad/iad-3.16.3
```

当前包装脚本默认使用这个本地可执行文件：

```text
third_party/iad/iad-3.16.3/src/iad.exe
```

该可执行文件由 v3.16.3 源码构建。上游项目使用 MIT 许可证，许可证文件保留在：

```text
third_party/iad/iad-3.16.3/License
```

## 默认 450 nm 反演

如果 `inputs/iad_inputs/511/samples.csv` 不存在，需要先从原始光谱重新生成 IAD 输入：

```powershell
python .\tools\prepare_iad_inputs.py
```

这一步要求原始文件仍在：

```text
inputs/511/
```

```powershell
python .\tools\run_iad.py
```

默认读取扣玻璃后的整理数据：

```text
inputs/iad_inputs/511/per_sample/*_iad_input.csv
```

默认输出：

```text
outputs/iad/511/iad_results_450nm.csv
```

默认参数：

- 波长：`450.0 nm`
- 各向异性因子：`g = 0.4`
- 样品折射率：不显式指定，使用 IAD 默认值；可用 `--sample-index` 覆盖
- Monte Carlo lost-light 修正：关闭，即传给 IAD 的参数为 `-M 0`

输出表同时包含 IAD 原生单位和项目中常用单位：

- `mua_1_per_mm`、`musp_1_per_mm`：IAD 原生 `1/mm`
- `mua_per_um`、`musp_per_um`、`mus_per_um`：换算后的 `1/um`
- `MR_raw_percent`、`glass_R_percent`、`MR_corrected_percent` 等列：记录玻璃扣除过程

这些输出不会被 OpticalMC 自动读取。

注意：`inputs/iad_inputs/511` 的 `MR` 和 `MT` 已经扣除了玻璃参考片。全谱数据里
有部分波长点扣玻璃后透射为负，450 nm 默认点是有效输入。

## 1.1 mm 玻璃基底样品手工批处理

如果新样品是附在 1.1 mm 玻璃上的单点数据，可以直接填写一个 CSV 表格：

```csv
配比,厚度,总透射,总反射,g
1-2,100,35.2,18.4,0.4
```

仓库内也有一个格式示例：

```text
tools/iad_glass_batch.example.csv
```

默认输入文件：

```powershell
python .\tools\run_iad_glass_batch.py --write-template inputs\iad_inputs\glass_1p1mm\samples.csv
```

填好数据后运行：

```powershell
python .\tools\run_iad_glass_batch.py
```

默认输出：

```text
outputs/iad/glass_1p1mm/iad_results.csv
```

默认假设：

- 样品厚度列单位为 `um`；如果输入的是 `mm`，运行时加 `--thickness-unit mm`，或把列名写成 `thickness_mm`/`厚度_mm`。
- `总透射` 和 `总反射` 自动识别单位：大于 1 的数按百分数处理，小于等于 1 的数按 0-1 分数处理。
- `配比` 按 BN:ZnS 质量比解析，例如 `1-2` 或 `1:2.5`。
- 样品折射率会按配比、密度和混合规则自动计算 `n_eff`，并传给 IAD 的 `-n`。
- 样品在一片 1.1 mm 底部玻璃上，即传给 IAD：`-G b -N 1.52 -D 1.1`。
- `bottom` 表示光先打到样品，再到玻璃。如果测量方向是光先过玻璃，再到样品，使用 `--slide-configuration top`。
- 玻璃折射率默认 `1.52`，可用 `--slide-index` 覆盖；玻璃厚度默认 `1.1 mm`，可用 `--slide-thickness-mm` 覆盖。

折射率默认使用与 `run_opticalmc_batch.py` 一致的参数：

```text
BN/ZnS 粉体 : PMMA : air = 64 : 21.6 : 14.4
n_BN = 1.80, n_ZnS = 2.36, n_PMMA = 1.49, n_air = 1.000293
rho_BN = 2.10 g/cm3, rho_ZnS = 4.09 g/cm3
混合规则：Lorentz-Lorenz
```

常用命令：

```powershell
python .\tools\run_iad_glass_batch.py --slide-configuration top
python .\tools\run_iad_glass_batch.py --n-model linear
python .\tools\run_iad_glass_batch.py --sample-index 1.6
python .\tools\run_iad_glass_batch.py --rt-unit percent
python .\tools\run_iad_glass_batch.py --dry-run
```

如果只测得总反射和总透射，并希望使用输运等效模型，固定
`g_eff = 0` 重新反演有效吸收和约化散射系数：

```powershell
python .\tools\run_iad_glass_batch.py --override-g 0 --quadrature 24 --extra-iad-args "-e 0.00015" --output-file iad_results_g0_full.csv --effective-output-file iad_results_g0_effective.csv
```

这会显式传给 IAD `-g 0`，不使用 `-j` 或 `-F` 固定散射系数。精简输出表使用：

```text
mua_eff_1_per_mm
musp_eff_1_per_mm
g_model
```

## 常用示例

只跑一个样品，使用默认 450 nm：

```powershell
python .\tools\run_iad.py --sample X-X-100
```

输出文件：

```text
outputs/iad/511/iad_results_X-X-100_450nm.csv
```

跑所有样品的 511 nm：

```powershell
python .\tools\run_iad.py --wavelength-nm 511
```

指定样品折射率：

```powershell
python .\tools\run_iad.py --sample-index 1.6
```

指定不同的 `g`：

```powershell
python .\tools\run_iad.py --g 0.75
```

跑完整光谱：

```powershell
python .\tools\run_iad.py --all-wavelengths
```

直接透传额外 IAD 参数：

```powershell
python .\tools\run_iad.py --extra-iad-args "-X -i 8"
```

只打印命令，不实际运行 IAD：

```powershell
python .\tools\run_iad.py --dry-run
```

显式指定输出 CSV：

```powershell
python .\tools\run_iad.py --output-file outputs\iad\511\custom_iad_results.csv
```
