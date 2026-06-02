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
