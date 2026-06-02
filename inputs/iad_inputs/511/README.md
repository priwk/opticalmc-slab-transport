# IAD 输入数据：511 日期批次

本目录保存从原始 `inputs/511` 文本文件整理出的漫反射/漫透射光谱数据，
供后续 IAD 反演使用。这里的 `511` 是数据集/日期批次标签，不是默认反演波长。

样品数据已经做玻璃参考扣除：同一波长下用样品的漫反射/漫透射百分数分别减去
`Glass-R-D.txt` 和 `Glass-T-D.txt` 的百分数，然后再除以 100 得到 IAD 输入用的
`MR` 和 `MT`。

## 重新生成

如果本目录被清空，先确认原始文件还在 `inputs/511`，然后运行：

```powershell
python .\tools\prepare_iad_inputs.py
```

## 文件说明

- `samples.csv`：每个样品/厚度一行，包含扣玻璃后的 450 nm 和 511 nm 快速检查值。
- `measurements.csv`：所有样品光谱合并后的长表。
- `per_sample/*_iad_input.csv`：每个样品一个清洗后的光谱文件，包含
  `wavelength_nm`、`MR`、`MT` 等列，供 IAD 包装脚本读取。
- `references/glass_diffuse_reference.csv`：玻璃参考片的漫反射/漫透射光谱。
- `manifest.json`：机器可读的数据集摘要。
