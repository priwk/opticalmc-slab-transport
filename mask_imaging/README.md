# Mask Imaging 说明

`mask_imaging` 是一个独立的 OpticalMC 掩膜成像包装脚本。它不修改
`src/OpticalMC.cpp`、`run_opticalmc_batch.py`，也不改变原始生成源文件格式。

它的核心逻辑是：把一张正方形黑白掩膜图映射到 StageB 里的 neutron source
位置 `source_x_um/source_y_um`，然后根据该位置对应的掩膜像素给这个 neutron /
capture event 一个透过率权重。后续 OpticalMC 继续正常生成光子，最后成像时使用
weighted photon histogram。

## 推荐运行命令

在项目根目录运行：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --samples-per-step 32 \
  --num-threads 40
```

默认会遍历所有同时具有 StageB 输入和 optical 参数的配比，并输出四个厚度：

```text
50 100 200 500 um
```

如果只跑指定配比或厚度：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --ratio 1-1 1-3 \
  --thickness 50 100 200 500 \
  --samples-per-step 32 \
  --num-threads 40
```

## 掩膜透过率模型

当前掩膜不是“黑区完全删除 neutron 信息”的硬阻挡模型，而是参数化透过率模型。

默认参数：

```text
white_transmission = 1.0
black_transmission = 0.05
edge_transmission  = 0.2
```

也就是：

```text
白区        -> 100% 透过
黑区内部    -> 5% 透过
黑区边缘    -> 20% 透过
```

这里的 5% 透过不是“随机保留 5% 的 event”，而是保留 event，但把它的
`trajectory_weight` 乘以 `0.05`。这是一种连续权重近似，优点是统计噪声更小，
缺点是没有模拟真实掩膜材料里的散射、能量变化和角度偏转。

可以通过参数修改：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --black-transmission 0.05 \
  --edge-transmission 0.2 \
  --white-transmission 1.0
```

如果希望黑区边缘也按黑区内部处理：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --block-mask-edge
```

如果希望恢复“黑区完全吸收”的旧近似：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --black-transmission 0 \
  --edge-transmission 0
```

## 掩膜坐标映射

掩膜图必须是正方形。彩色图会先转换成亮度：

```text
luminance = 0.2126 R + 0.7152 G + 0.0722 B
```

默认阈值：

```text
--mask-threshold 0.5
```

因此：

```text
luminance < 0.5  -> 黑区
luminance >= 0.5 -> 白区
```

掩膜默认覆盖 StageB 源坐标：

```text
source_x_um: -5000 到 +5000 um
source_y_um: -5000 到 +5000 um
```

对应参数：

```text
--source-range-um 5000
```

注意：掩膜作用在 neutron source/capture anchor 的源位置上，不是直接作用在最终
readout 平面的 `readout_x_um/readout_y_um` 上。

如果 source 坐标落在掩膜范围外，默认允许通过：

```text
--outside-mask allow
```

如果希望范围外直接阻挡：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --outside-mask block
```

## 代码流程

每个厚度和配比的主要流程：

1. 读取 StageB 的 `<thickness>_capture_anchors.csv`。
2. 按 `source_x_um/source_y_um` 查找对应掩膜像素。
3. 得到该 source event 的 `mask_transmission`。
4. 如果 `mask_transmission <= 0`，这个 physical event 被删除。
5. 如果 `mask_transmission > 0`，这个 physical event 保留，并写入：

```text
trajectory_weight = 原 trajectory_weight * mask_transmission
mask_transmission = 当前掩膜透过率
```

6. 根据保留下来的 `physical_event_uid/source_event_uid` 同步过滤 StageB 相关文件。
7. 使用临时过滤后的 StageB 输入运行 OpticalMC。
8. 读取 detected photon，按 photon `weight` 累加成 radiograph histogram。

## 输出位置

默认输出根目录：

```text
outputs/mask_imaging/
```

注意：`masks/` 只是输入掩膜图片目录，不会作为默认输出目录。

单个配比和厚度的输出：

```text
outputs/mask_imaging/<ratio>/<thickness>/
  detected_photon_positions.csv
  photon_hit_map.png
  radiograph_histogram.npy
  simulated_radiograph.png
  simulated_radiograph_gamma.png
  simulated_radiograph_contrast.png
  image_metrics.csv
  mask_transmission_summary.csv
  mask_transmission_summary.json
```

每个配比的四厚度横向拼图：

```text
outputs/mask_imaging/<ratio>/
  radiograph_strip_50_100_200_500um.png
  radiograph_strip_50_100_200_500um_gamma.png
  radiograph_strip_50_100_200_500um_contrast.png
  radiograph_strip_50_100_200_500um_ratio_to_50um.png
  radiograph_strip_50_100_200_500um_relative_difference_to_50um.png
  metrics_50_100_200_500um.csv
```

自定义输出目录：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --output-dir mask_imaging_results \
  --samples-per-step 32 \
  --num-threads 40
```

## 应该看哪张图

`simulated_radiograph.png`

固定 `--display-max` 显示，不按每张图单独归一化。它用于查看同一显示尺度下的
线性亮度差，但不等于完整探测器响应图像。

`simulated_radiograph_gamma.png`

推荐优先看。它同样使用固定 `--display-max`，然后施加固定 `--display-gamma`
提亮暗部。因此它不会像局部归一化那样抹掉不同厚度之间的整体亮度差。

`simulated_radiograph_contrast.png`

局部对比增强图。它会按每张图自己的 1% 到 99% 分位数重新拉伸，因此适合看结构
细节，但不适合比较不同厚度的整体亮度。

四厚度横向条带同理：

```text
radiograph_strip_50_100_200_500um.png
radiograph_strip_50_100_200_500um_gamma.png
radiograph_strip_50_100_200_500um_contrast.png
```

其中推荐看：

```text
radiograph_strip_50_100_200_500um_gamma.png
```

## 显示参数

固定亮度上限：

```text
--display-max 5000
```

固定 gamma：

```text
--display-gamma 0.5
```

如果图仍然太暗，可以降低 gamma，例如：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --display-gamma 0.35
```

如果图整体过曝，可以提高 `--display-max`；如果图整体太暗，可以降低
`--display-max` 或降低 `--display-gamma`。

## 常用参数

```text
--source-range-um 5000     掩膜覆盖 source_x/source_y 的半宽范围。
--image-range-um 5000      输出图像覆盖 readout_x/readout_y 的半宽范围。
--raw-bin-size-um 25       photon_hit_map.png 的 bin 尺寸。
--image-pixels 512         simulated_radiograph.png 的像素宽高。
--blur-sigma-px 1.2        radiograph 显示图使用的 Gaussian blur。
--display-max 5000         PNG 固定显示上限；不做每图归一化。
--display-gamma 0.5        *_gamma.png 的固定 gamma；不做每图归一化。
--mask-threshold 0.5       黑白阈值。
--white-transmission 1.0   白区透过率。
--black-transmission 0.05  黑区内部透过率。
--edge-transmission 0.2    黑区边缘透过率。
--block-mask-edge          黑区边缘也使用 black-transmission。
--outside-mask block       阻挡落在掩膜坐标范围外的 source event。
```

## 中间文件

默认情况下，大的临时 source/event 文件和 OpticalMC 原始中间结果会在运行结束后删除。
如果需要调试临时 StageB 输入，可以加：

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --keep-intermediates
```
