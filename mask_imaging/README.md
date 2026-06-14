# Mask Imaging

This directory contains an independent OpticalMC wrapper for source-position mask imaging.
It does not modify `src/OpticalMC.cpp`, `run_opticalmc_batch.py`, or the generated source
format used by the original workflow.

The mask image must be square and black/white. It is mapped onto neutron source positions
`source_x_um/source_y_um` in the StageB `<thickness>_capture_anchors.csv` file. By default,
the image covers `x = +/-5000 um` and `y = +/-5000 um`; white pixels pass, black interior
pixels block, and black edge pixels are treated as usable boundary events.

Single-run example:

```powershell
python .\mask_imaging\run_mask_imaging.py `
  --ratio 1-1 `
  --thickness 100 `
  --mask-image .\masks\pattern.png `
  --samples-per-step 16 `
  --num-threads 8
```

Batch example for every ratio with available StageB and optical-parameter inputs:

```powershell
python .\mask_imaging\run_mask_imaging.py `
  --mask-image .\masks\SYSU.png `
  --samples-per-step 16 `
  --num-threads 8
```

The batch default thicknesses are:

```text
50 100 200 500 um
```

You can override them:

```powershell
python .\mask_imaging\run_mask_imaging.py `
  --mask-image .\masks\SYSU.png `
  --ratio 1-1 1-3 `
  --thickness 50 100 200 500
```

Main outputs are written to:

```text
outputs/mask_imaging/<ratio>/<thickness>/
  detected_photon_positions.csv
  photon_hit_map.png
  simulated_radiograph.png

outputs/mask_imaging/<ratio>/
  radiograph_strip_50_100_200_500um.png
```

The strip image places the four simulated radiographs in one horizontal row. Each panel is
drawn as a square.

`masks/` is only the input-mask directory. It is not used for outputs unless you explicitly
pass `--output-dir`.

On Linux, for example:

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --samples-per-step 16 \
  --num-threads 8
```

The default output root is:

```text
outputs/mask_imaging/
```

Use this if you want a custom output location:

```bash
python mask_imaging/run_mask_imaging.py \
  --mask-image masks/SYSU.png \
  --output-dir mask_imaging_results \
  --samples-per-step 16 \
  --num-threads 8
```

Large intermediate source/event files and raw OpticalMC summary files are created in a
temporary directory and removed when the run finishes. Add `--keep-intermediates` only when
debugging.

Useful options:

```text
--source-range-um 5000     Mask coordinate half-range for source_x/source_y.
--image-range-um 5000      Output image half-range for readout_x/readout_y.
--raw-bin-size-um 25       Bin size for photon_hit_map.png.
--image-pixels 512         Pixel width/height for simulated_radiograph.png.
--blur-sigma-px 1.2        Gaussian blur applied to simulated_radiograph.png.
--block-mask-edge          Treat every black pixel as blocked, including edges.
--outside-mask block       Block source positions outside the mask coordinate range.
```
