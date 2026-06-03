# Mask Imaging

This directory contains an independent OpticalMC wrapper for source-position mask imaging.
It does not modify `src/OpticalMC.cpp`, `run_opticalmc_batch.py`, or the generated source
format used by the original workflow.

The mask image must be square and black/white. It is mapped onto neutron source positions
`source_x_um/source_y_um` in the StageB `<thickness>_capture_anchors.csv` file. By default,
the image covers `x = +/-5000 um` and `y = +/-5000 um`; white pixels pass, black interior
pixels block, and black edge pixels are treated as usable boundary events.

Example:

```powershell
python .\mask_imaging\run_mask_imaging.py `
  --ratio 1-1 `
  --thickness 100 `
  --mask-image .\masks\pattern.png `
  --samples-per-step 16 `
  --num-threads 8
```

Main outputs are written to:

```text
outputs/mask_imaging/<ratio>/<thickness>/
  detected_photon_positions.csv
  photon_hit_map.png
  simulated_radiograph.png
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
