# MERFISH_analysis Folder 08 Translation: `analysis/SLURM_scripts/`

This ZIP contains Python translations for the 16 MATLAB files in `analysis/SLURM_scripts/`.

## Files translated

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `CalculateDoubletScore.m` | `CalculateDoubletScore.py` | Calculate per-feature barcode counts, centers of mass, and coordinate variances used for doublet/segmentation-error scoring. |
| `CalculateNumbers.m` | `CalculateNumbers.py` | Calculate barcode counts per feature and generate the feature-count report. |
| `CombineBoundaries.m` | `CombineBoundaries.py` | Combine segmentation boundaries/features after FOV-level segmentation. |
| `CombineFoundFeatures.m` | `CombineFoundFeatures.py` | Combine found features, generate a report, and export feature CSV. |
| `CombineSum.m` | `CombineSum.py` | Combine raw-signal summation outputs and generate a summation report. |
| `DecodeFOV.m` | `DecodeFOV.py` | Decode one MERFISH field of view. |
| `ExportBarcodeMetadata.m` | `ExportBarcodeMetadata.py` | Export barcode metadata and compute doublet-score metadata. |
| `LowResMosaic.m` | `LowResMosaic.py` | Generate low-resolution MERFISH mosaic outputs. |
| `MERFISHScheduler.m` | `MERFISHScheduler.py` | Create/coordinate SLURM workflow scripts for the MERFISH analysis pipeline. |
| `Optimize.m` | `Optimize.py` | Generate warp report, initialize scale factors, and optimize decoding scale factors. |
| `ParseFOV.m` | `ParseFOV.py` | Parse decoded barcodes into feature boundaries for one FOV. |
| `Performance.m` | `Performance.py` | Validate barcode files and run MERFISH performance metrics. |
| `ProcessFOV.m` | `ProcessFOV.py` | Warp and preprocess one FOV. |
| `Segment.m` | `Segment.py` | Segment all FOVs. |
| `SegmentFOV.m` | `SegmentFOV.py` | Segment one FOV. |
| `SumFOV.m` | `SumFOV.py` | Sum raw signal inside features for one FOV. |

## Notes

These MATLAB files are SLURM driver scripts. Their main behavior is loading a `MERFISHDecoder` object and calling decoder methods. The Python versions preserve that behavior using `decoder_loader=...` or an importable translated `MERFISHDecoder` class.

`CalculateDoubletScore.py` includes direct pandas/numpy implementation for the barcode count, center-of-mass, and variance CSV outputs.

`MERFISHScheduler.py` writes Python-based `.slurm` scripts using the translated driver modules.

## Checks

- Python compile check: success
- Unit tests: 5/5 success
- `matlab_stmt`: 0
- `NotImplementedError`: 0
- Python `pass` statements: 0
