# MERFISH_analysis Folder 07: `analysis/functions/` Python Translation

This ZIP contains the Python translation for the GitHub folder:

`MERFISH_analysis/analysis/functions/`

Translated MATLAB files:

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `GenerateGeoTransformReport.m` | `GenerateGeoTransformReport.py` | Generates geometric transform QC summaries and optional matplotlib figures from transform residuals. |
| `MERFISHPerformanceMetrics.m` | `MERFISHPerformanceMetrics.py` | Computes MERFISH barcode performance metrics: exact/corrected counts, bit-error counts, brightness-area histograms, barcode-density maps, optional FPKM correlation reports. |
| `MLists2Transform.m` | `MLists2Transform.py` | Builds geometric transforms from reference and moving molecule lists per frame and applies them to `xc`/`yc`. |

## Included helper

`geometric_transforms.py` implements Python equivalents for MATLAB `affine2d`, `fitgeotrans`, and `transformPointsForward` behavior used by this folder.

Supported transform types:

- `nonreflectivesimilarity`
- `similarity`
- `affine`
- `projective`
- `polynomial`

## Verification

- Python compile check: passed
- Unit tests: 3/3 passed
- `matlab_stmt`: 0
- `placeholder`: 0
- `TODO`: 0
- `NotImplementedError`: 0
- empty `pass` stubs: 0

## Important boundary

This ZIP translates only `analysis/functions/`. The rest of the large top-level `analysis/` folder remains to be handled separately:

- `analysis/classes/`
- `analysis/SLURM_scripts/`
- `analysis/deprecated/`

`MERFISHPerformanceMetrics.py` is fully implemented for Python-readable tables (`csv`, `tsv`, `parquet`, `json`, `npz`, `pickle`). If you want to run it on the original MATLAB custom binary `.bin` outputs, combine it with the translated `fileIO/` package or convert those binary barcode lists first.
