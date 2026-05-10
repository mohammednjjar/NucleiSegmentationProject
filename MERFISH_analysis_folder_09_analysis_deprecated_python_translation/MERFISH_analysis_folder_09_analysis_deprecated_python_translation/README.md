# MERFISH_analysis Folder 09 Translation: `analysis/deprecated/`

This ZIP contains Python translations for the 7 MATLAB files currently listed in GitHub `analysis/deprecated/`.

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `AlignFiducials.m` | `python_translation/analysis/deprecated/AlignFiducials.py` | Align fiducial molecule lists across hybridization rounds and attach transforms/error metadata. |
| `AnalyzeMERFISH.m` | `python_translation/analysis/deprecated/AnalyzeMERFISH.py` | Legacy end-to-end MERFISH analysis orchestration: discover files, load molecule lists, align fiducials, create/decode words, and collect outputs. |
| `CreateWords.m` | `python_translation/analysis/deprecated/CreateWords.py` | Construct MERFISH word structures from aligned molecule lists using common-centroid or per-localization logic. |
| `DecodeWords.m` | `python_translation/analysis/deprecated/DecodeWords.py` | Decode measured codewords using exact and correctable codebook maps. |
| `StripWords.m` | `python_translation/analysis/deprecated/StripWords.py` | Reduce word structures to selected output fields. |
| `TransformImageData.m` | `python_translation/analysis/deprecated/TransformImageData.py` | Apply fiducial inverse transforms to image localizations and copy fiducial metadata. |
| `Warp2BestPair.m` | `python_translation/analysis/deprecated/Warp2BestPair.py` | Compute fiducial alignment transform from the most mutually consistent bead-pair shifts. |

Checks included in this ZIP: compile test, pytest tests, and source scan for scaffold markers.

External STORM/MATLAB file readers are replaced by Python-readable CSV/TSV/JSON/NPY/NPZ/Pickle readers plus injection hooks in `AnalyzeMERFISH.py`.
