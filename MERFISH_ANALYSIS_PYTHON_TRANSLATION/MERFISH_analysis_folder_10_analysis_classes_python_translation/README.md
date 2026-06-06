# MERFISH_analysis folder 10: analysis/classes Python translation

This package translates the repository folder `analysis/classes/` into Python.

Included MATLAB-to-Python files:

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `FoundFeature.m` | `FoundFeature.py` | Segmented feature/cell container; boundaries, morphology, overlap, masks, table export |
| `MERFISHDecoder.m` | `MERFISHDecoder.py` | Main MERFISH dataset controller; metadata, preprocessing, decoding, segmentation, reporting, save/load |
| `Morpheus.m` | `Morpheus.py` | Email/status notification helper for jobs/job arrays |
| `SLURMJob.m` | `SLURMJob.py` | Single SLURM job wrapper; script creation, submission, status, cancel/requeue |
| `SLURMJobArray.m` | `SLURMJobArray.py` | Coordinates multiple SLURMJob objects and tracks array status |
| `Vrykolakas.m` | `Vrykolakas.py` | Timer that prints keep-alive messages for shell/SSH sessions |

Checks performed:

- Python compile check passed.
- Basic unit tests passed.
- No `matlab_stmt` wrappers.
- No `NotImplementedError` stubs.

Scientific/runtime note: this folder contains the largest MATLAB pipeline class (`MERFISHDecoder.m`). The Python translation implements portable behavior with numpy, pandas, scikit-image, scipy, tifffile, and shapely where available. Exact equivalence to the original MATLAB runtime still depends on the original raw MERFISH image formats and external STORM/MATLAB environment used by the Zhuang Lab pipeline.
