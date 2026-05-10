# Folder 04 — `startup/` Python translation

This package translates the MERFISH_analysis `startup/` folder from MATLAB to Python.

## MATLAB source covered

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `startup/merfish_startup.m` | `startup/merfish_startup.py` | Initializes local paths and environment variables for MERFISH_analysis, matlab-storm, storm-analysis, legacy BLAST/OligoArray, and export_fig. |

## Use

```python
from startup import MerfishStartupConfig, configure_merfish_environment

cfg = MerfishStartupConfig(
    merfish_analysis_path="/path/to/MERFISH_analysis",
    matlab_storm_path="/path/to/matlab-storm",
    storm_analysis_path="/path/to/storm-analysis",
    scratch_path="/tmp/merfish_scratch",
)

report = configure_merfish_environment(cfg, strict=False, verbose=True)
```

## Translation notes

- MATLAB `addpath(path)` is translated to insertion into `sys.path`.
- MATLAB `addpath(genpath(path), '-begin')` is translated to recursive insertion of all subdirectories into `sys.path`.
- MATLAB `global` path variables are translated to `os.environ` variables plus a returned report dictionary.
- MATLAB `close all`, `clear all`, and `clc` have no safe Python runtime equivalent and are intentionally omitted.

## Checks

- Python compile: passed
- Test script: passed
- `matlab_stmt`: 0
- placeholders: 0
- `NotImplementedError`: 0
