# MERFISH_analysis folder 03: `probe_construction/` Python translation

This package translates the MATLAB files in `ZhuangLab/MERFISH_analysis/probe_construction` into Python.

## Included MATLAB-to-Python files

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `OTMap.m` | `probe_construction/OTMap.py` | Array-backed key/value map for off-target seed penalties. |
| `OTMap2.m` | `probe_construction/OTMap2.py` | Dictionary-backed key/value map for faster lookup. |
| `OTTable.m` | `probe_construction/OTTable.py` | Builds n-mer off-target tables and calculates penalties. |
| `PrimerDesigner.m` | `probe_construction/PrimerDesigner.py` | Generates and filters orthogonal primers. |
| `TRDesigner.m` | `probe_construction/TRDesigner.py` | Designs target regions using GC, Tm, penalties, specificity, and tiling. |
| `TargetRegions.m` | `probe_construction/TargetRegions.py` | Stores/writes/saves target-region records. |
| `Transcriptome.m` | `probe_construction/Transcriptome.py` | Loads/indexes transcriptomes and sequence metadata. |

## Extra helper file

`probe_construction/utils.py` implements Python equivalents for MATLAB/Bioinformatics Toolbox helpers used by this folder:

- `nt2int` / `int2nt` style conversion
- reverse complement
- rolling k-mer hashing
- FASTA read/write
- pickle-based save/load replacing MATLAB byte streams

## Verification included

- `PYTHON_COMPILE_TEST.txt`
- `UNIT_TEST_RESULTS.txt`
- `PLACEHOLDER_SCAN.txt`
- `tests/test_probe_construction_folder.py`

## Basic usage

```python
from probe_construction import Transcriptome, OTTable, TRDesigner

tr = Transcriptome([
    ["id1"],
    ["geneA"],
    ["ACGTACGTACGT"],
    [1.0],
])

ot = OTTable(tr, seedLength=3)
designer = TRDesigner(
    transcriptome=tr,
    OTTables=[ot],
    OTTableNames=["self"],
    specificityTable=ot,
    verbose=False,
)
regions = designer.DesignTargetRegions(regionLength=4, GC=[0.0, 1.0], Tm=[-100, 100], specificity=[0, 10])
```

## Note

This folder has been translated as executable Python logic with tests. MATLAB-specific serialization is represented with Python pickle files because MATLAB byte-stream `.matb` files are not a Python-native format.
