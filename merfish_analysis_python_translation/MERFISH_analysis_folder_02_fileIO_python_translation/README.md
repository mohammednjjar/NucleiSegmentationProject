# MERFISH_analysis folder 02: `fileIO/` Python translation

This package translates the repository folder `fileIO/` and `fileIO/deprecated/` from MATLAB to Python.

## Translated files

| MATLAB source | Python translation | Purpose/use |
|---|---|---|
| `fileIO/BuildFileStructure.m` | `fileIO/BuildFileStructure.py` | Parse file names into records using delimiters/regex and conversion functions. |
| `fileIO/CatBinaryFiles.m` | `fileIO/CatBinaryFiles.py` | Append one MERFISH custom binary file to another after header/layout validation. |
| `fileIO/LoadByteStream.m` | `fileIO/LoadByteStream.py` | Load a serialized `.matb` byte-stream object. Python equivalent uses `pickle`. |
| `fileIO/LoadCodebook.m` | `fileIO/LoadCodebook.py` | Parse MERFISH codebook CSV into header + barcode records. |
| `fileIO/LoadDataOrganization.m` | `fileIO/LoadDataOrganization.py` | Parse MERFISH data-organization CSV and return metadata. |
| `fileIO/LoadSplitByteStream.m` | `fileIO/LoadSplitByteStream.py` | Load split `.matb` byte-stream chunks and reshape. |
| `fileIO/LoadSplitByteStreamHeader.m` | `fileIO/LoadSplitByteStreamHeader.py` | Load and validate split-byte-stream metadata. |
| `fileIO/ReadBinaryFile.m` | `fileIO/ReadBinaryFile.py` | Read MERFISH custom binary format into Python dictionaries. |
| `fileIO/ReadBinaryFileHeader.m` | `fileIO/ReadBinaryFileHeader.py` | Read fixed and variable binary-file headers. |
| `fileIO/SaveAsByteStream.m` | `fileIO/SaveAsByteStream.py` | Save Python objects to `.matb`-style files using `pickle`. |
| `fileIO/SaveFigure.m` | `fileIO/SaveFigure.py` | Save Matplotlib figures in selected formats; `.fig` is stored as a pickled Matplotlib figure. |
| `fileIO/SaveSplitByteStream.m` | `fileIO/SaveSplitByteStream.py` | Save large array-like objects as split serialized byte-stream blocks. |
| `fileIO/SetFigureSavePath.m` | `fileIO/SetFigureSavePath.py` | Set/create global figure-save directory. |
| `fileIO/WriteBinaryFile.m` | `fileIO/WriteBinaryFile.py` | Write dictionaries to MERFISH custom binary format. |
| `fileIO/WriteCodebook.m` | `fileIO/WriteCodebook.py` | Write MERFISH codebook CSV from barcodes/readouts/names/ids. |
| `fileIO/deprecated/BuildImageDataStructures.m` | `fileIO/deprecated/BuildImageDataStructures.py` | Legacy wrapper that builds imageData records from parsed filenames. |
| `fileIO/deprecated/CreateImageDataStructure.m` | `fileIO/deprecated/CreateImageDataStructure.py` | Legacy empty imageData record generator. |
| `fileIO/deprecated/CreateWordsStructure.m` | `fileIO/deprecated/CreateWordsStructure.py` | Legacy empty decoded-word record generator. |

## Translation notes

- No `matlab_stmt(...)` wrappers are used.
- No placeholder method stubs are used.
- The MATLAB byte-stream functions use undocumented MATLAB serialization; the Python translation uses `pickle`, so it is Python-to-Python compatible, not MATLAB-bytestream-compatible.
- `WriteCodebook.m` in the source references undefined internal variables (`finalBarcodes`, `finalGenes`, `finalTargetRegions`). The Python version follows the documented function signature and uses `barcodes`, `readouts`, `names`, and `ids`.
- The deprecated image/word structures depended on `CreateMoleculeList` and MATLAB affine transform objects from external MATLAB-STORM code; the Python translation supplies NumPy-compatible equivalents.
