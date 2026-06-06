# Folder 01 — `codes/` Python translation

This folder contains the Python translation of the MATLAB files in `ZhuangLab/MERFISH_analysis/codes`.

The original repository describes itself as MATLAB software for MERFISH analysis and probe construction. This translated folder covers the barcode/codebook construction utilities only.

## Translated files

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `CodebookToMap.m` | `codes/CodebookToMap.py` | Reads a MERFISH codebook from FASTA-like records or file and builds a mapping from codeword to gene/object name; can include SECDED-correctable words. |
| `GenSECDED.m` | `codes/GenSECDED.py` | Generates SECDED codewords from binary data words and filters by number of on bits. |
| `GenerateExtendedHammingWords.m` | `codes/GenerateExtendedHammingWords.py` | Generates extended Hamming codewords and the generator matrix for a requested number of data bits. |
| `GenerateSurroundingCodewords.m` | `codes/GenerateSurroundingCodewords.py` | Generates all binary codewords at an exact Hamming distance from a central codeword. |
| `SECDEDCorrectableWords.m` | `codes/SECDEDCorrectableWords.py` | Returns all single-bit-neighbor words correctable to a SECDED codeword. |

## Usage

```python
import numpy as np
from codes import GenerateSurroundingCodewords, SECDEDCorrectableWords, GenSECDED

neighbors = GenerateSurroundingCodewords(np.array([1, 0, 1, 0], dtype=bool), 1)
correctable = SECDEDCorrectableWords(np.array([1, 0, 1, 0], dtype=bool))
words = GenSECDED(numLetters=16, numDataBits=11, onBits=4)
```

## Test

```bash
python tests/test_codes_folder.py
```
