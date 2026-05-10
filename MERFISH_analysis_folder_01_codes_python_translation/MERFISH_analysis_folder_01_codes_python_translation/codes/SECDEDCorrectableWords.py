"""Python translation of `codes/SECDEDCorrectableWords.m`."""
from __future__ import annotations

from typing import Sequence

import numpy as np

try:
    from .GenerateSurroundingCodewords import GenerateSurroundingCodewords
except ImportError:  # direct script execution
    from GenerateSurroundingCodewords import GenerateSurroundingCodewords


def SECDEDCorrectableWords(codeword: Sequence[bool] | np.ndarray):
    """Return all single-bit-error words corrected to `codeword` by SECDED.

    This is the MATLAB function's behavior: it calls
    `GenerateSurroundingCodewords(codeword, 1)`.
    """
    if codeword is None:
        raise ValueError("A codeword must be provided.")
    arr = np.asarray(codeword)
    if arr.dtype != np.bool_:
        if not np.all((arr == 0) | (arr == 1)):
            raise TypeError("Codeword must be logical/binary.")
        arr = arr.astype(bool)
    return GenerateSurroundingCodewords(arr, 1)
