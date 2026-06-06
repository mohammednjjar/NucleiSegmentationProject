"""Python translation of `codes/GenerateSurroundingCodewords.m`."""
from __future__ import annotations

from itertools import combinations
from typing import Iterable, List, Sequence

import numpy as np


def GenerateSurroundingCodewords(codeword: Sequence[bool] | np.ndarray,
                                 hammDist: int,
                                 logical: bool = False):
    """Generate all codewords exactly `hammDist` bits away from `codeword`.

    Parameters
    ----------
    codeword:
        One-dimensional boolean/0-1 vector.
    hammDist:
        Exact Hamming distance from the input codeword.
    logical:
        MATLAB option equivalent. If `True`, return a 2-D NumPy boolean
        array. If `False`, return a list of 1-D boolean arrays, matching the
        MATLAB cell-array output.

    Returns
    -------
    list[np.ndarray] | np.ndarray
        Surrounding codewords.
    """
    if codeword is None:
        raise ValueError("Both a codeword and a hamming distance are required.")
    arr = np.asarray(codeword)
    if arr.ndim != 1:
        raise ValueError("Codeword must be one-dimensional.")
    if arr.dtype != np.bool_:
        # MATLAB required logical input. For practical Python use, accept 0/1
        # and convert, but reject non-binary values.
        if not np.all((arr == 0) | (arr == 1)):
            raise TypeError("Codeword must be logical/binary.")
        arr = arr.astype(bool)
    if not isinstance(hammDist, (int, np.integer)):
        raise ValueError("Incorrect hamming distance.")
    hammDist = int(hammDist)
    if hammDist < 0 or hammDist > arr.size:
        raise ValueError("Incorrect hamming distance.")

    combos = list(combinations(range(arr.size), hammDist))
    out = np.tile(arr, (len(combos), 1))
    for i, idxs in enumerate(combos):
        out[i, list(idxs)] = ~out[i, list(idxs)]

    if logical:
        return out
    return [out[i, :].copy() for i in range(out.shape[0])]
