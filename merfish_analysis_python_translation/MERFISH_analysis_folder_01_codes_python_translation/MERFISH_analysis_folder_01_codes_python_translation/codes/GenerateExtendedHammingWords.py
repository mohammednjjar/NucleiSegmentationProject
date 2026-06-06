"""Python translation of `codes/GenerateExtendedHammingWords.m`."""
from __future__ import annotations

import numpy as np

try:
    from .hamming_utils import all_binary_messages, extended_hamming_generator
except ImportError:  # direct script execution
    from hamming_utils import all_binary_messages, extended_hamming_generator


def GenerateExtendedHammingWords(numDataBits: int,
                                 numOn: int = 0,
                                 parallel=None):
    """Return extended-Hamming codewords, generator matrix, and parity count.

    Parameters
    ----------
    numDataBits:
        Number of data bits.
    numOn:
        If nonzero, keep only codewords with exactly this number of 1 bits.
    parallel:
        Accepted for API compatibility with MATLAB. Computation is vectorized
        in NumPy, so this parameter is ignored.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, int]
        `(words, gen, numParityBits)`.
    """
    if numDataBits is None or int(numDataBits) < 1:
        raise ValueError("A valid number of data bits is required.")
    numDataBits = int(numDataBits)
    numOn = int(numOn or 0)
    if numOn < 0:
        raise ValueError("numOn must be non-negative.")

    gen, numParityBits = extended_hamming_generator(numDataBits)
    messages = all_binary_messages(numDataBits)
    words = np.mod(messages @ gen, 2).astype(np.uint8)

    if numOn != 0:
        words = words[words.sum(axis=1) == numOn]
    return words, gen.astype(np.uint8), int(numParityBits)
