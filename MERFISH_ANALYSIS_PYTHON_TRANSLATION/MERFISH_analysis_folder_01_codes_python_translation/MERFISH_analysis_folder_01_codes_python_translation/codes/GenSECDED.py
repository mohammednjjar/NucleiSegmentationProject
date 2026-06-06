"""Python translation of `codes/GenSECDED.m`."""
from __future__ import annotations

import numpy as np

try:
    from .hamming_utils import all_binary_messages, hamming_encode_messages, secded_from_hamming
except ImportError:  # direct script execution
    from hamming_utils import all_binary_messages, hamming_encode_messages, secded_from_hamming


def _generate_secded(num_letters: int, num_data_bits: int) -> np.ndarray:
    messages = all_binary_messages(num_data_bits)
    hamming_words = hamming_encode_messages(messages, code_length=num_letters - 1)
    return secded_from_hamming(hamming_words, parity_first=True)


def GenSECDED(numLetters: int, numDataBits: int, onBits=None):
    """Generate SECDED codewords.

    Parameters
    ----------
    numLetters:
        Total SECDED codeword length.
    numDataBits:
        Number of data bits required.
    onBits:
        If provided/non-empty, keep only codewords with this exact number of
        1 bits. If `None`, remove very short or very long words exactly as the
        MATLAB code does.

    Returns
    -------
    np.ndarray
        Matrix of selected SECDED codewords.
    """
    numLetters = int(numLetters)
    numDataBits = int(numDataBits)
    if numLetters < numDataBits:
        raise ValueError("Data bits cannot be more than total bits")
    if onBits is not None and numDataBits < int(onBits):
        raise ValueError("Used bits cannot be more than data bits")

    # MATLAB has a special hard-coded branch for `numLetters == 12` that first
    # builds 16-bit SECDED words from 11 data bits, then takes rows 1:16:end and
    # the first 12 columns.
    if numLetters == 12:
        secded16 = _generate_secded(16, 11)
        twelve_bit_secded = secded16[0::16, :12]
        if onBits is None:
            return twelve_bit_secded
        return twelve_bit_secded[twelve_bit_secded.sum(axis=1) == int(onBits)]

    secded = _generate_secded(numLetters, numDataBits)
    if onBits is not None:
        return secded[secded.sum(axis=1) == int(onBits)]

    mask = ~((secded.sum(axis=1) <= 4) | (secded.sum(axis=1) >= numLetters - 4))
    return secded[mask]
