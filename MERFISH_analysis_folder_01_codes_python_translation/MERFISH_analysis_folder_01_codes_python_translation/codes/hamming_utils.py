"""Binary Hamming/SECDED utilities used by the translated MERFISH `codes/` folder.

The original MATLAB code relies on Communications Toolbox functions such as
`encode`, `hammgen`, `gen2par`, and `de2bi`.  This module reimplements the
needed binary linear-code operations with NumPy only.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import ceil, log2
from typing import Iterable, List, Sequence, Tuple

import numpy as np


def _as_binary_array(x: Sequence[int] | np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=np.uint8)
    if arr.ndim != 1:
        raise ValueError("codeword must be one-dimensional")
    if not np.all((arr == 0) | (arr == 1)):
        raise ValueError("codeword must contain only 0/1 values")
    return arr


def int_to_bits(value: int, width: int, *, left_msb: bool = True) -> np.ndarray:
    """Return a 0/1 array representation of `value` with fixed `width` bits."""
    if value < 0:
        raise ValueError("value must be non-negative")
    if width < 1:
        raise ValueError("width must be positive")
    bits = [(value >> i) & 1 for i in range(width)]
    if left_msb:
        bits = bits[::-1]
    return np.asarray(bits, dtype=np.uint8)


def bits_to_int(bits: Sequence[int] | np.ndarray, *, left_msb: bool = True) -> int:
    """Convert a binary vector to an integer."""
    arr = _as_binary_array(bits)
    if not left_msb:
        arr = arr[::-1]
    out = 0
    for bit in arr:
        out = (out << 1) | int(bit)
    return int(out)


def all_binary_messages(num_data_bits: int) -> np.ndarray:
    """All data words ordered like MATLAB `0:2^k-1` with MSB first columns."""
    if num_data_bits < 1:
        raise ValueError("num_data_bits must be >= 1")
    return np.vstack([int_to_bits(i, num_data_bits, left_msb=True) for i in range(2 ** num_data_bits)])


def required_hamming_parity_bits(num_data_bits: int) -> int:
    """Smallest r satisfying 2^r >= k + r + 1."""
    if num_data_bits < 1:
        raise ValueError("num_data_bits must be >= 1")
    r = 1
    while 2**r < num_data_bits + r + 1:
        r += 1
    return r


def hamming_encode_messages(messages: np.ndarray, code_length: int | None = None) -> np.ndarray:
    """Encode binary messages using a systematic binary Hamming code.

    Parameters
    ----------
    messages:
        Matrix of shape `(n_messages, k)`, with one message per row.
    code_length:
        Optional Hamming codeword length. If omitted, uses `2^r - 1` for the
        smallest valid `r`. If provided, it must be at least `k + r`, where
        `r` is the number of parity bits required for `k`.

    Returns
    -------
    np.ndarray
        Matrix of shape `(n_messages, n)` containing Hamming codewords.

    Notes
    -----
    This is a pure-Python replacement for MATLAB's
    `encode(..., 'hamming/binary')`. It uses parity positions 1,2,4,8,...
    and data in the remaining positions.
    """
    msg = np.asarray(messages, dtype=np.uint8)
    if msg.ndim == 1:
        msg = msg.reshape(1, -1)
    if not np.all((msg == 0) | (msg == 1)):
        raise ValueError("messages must contain only 0/1 values")

    k = msg.shape[1]
    r = required_hamming_parity_bits(k)
    min_n = k + r
    if code_length is None:
        code_length = 2**r - 1
    if code_length < min_n:
        raise ValueError(f"code_length={code_length} is too small for k={k}; minimum is {min_n}")

    parity_positions = {2**i for i in range(r)}
    data_positions = [pos for pos in range(1, code_length + 1) if pos not in parity_positions]
    if len(data_positions) < k:
        raise ValueError("not enough non-parity positions for message bits")

    code = np.zeros((msg.shape[0], code_length), dtype=np.uint8)
    # Place data bits in non-parity positions, using the first k slots.
    for j, pos in enumerate(data_positions[:k]):
        code[:, pos - 1] = msg[:, j]

    # Fill parity bits so every parity-check group has even parity.
    for p in parity_positions:
        covered = [idx for idx in range(1, code_length + 1) if idx & p]
        non_parity_covered = [idx for idx in covered if idx != p]
        parity = np.mod(code[:, np.asarray(non_parity_covered) - 1].sum(axis=1), 2)
        code[:, p - 1] = parity
    return code


def secded_from_hamming(hamming_words: np.ndarray, *, parity_first: bool = True) -> np.ndarray:
    """Add one global parity bit to Hamming words to produce SECDED words."""
    hw = np.asarray(hamming_words, dtype=np.uint8)
    if hw.ndim == 1:
        hw = hw.reshape(1, -1)
    global_parity = np.mod(hw.sum(axis=1), 2).astype(np.uint8).reshape(-1, 1)
    if parity_first:
        return np.concatenate([global_parity, hw], axis=1)
    return np.concatenate([hw, global_parity], axis=1)


def rref_gf2(matrix: np.ndarray) -> Tuple[np.ndarray, List[int]]:
    """Reduced row-echelon form over GF(2), returning `(rref, pivot_columns)`."""
    a = np.asarray(matrix, dtype=np.uint8).copy() % 2
    rows, cols = a.shape
    pivot_cols: List[int] = []
    r = 0
    for c in range(cols):
        pivot = None
        for rr in range(r, rows):
            if a[rr, c]:
                pivot = rr
                break
        if pivot is None:
            continue
        if pivot != r:
            a[[r, pivot]] = a[[pivot, r]]
        for rr in range(rows):
            if rr != r and a[rr, c]:
                a[rr, :] ^= a[r, :]
        pivot_cols.append(c)
        r += 1
        if r == rows:
            break
    return a, pivot_cols


def parity_check_matrix_hamming(num_parity_bits: int) -> np.ndarray:
    """Standard Hamming parity-check matrix with columns 1..2^r-1, MSB rows."""
    if num_parity_bits < 2:
        raise ValueError("num_parity_bits must be >= 2")
    n = 2**num_parity_bits - 1
    cols = [int_to_bits(i, num_parity_bits, left_msb=True) for i in range(1, n + 1)]
    return np.asarray(cols, dtype=np.uint8).T


def extended_hamming_generator(num_data_bits: int) -> tuple[np.ndarray, int]:
    """Return a shortened extended-Hamming generator matrix.

    This follows the MATLAB algorithm conceptually:
    1. create a Hamming parity-check matrix,
    2. add an overall parity bit/row,
    3. reduce over GF(2),
    4. convert parity-check form to a generator matrix,
    5. shorten to the requested number of data bits.
    """
    r = required_hamming_parity_bits(num_data_bits)
    h = parity_check_matrix_hamming(max(ceil(log2(num_data_bits + r)), 3))
    # Add the extra global parity bit column and global parity row.
    h = np.vstack([np.hstack([h, np.zeros((h.shape[0], 1), dtype=np.uint8)]),
                   np.ones((1, h.shape[1] + 1), dtype=np.uint8)])
    h_rref, pivots = rref_gf2(h)
    pivot_set = set(pivots)
    free_cols = [c for c in range(h_rref.shape[1]) if c not in pivot_set]
    k_full = len(free_cols)
    n_full = h_rref.shape[1]
    gen = np.zeros((k_full, n_full), dtype=np.uint8)
    for row, free_col in enumerate(free_cols):
        gen[row, free_col] = 1
        for pivot_row, pivot_col in enumerate(pivots):
            gen[row, pivot_col] = h_rref[pivot_row, free_col]
    # Match MATLAB shortening step: remove equal number of trailing rows/cols.
    excess = gen.shape[0] - num_data_bits
    if excess > 0:
        gen = gen[:-excess, :-excess]
    return gen % 2, r
