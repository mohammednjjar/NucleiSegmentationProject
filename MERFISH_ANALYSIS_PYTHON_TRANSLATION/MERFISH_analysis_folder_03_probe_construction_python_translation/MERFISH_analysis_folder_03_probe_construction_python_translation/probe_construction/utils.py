"""Shared utilities for the Python translation of MERFISH_analysis/probe_construction.

The original MATLAB code relies on Bioinformatics Toolbox helpers such as nt2int,
int2nt, fastaread, fastawrite, and custom byte-stream files.  This module provides
pure-Python equivalents used by the translated classes.
"""
from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple
import math
import pickle
import re

import numpy as np

IUPAC_ORDER = "ACGTRYKMSWBDHVN-"
NT_TO_INT: Dict[str, int] = {base: i for i, base in enumerate(IUPAC_ORDER)}
NT_TO_INT.update({"U": 3})
INT_TO_NT: Dict[int, str] = {i: base for i, base in enumerate(IUPAC_ORDER)}
STRICT_NT_TO_INT: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_pickle(path: str | Path, obj: Any) -> None:
    path = Path(path)
    if path.parent:
        path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_pickle(path: str | Path) -> Any:
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def is_sequence_string(value: Any) -> bool:
    return isinstance(value, str) and len(value) > 0 and re.fullmatch(r"[A-Za-z\-\*]+", value) is not None


def sequence_to_int(seq: str | Sequence[int] | np.ndarray, acgt_only: bool = False) -> np.ndarray:
    """Convert nucleotide characters to 0-based integer encoding.

    A=0, C=1, G=2, T/U=3.  When acgt_only=False, IUPAC ambiguity codes are
    retained as values >3.  Unknown characters become -1.
    """
    if isinstance(seq, np.ndarray):
        return seq.astype(np.int16, copy=False)
    if isinstance(seq, (list, tuple)) and not isinstance(seq, str):
        return np.asarray(seq, dtype=np.int16)
    text = str(seq).upper().replace(" ", "").replace("\n", "").replace("\r", "")
    mapper = STRICT_NT_TO_INT if acgt_only else NT_TO_INT
    return np.asarray([mapper.get(ch, -1) for ch in text], dtype=np.int16)


def int_to_sequence(seq: Sequence[int] | np.ndarray) -> str:
    arr = np.asarray(seq, dtype=np.int16).ravel()
    return "".join(INT_TO_NT.get(int(v), "*") for v in arr)


def reverse_complement_int(seq: Sequence[int] | np.ndarray) -> np.ndarray:
    arr = np.asarray(seq, dtype=np.int16)
    out = np.full(arr.shape, -1, dtype=np.int16)
    valid = (arr >= 0) & (arr <= 3)
    out[valid] = 3 - arr[valid]
    return out[::-1]


def reverse_complement(seq: str | Sequence[int] | np.ndarray) -> str | np.ndarray:
    if isinstance(seq, str):
        return int_to_sequence(reverse_complement_int(sequence_to_int(seq, acgt_only=True)))
    return reverse_complement_int(seq)


def rolling_hash(seq: str | Sequence[int] | np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return base-4 hashes for every k-mer and a boolean validity mask."""
    arr = sequence_to_int(seq, acgt_only=False)
    if k <= 0:
        raise ValueError("k must be positive")
    if arr.size < k:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=bool)
    powers = (4 ** np.arange(k - 1, -1, -1, dtype=np.int64)).astype(np.int64)
    n = arr.size - k + 1
    hashes = np.empty(n, dtype=np.int64)
    valid = np.empty(n, dtype=bool)
    for i in range(n):
        window = arr[i : i + k]
        valid[i] = bool(np.all((window >= 0) & (window <= 3)))
        hashes[i] = int(np.sum(window.astype(np.int64) * powers)) if valid[i] else -1
    return hashes, valid


def sliding_sum(values: Sequence[float] | np.ndarray, window: int) -> np.ndarray:
    arr = np.asarray(values)
    if window <= 0:
        raise ValueError("window must be positive")
    if arr.size < window:
        return np.asarray([], dtype=float)
    return np.convolve(arr.astype(float), np.ones(window, dtype=float), mode="valid")


def sliding_mean(values: Sequence[float] | np.ndarray, window: int) -> np.ndarray:
    return sliding_sum(values, window) / float(window)


def fasta_read(path: str | Path) -> List[Dict[str, str]]:
    records: List[Dict[str, str]] = []
    header: str | None = None
    chunks: List[str] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"Header": header, "Sequence": "".join(chunks)})
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append({"Header": header, "Sequence": "".join(chunks)})
    return records


def fasta_write(path: str | Path, headers: Sequence[str], sequences: Sequence[str], append: bool = True, line_width: int = 80) -> None:
    path = Path(path)
    if path.parent:
        path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append and path.exists() else "w"
    with path.open(mode, encoding="utf-8") as handle:
        for header, seq in zip(headers, sequences):
            handle.write(f">{header}\n")
            seq = str(seq)
            for i in range(0, len(seq), line_width):
                handle.write(seq[i : i + line_width] + "\n")


def as_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return [value]


def flatten(list_of_lists: Iterable[Iterable[Any]]) -> List[Any]:
    out: List[Any] = []
    for x in list_of_lists:
        out.extend(list(x))
    return out


def to_jsonable(value: Any) -> Any:
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {k: to_jsonable(v) for k, v in value.items()}
    return value


def from_possible_fasta_struct(records: Any) -> List[str]:
    """Extract sequences from MATLAB fastaread-like structs, dicts, or lists."""
    if records is None:
        return []
    if hasattr(records, "intSequences"):
        return [np.asarray(s, dtype=np.int16) for s in records.intSequences]
    if isinstance(records, dict):
        if "Sequence" in records:
            seq = records["Sequence"]
            return seq if isinstance(seq, list) else [seq]
        return list(records.values())
    if isinstance(records, list):
        if len(records) == 0:
            return []
        if isinstance(records[0], dict) and "Sequence" in records[0]:
            return [r["Sequence"] for r in records]
        return records
    if isinstance(records, tuple):
        return list(records)
    raise TypeError("target sequences must be a Transcriptome, fasta-record list, dict, or sequence list")
