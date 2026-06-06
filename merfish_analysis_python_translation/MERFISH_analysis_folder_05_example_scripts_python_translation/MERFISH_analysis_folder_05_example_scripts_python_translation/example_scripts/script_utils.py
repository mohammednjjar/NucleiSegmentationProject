"""Shared utilities for translated MERFISH example scripts."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence
import importlib
import math
import os
import pickle
import random
import shutil
import sys
import time

import numpy as np


@dataclass
class FastaRecord:
    Header: str
    Sequence: str

    @property
    def header(self) -> str:
        return self.Header

    @property
    def sequence(self) -> str:
        return self.Sequence


def page_break() -> None:
    print("-" * 66)


def tic() -> float:
    return time.perf_counter()


def toc(start: float) -> float:
    return time.perf_counter() - float(start)


def ensure_dir(path: str | os.PathLike[str]) -> Path:
    out = Path(path).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    return out


def set_figure_save_path(path: str | os.PathLike[str], make_dir: bool = True) -> str:
    out = Path(path).expanduser()
    if make_dir:
        out.mkdir(parents=True, exist_ok=True)
    return str(out) + ("/" if str(out)[-1:] not in ("/", "\\") else "")


def read_fasta(path: str | os.PathLike[str]) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    seq_parts: list[str] = []
    with Path(path).expanduser().open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.replace(" ", ""))
    if header is not None:
        records.append(FastaRecord(header, "".join(seq_parts)))
    return records


def write_fasta(path: str | os.PathLike[str], records: Sequence[FastaRecord | dict[str, str]]) -> None:
    out = Path(path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for rec in records:
            header = getattr(rec, "Header", None) if not isinstance(rec, dict) else rec.get("Header", "")
            seq = getattr(rec, "Sequence", None) if not isinstance(rec, dict) else rec.get("Sequence", "")
            handle.write(f">{header}\n")
            sequence = str(seq)
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return str(seq).translate(table)[::-1]


def strip_spaces(seq: str) -> str:
    return "".join(str(seq).split())


def load_bytestream(path: str | os.PathLike[str]) -> Any:
    with Path(path).expanduser().open("rb") as handle:
        return pickle.load(handle)


def save_bytestream(path: str | os.PathLike[str], obj: Any) -> None:
    out = Path(path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_mat_file(path: str | os.PathLike[str]) -> dict[str, Any]:
    try:
        from scipy.io import loadmat
    except Exception as exc:
        raise RuntimeError("scipy is required to read MATLAB .mat files") from exc
    return loadmat(str(Path(path).expanduser()), squeeze_me=True, struct_as_record=False)


def save_mat_file(path: str | os.PathLike[str], payload: dict[str, Any]) -> None:
    try:
        from scipy.io import savemat
    except Exception as exc:
        raise RuntimeError("scipy is required to write MATLAB .mat files") from exc
    out = Path(path).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    savemat(str(out), payload)


def import_symbol(candidates: Sequence[str]) -> Any:
    errors: list[str] = []
    for qualified_name in candidates:
        module_name, symbol_name = qualified_name.rsplit(".", 1)
        try:
            module = importlib.import_module(module_name)
            return getattr(module, symbol_name)
        except Exception as exc:
            errors.append(f"{qualified_name}: {exc}")
    raise RuntimeError("Required translated dependency is unavailable. Tried: " + "; ".join(errors))


def pairwise_hamming_matrix(words: np.ndarray) -> np.ndarray:
    arr = np.asarray(words, dtype=np.int16)
    n = arr.shape[0]
    distances = np.full((n, n), np.inf, dtype=float)
    for i in range(n):
        distances[i, :] = np.sum(np.abs(arr[i] - arr), axis=1)
        distances[i, i] = np.inf
    return distances


def combinations_to_binary(num_bits: int, on_bits: int) -> np.ndarray:
    import itertools

    combos = list(itertools.combinations(range(num_bits), on_bits))
    words = np.zeros((len(combos), num_bits), dtype=np.uint8)
    for row, combo in enumerate(combos):
        words[row, list(combo)] = 1
    return words


def select_random(records: Sequence[Any], k: int) -> list[Any]:
    k = min(int(k), len(records))
    if k <= 0:
        return []
    inds = list(range(len(records)))
    random.shuffle(inds)
    return [records[i] for i in inds[:k]]


def copy_script_to(src_file: str | os.PathLike[str], dst_dir: str | os.PathLike[str]) -> None:
    src = Path(src_file).expanduser()
    dst = Path(dst_dir).expanduser() / src.name
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
