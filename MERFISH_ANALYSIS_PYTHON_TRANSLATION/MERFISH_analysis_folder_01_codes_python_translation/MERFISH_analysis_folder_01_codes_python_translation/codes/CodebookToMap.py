"""Python translation of `codes/CodebookToMap.m`."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np

try:
    from .hamming_utils import bits_to_int
except ImportError:  # direct script execution
    from hamming_utils import bits_to_int


@dataclass(frozen=True)
class FastaRecord:
    """Minimal FASTA record matching MATLAB fastaread fields."""
    Header: str
    Sequence: str


def _read_fasta(path: str | Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    seq_parts: list[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            records.append(FastaRecord(header, "".join(seq_parts)))
    if not records:
        raise ValueError("The provided path is not a valid FASTA file.")
    return records


def _normalise_codebook(codebook) -> list[FastaRecord]:
    if isinstance(codebook, (str, Path)):
        path = Path(codebook)
        if not path.is_file():
            raise ValueError("The provided path is not to a valid file.")
        return _read_fasta(path)

    records: list[FastaRecord] = []
    for item in codebook:
        if isinstance(item, FastaRecord):
            records.append(item)
        elif isinstance(item, Mapping):
            if "Header" not in item or "Sequence" not in item:
                raise ValueError("The provided structure does not have Header and Sequence fields.")
            records.append(FastaRecord(str(item["Header"]), str(item["Sequence"])))
        else:
            if not hasattr(item, "Header") or not hasattr(item, "Sequence"):
                raise ValueError("The provided structure does not have Header and Sequence fields.")
            records.append(FastaRecord(str(item.Header), str(item.Sequence)))
    return records


def _remove_ws(text: str) -> str:
    return "".join(str(text).split())


def _key_converter(keyType: str):
    if keyType == "int":
        def conv(x) -> int:
            s = _remove_ws("".join(map(str, x)) if not isinstance(x, str) else x)
            if not set(s).issubset({"0", "1"}):
                raise ValueError(f"Invalid binary codeword: {x!r}")
            return bits_to_int([int(ch) for ch in s], left_msb=True)
        return conv
    if keyType == "binStr":
        return _remove_ws
    raise ValueError("keyType must be 'int' or 'binStr'.")


def _gene_name_from_sequence(sequence: str) -> str:
    parts = sequence.split()
    if not parts:
        raise ValueError("Codebook sequence is empty; cannot extract gene name.")
    return parts[0]


def CodebookToMap(codebook,
                  errCorrFunc: Callable[[np.ndarray], Iterable[Sequence[bool]]] | None = None,
                  keyType: str = "int",
                  mapContents: str = "all"):
    """Build a Python dict from a MERFISH codebook.

    Parameters
    ----------
    codebook:
        FASTA path, or iterable of objects/dicts with `Header` and `Sequence`.
        Headers contain binary codewords; the first whitespace-separated token
        in `Sequence` is used as the gene/object name.
    errCorrFunc:
        Optional function that receives a boolean codeword and returns
        correctable codewords. Equivalent to MATLAB `errCorrFunc`.
    keyType:
        `'int'` returns integer keys; `'binStr'` returns binary-string keys.
    mapContents:
        `'exact'`, `'correctable'`, or `'all'`.

    Returns
    -------
    tuple
        `(mapping, geneNames, codewords, parameters)`.
    """
    if codebook is None:
        raise ValueError("A codebook is required.")
    if mapContents not in {"exact", "all", "correctable"}:
        raise ValueError("mapContents must be 'exact', 'all', or 'correctable'.")

    records = _normalise_codebook(codebook)
    keyConv = _key_converter(keyType)

    exact_map: Dict[int | str, str] = {}
    geneNames: list[str] = []
    codewords: list[int | str] = []

    for rec in records:
        key = keyConv(rec.Header)
        value = _gene_name_from_sequence(rec.Sequence)
        exact_map[key] = value
        codewords.append(key)
        geneNames.append(value)

    correctable_map: Dict[int | str, str] = {}
    if errCorrFunc is not None:
        for rec in records:
            logical_word = np.asarray([ch == "1" for ch in _remove_ws(rec.Header)], dtype=bool)
            value = _gene_name_from_sequence(rec.Sequence)
            for new_word in errCorrFunc(logical_word):
                new_key = keyConv("".join("1" if b else "0" for b in np.asarray(new_word, dtype=bool)))
                correctable_map[new_key] = value

    if errCorrFunc is None or mapContents == "exact":
        mapping = exact_map
    elif mapContents == "correctable":
        mapping = correctable_map
    else:
        mapping = {**exact_map, **correctable_map}

    parameters = {
        "errCorrFunc": errCorrFunc,
        "keyType": keyType,
        "mapContents": mapContents,
    }
    return mapping, geneNames, codewords, parameters
