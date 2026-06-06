"""Shared Python utilities for the translated MERFISH_analysis fileIO folder."""
from __future__ import annotations

import csv
import os
import pickle
import re
import shutil
import struct
import warnings
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Iterable, Sequence
from uuid import uuid4

import numpy as np


FIGURE_SAVE_PATH: str = ""

_DTYPE_ALIASES = {
    "single": "float32",
    "double": "float64",
    "float": "float64",
    "logical": "uint8",
    "bool": "uint8",
    "boolean": "uint8",
    "char": "uint8",
    "int8": "int8",
    "int16": "int16",
    "int32": "int32",
    "int64": "int64",
    "uint8": "uint8",
    "uint16": "uint16",
    "uint32": "uint32",
    "uint64": "uint64",
}

_SUPPORTED_BINARY_DTYPES = {
    "uint8", "uint16", "uint32", "uint64",
    "int8", "int16", "int32", "int64",
    "single", "double",
}

_NUMPY_TO_MATLAB_DTYPE = {
    "uint8": "uint8",
    "uint16": "uint16",
    "uint32": "uint32",
    "uint64": "uint64",
    "int8": "int8",
    "int16": "int16",
    "int32": "int32",
    "int64": "int64",
    "float32": "single",
    "float64": "double",
    "bool": "uint8",
}

_MATLAB_TO_NUMPY_DTYPE = {
    "uint8": np.uint8,
    "uint16": np.uint16,
    "uint32": np.uint32,
    "uint64": np.uint64,
    "int8": np.int8,
    "int16": np.int16,
    "int32": np.int32,
    "int64": np.int64,
    "single": np.float32,
    "double": np.float64,
}

_TYPE_SIZES = {
    "uint8": 1, "uint16": 2, "uint32": 4, "uint64": 8,
    "int8": 1, "int16": 2, "int32": 4, "int64": 8,
    "single": 4, "double": 8,
}


def make_parameters(defaults: dict[str, Any], overrides: dict[str, Any] | None = None) -> SimpleNamespace:
    params = dict(defaults)
    if overrides:
        if "parameters" in overrides:
            supplied = overrides.pop("parameters")
            if isinstance(supplied, SimpleNamespace):
                params.update(vars(supplied))
            elif isinstance(supplied, dict):
                params.update(supplied)
            else:
                params.update(vars(supplied))
        params.update(overrides)
    return SimpleNamespace(**params)


def page_break() -> None:
    print("=" * 80)


def ensure_dir(path: str | os.PathLike[str]) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def increment_save_name(file_path: str | os.PathLike[str]) -> str:
    path = Path(file_path)
    base = path.with_suffix("")
    suffix = path.suffix
    count = 1
    while True:
        candidate = Path(f"{base}_{count}{suffix}")
        if not candidate.exists():
            return str(candidate)
        count += 1


def _matlab_shape(value: Any) -> tuple[int, ...]:
    arr = np.asarray(value)
    if arr.ndim == 0:
        return (1,)
    return tuple(int(x) for x in arr.shape)


def _shape_to_string(shape: Sequence[int]) -> str:
    return " ".join(str(int(x)) for x in shape)


def _parse_shape(text: str) -> tuple[int, ...]:
    nums = re.findall(r"[-+]?\d+", text)
    if not nums:
        return (1,)
    return tuple(int(x) for x in nums)


def _field_dtype_name(value: Any) -> str:
    arr = np.asarray(value)
    dtype_name = arr.dtype.name
    if dtype_name not in _NUMPY_TO_MATLAB_DTYPE:
        raise TypeError(f"Unsupported field dtype {arr.dtype!s}; use numeric NumPy-compatible fields only.")
    return _NUMPY_TO_MATLAB_DTYPE[dtype_name]


def _as_field_array(value: Any, dtype_name: str) -> np.ndarray:
    np_dtype = _MATLAB_TO_NUMPY_DTYPE[dtype_name]
    arr = np.asarray(value, dtype=np_dtype)
    if arr.ndim == 0:
        arr = arr.reshape((1,))
    return arr


def _copy_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.copy()
    if isinstance(value, list):
        return [_copy_value(v) for v in value]
    if isinstance(value, dict):
        return {k: _copy_value(v) for k, v in value.items()}
    return value


def _get_header_value(readout: Any) -> str:
    if isinstance(readout, dict):
        if "Header" not in readout:
            raise ValueError("Each readout dict must contain a 'Header' key.")
        return str(readout["Header"])
    if hasattr(readout, "Header"):
        return str(getattr(readout, "Header"))
    raise ValueError("Each readout must be a dict/object with a Header field.")


def _split_preserve_delimiters(text: str, delimiters: Sequence[str]) -> tuple[list[str], list[str]]:
    if not delimiters:
        return [text], []
    pattern = "(" + "|".join(re.escape(d) for d in delimiters) + ")"
    parts = re.split(pattern, text)
    tokens = parts[0::2]
    found = parts[1::2]
    if tokens and tokens[-1] == "":
        tokens = tokens[:-1]
    return tokens, found


def _apply_conv(value: str, conv: Callable[[str], Any] | type | str | None) -> Any:
    if conv is None:
        return value
    if conv in (str, "char", "string"):
        return str(value)
    if conv in (int, "int", "integer"):
        return int(value)
    if conv in (float, "float", "double", "single"):
        return float(value)
    if conv in (bool, "bool", "boolean"):
        lowered = str(value).strip().lower()
        return lowered in {"1", "true", "t", "yes", "y"}
    return conv(value)


def _normalize_records(struct_array: Sequence[dict[str, Any]] | dict[str, Any]) -> list[dict[str, Any]]:
    if isinstance(struct_array, dict):
        return [struct_array]
    records = [dict(item) for item in struct_array]
    if not records:
        raise ValueError("structArray must contain at least one record.")
    return records


def _default_molecule_list(num_elements: int) -> dict[str, np.ndarray]:
    if num_elements < 0:
        raise ValueError("numElements must be non-negative.")
    float_fields = [
        "x", "y", "xc", "yc", "h", "a", "w", "phi", "ax", "bg", "i",
        "c", "density", "photons", "z", "zIndex",
    ]
    int_fields = ["frame", "length", "link", "category"]
    out: dict[str, np.ndarray] = {k: np.full(num_elements, np.nan, dtype=float) for k in float_fields}
    out.update({k: np.zeros(num_elements, dtype=np.int32) for k in int_fields})
    return out


def _write_pickle(path: str | os.PathLike[str], obj: Any) -> None:
    with open(path, "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)


def _read_pickle(path: str | os.PathLike[str]) -> Any:
    with open(path, "rb") as fh:
        return pickle.load(fh)


def _parse_codebook_header_rows(rows: list[list[str]]) -> tuple[dict[str, Any], int]:
    header: dict[str, Any] = {}
    data_start = -1
    for idx, row in enumerate(rows):
        parts = [part.strip() for part in row]
        if set(["name", "id", "barcode"]).issubset(set(parts)) and len(parts) >= 3:
            data_start = idx + 1
            break
        if len(parts) == 2:
            header[parts[0]] = parts[1]
        elif parts:
            header[parts[0]] = parts[1:]
    if data_start < 0:
        raise ValueError("Codebook header row 'name, id, barcode' was not found.")
    if "version" not in header or "bit_names" not in header:
        raise ValueError("The codebook is corrupt: both version and bit_names are required.")
    return header, data_start


def _read_csv_rows(path: str | os.PathLike[str]) -> list[list[str]]:
    with open(path, newline="") as fh:
        return [row for row in csv.reader(fh)]


def _write_csv_rows(path: str | os.PathLike[str], rows: Iterable[Sequence[Any]]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)
