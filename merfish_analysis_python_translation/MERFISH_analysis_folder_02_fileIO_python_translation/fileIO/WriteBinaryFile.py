from __future__ import annotations

import os
import struct
from pathlib import Path
from typing import Any

import numpy as np

from .ReadBinaryFileHeader import ReadBinaryFileHeader
from ._utils import (
    _SUPPORTED_BINARY_DTYPES,
    _as_field_array,
    _field_dtype_name,
    _matlab_shape,
    _shape_to_string,
    make_parameters,
)


def _layout_from_records(records: list[dict[str, Any]]):
    fields = list(records[0].keys())
    layout = []
    for field in fields:
        dtype_name = _field_dtype_name(records[0][field])
        if dtype_name not in _SUPPORTED_BINARY_DTYPES:
            raise TypeError(f"Unsupported data type for field {field}: {dtype_name}")
        shape = _matlab_shape(records[0][field])
        for rec in records[1:]:
            if field not in rec:
                raise ValueError(f"Field {field} is missing from one record.")
            if _field_dtype_name(rec[field]) != dtype_name or _matlab_shape(rec[field]) != shape:
                raise ValueError(f"Field {field} must have the same dtype and shape in every record.")
        layout.append({"field": field, "shape": shape, "dtype": dtype_name})
    return layout


def WriteBinaryFile(filePath: str, structArray, **kwargs: Any):
    """Write a list of dictionaries to the MERFISH custom binary format."""
    defaults = {"verbose": False, "overwrite": True, "append": False}
    parameters = make_parameters(defaults, kwargs)

    if structArray is None:
        raise ValueError("A path and a structure array must be provided.")
    records = [dict(structArray)] if isinstance(structArray, dict) else [dict(x) for x in structArray]
    if not records:
        raise ValueError("A structure array must contain at least one entry.")

    version = 1
    layout = _layout_from_records(records)
    header_parts = []
    for item in layout:
        header_parts.extend([item["field"], _shape_to_string(item["shape"]), item["dtype"]])
    header_string = ",".join(header_parts)
    header_bytes = header_string.encode("ascii")
    path = Path(filePath)

    is_append = False
    old_num_entries = 0
    if path.exists():
        if parameters.overwrite and not parameters.append:
            path.unlink()
            if parameters.verbose:
                print(f"Deleting existing file: {filePath}")
        elif not parameters.append:
            raise FileExistsError("Found existing file.")
        else:
            flat_header, old_layout = ReadBinaryFileHeader(str(path))
            if flat_header["isCorrupt"]:
                raise ValueError("File appears to be corrupt")
            if len(old_layout) != len(layout):
                raise ValueError("Cannot append binary data with a different organization")
            for a, b in zip(old_layout, layout):
                if a["field"] != b["field"] or tuple(a["shape"]) != tuple(b["shape"]) or a["dtype"] != b["dtype"]:
                    raise ValueError("Cannot append binary data with a different organization")
            is_append = True
            old_num_entries = flat_header["numEntries"]

    if not is_append:
        with path.open("wb") as fh:
            fh.write(struct.pack("<BBII", version, 1, 0, len(header_bytes)))
            fh.write(header_bytes)
    else:
        with path.open("r+b") as fh:
            fh.seek(1)
            fh.write(struct.pack("<B", 1))

    with path.open("ab") as fh:
        for rec in records:
            for item in layout:
                arr = _as_field_array(rec[item["field"]], item["dtype"])
                fh.write(np.asfortranarray(arr).tobytes(order="F"))

    total_entries = old_num_entries + len(records)
    with path.open("r+b") as fh:
        fh.seek(1)
        fh.write(struct.pack("<BI", 0, total_entries))
    return parameters
