from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .ReadBinaryFileHeader import ReadBinaryFileHeader
from ._utils import _MATLAB_TO_NUMPY_DTYPE, _TYPE_SIZES, make_parameters, page_break


def _block_layout(file_layout):
    block_sizes = []
    offsets = [0]
    for item in file_layout:
        dtype_name = item["dtype"]
        if dtype_name not in _MATLAB_TO_NUMPY_DTYPE:
            raise TypeError(f"Found an invalid data type: {dtype_name}")
        count = int(np.prod(item["shape"]))
        size = _TYPE_SIZES[dtype_name] * count
        block_sizes.append(size)
        offsets.append(offsets[-1] + size)
    return block_sizes, offsets[:-1], offsets[-1]


def ReadBinaryFile(filePath: str, **kwargs: Any):
    """Read a MERFISH custom binary file into a list of dictionaries."""
    defaults = {"verbose": False, "format": "array", "first": None, "last": None, "fieldsToLoad": []}
    parameters = make_parameters(defaults, kwargs)

    path = Path(filePath)
    if not path.exists():
        raise ValueError("A valid path must be provided.")

    flat_header, file_layout = ReadBinaryFileHeader(str(path))
    if flat_header["isCorrupt"]:
        raise ValueError("File appears to be corrupt")

    first = 1 if parameters.first is None else int(parameters.first)
    last = flat_header["numEntries"] if parameters.last is None else int(parameters.last)
    if first < 1:
        raise ValueError("The first load index cannot be smaller than 1")
    if last > flat_header["numEntries"]:
        raise ValueError("The last load index cannot be larger than the number of entries in the file")
    if last < first:
        return [], flat_header, parameters

    fields_to_load = list(parameters.fieldsToLoad or [])
    if fields_to_load:
        available = {item["field"] for item in file_layout}
        missing = set(fields_to_load) - available
        if missing:
            raise ValueError(f"Some requested fields are absent from the binary file: {sorted(missing)}")

    if parameters.format not in {"array", "structureArray"}:
        raise ValueError("Only array and structureArray formats are supported in the Python translation.")

    block_sizes, block_offsets, total_block_size = _block_layout(file_layout)
    if parameters.verbose:
        page_break()
        print(f"Loading: {filePath}")
        print(f"    Writer version: {flat_header['version']}")
        print(f"    Number of entries: {flat_header['numEntries']}")
        print(f"    Loading: {first} to {last}")
        print("    Found fields:")
        for item in file_layout:
            print(f"        {item['field']}: {item['dtype']} of size {item['shape']}")

    data_start = 10 + flat_header["headerLength"]
    records = []
    with path.open("rb") as fh:
        for entry_idx in range(first - 1, last):
            rec: dict[str, Any] = {}
            base = data_start + total_block_size * entry_idx
            for field_idx, item in enumerate(file_layout):
                field = item["field"]
                if fields_to_load and field not in fields_to_load:
                    continue
                dtype = _MATLAB_TO_NUMPY_DTYPE[item["dtype"]]
                count = int(np.prod(item["shape"]))
                fh.seek(base + block_offsets[field_idx])
                raw = fh.read(block_sizes[field_idx])
                if len(raw) != block_sizes[field_idx]:
                    raise ValueError("Binary file ended before all requested data were read.")
                arr = np.frombuffer(raw, dtype=dtype, count=count).copy().reshape(item["shape"], order="F")
                rec[field] = arr.item() if arr.size == 1 else arr
            records.append(rec)

    if parameters.format == "structureArray":
        combined = {field: [] for field in (fields_to_load or [x["field"] for x in file_layout])}
        for rec in records:
            for field, value in rec.items():
                combined[field].append(value)
        return combined, flat_header, parameters
    return records, flat_header, parameters
