from __future__ import annotations

import struct
from pathlib import Path
from typing import Any

from ._utils import _parse_shape


def ReadBinaryFileHeader(filePath: str):
    """Read the custom MERFISH binary-file header written by WriteBinaryFile."""
    path = Path(filePath)
    if not path.exists():
        raise ValueError("A valid path must be provided.")
    with path.open("rb") as fh:
        fixed = fh.read(10)
        if len(fixed) != 10:
            raise ValueError("File is too small to contain a valid binary header.")
        version, is_corrupt, num_entries, header_length = struct.unpack("<BBII", fixed)
        if version != 1:
            raise ValueError("The version is not supported.")
        header_bytes = fh.read(header_length)
        if len(header_bytes) != header_length:
            raise ValueError("File ended before the full variable header was read.")
    layout_string = header_bytes.decode("ascii")
    parts = [part.strip() for part in layout_string.split(",") if part.strip() != ""]
    if len(parts) % 3 != 0:
        raise ValueError("Invalid binary file layout header.")
    file_layout: list[dict[str, Any]] = []
    for i in range(0, len(parts), 3):
        field_name, shape_text, dtype_name = parts[i], parts[i + 1], parts[i + 2]
        file_layout.append({"dtype": dtype_name, "shape": _parse_shape(shape_text), "field": field_name})
    flat_header = {
        "version": version,
        "isCorrupt": bool(is_corrupt),
        "numEntries": int(num_entries),
        "headerLength": int(header_length),
    }
    return flat_header, file_layout
