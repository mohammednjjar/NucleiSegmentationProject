from __future__ import annotations

import shutil
import struct
from pathlib import Path
from typing import Any

from .ReadBinaryFileHeader import ReadBinaryFileHeader
from ._utils import _TYPE_SIZES, make_parameters


def CatBinaryFiles(filePath1: str, filePath2: str, **kwargs: Any):
    """Append the records in binary file 2 to binary file 1 after layout validation."""
    parameters = make_parameters({"verbose": False}, kwargs)
    src = Path(filePath2)
    dst = Path(filePath1)
    if not src.exists():
        raise ValueError("Two paths must be provided and the path to the file to add must be valid")
    if not dst.exists():
        shutil.copyfile(src, dst)
        if parameters.verbose:
            print("File 1 does not exist. So File 2 was copied to File 1")
        return parameters

    header1, layout1 = ReadBinaryFileHeader(str(dst))
    header2, layout2 = ReadBinaryFileHeader(str(src))
    if header1["isCorrupt"] or header2["isCorrupt"]:
        raise ValueError("One of the binary files is corrupt")
    if header1["version"] != header2["version"]:
        raise ValueError("The two binary files do not have the same version")
    if header1["version"] != 1:
        raise ValueError("This version is not yet supported.")
    if header1["headerLength"] != header2["headerLength"] or layout1 != layout2:
        raise ValueError("The two binary files do not have the same layout")

    with dst.open("r+b") as fh:
        fh.seek(1)
        fh.write(struct.pack("<B", 1))
    data_offset = 10 + header2["headerLength"]
    with src.open("rb") as fh:
        fh.seek(data_offset)
        data = fh.read()
    with dst.open("ab") as fh:
        fh.write(data)
    with dst.open("r+b") as fh:
        fh.seek(1)
        fh.write(struct.pack("<BI", 0, header1["numEntries"] + header2["numEntries"]))
    return parameters
