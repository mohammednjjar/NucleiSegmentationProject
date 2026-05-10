from __future__ import annotations

from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np

from .SaveAsByteStream import SaveAsByteStream
from ._utils import make_parameters


def SaveSplitByteStream(filePath: str, variable: Any, blockSize: int, **kwargs: Any):
    """Save a flattened array-like object as several .matb byte-stream blocks."""
    parameters = make_parameters({"verbose": True, "overwrite": True, "splitPostFix": "_"}, kwargs)
    if not isinstance(filePath, str):
        raise ValueError("The first argument must be a file path")
    if blockSize <= 0:
        raise ValueError("blockSize must be positive")

    path = Path(filePath)
    ext = path.suffix or ".matb"
    base_dir = path.parent if str(path.parent) else Path(".")
    file_name = path.stem if path.suffix else path.name
    base_dir.mkdir(parents=True, exist_ok=True)

    arr = np.asarray(variable)
    original_shape = arr.shape
    flat = arr.reshape(-1, order="F")
    indices = list(range(0, flat.size, int(blockSize))) + [flat.size]
    if parameters.verbose:
        print(f"Splitting file into {len(indices) - 1} chunks")

    sbs_info = {
        "size": original_shape,
        "numBlocks": len(indices) - 1,
        "splitPostFix": parameters.splitPostFix,
        "version": "1.1",
        "isByteStream": "Y",
        "fileNames": [],
        "uuID": str(uuid4()),
    }

    save_base = f"{file_name}{parameters.splitPostFix}"
    for i in range(len(indices) - 1):
        local_name = f"{save_base}{i + 1}{ext}"
        sbs_info["fileNames"].append(local_name)
        SaveAsByteStream(str(base_dir / local_name), flat[indices[i]:indices[i + 1]], verbose=parameters.verbose, overwrite=parameters.overwrite)

    SaveAsByteStream(str(base_dir / f"{file_name}{ext}"), sbs_info, verbose=parameters.verbose, overwrite=parameters.overwrite)
    return save_base, parameters
