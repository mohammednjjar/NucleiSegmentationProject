from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .LoadByteStream import LoadByteStream
from .LoadSplitByteStreamHeader import LoadSplitByteStreamHeader
from ._utils import make_parameters


def LoadSplitByteStream(filePath: str, **kwargs: Any):
    """Load split bytestream blocks saved by SaveSplitByteStream."""
    parameters = make_parameters({"verbose": True}, kwargs)
    if parameters.verbose:
        print(f"Loading split bytestream info from {filePath}")
    sbs_info, _ = LoadSplitByteStreamHeader(filePath)
    if sbs_info.get("version") not in {"1.0", "1.1"}:
        raise ValueError("The specified split bytestream version is not valid")
    base_path = Path(filePath).parent
    chunks = []
    for name in sbs_info["fileNames"]:
        block, _ = LoadByteStream(str(base_path / name), verbose=parameters.verbose)
        chunks.append(np.asarray(block))
    data = np.concatenate(chunks, axis=0) if chunks else np.array([])
    data = data.reshape(tuple(sbs_info["size"]), order="F")
    return data, sbs_info, parameters
