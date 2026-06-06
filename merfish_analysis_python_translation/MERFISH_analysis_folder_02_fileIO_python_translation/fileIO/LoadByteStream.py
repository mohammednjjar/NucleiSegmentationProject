from __future__ import annotations

from pathlib import Path
from typing import Any

from ._utils import _read_pickle, make_parameters


def LoadByteStream(filePath: str, **kwargs: Any):
    """Load a Python object saved by SaveAsByteStream."""
    parameters = make_parameters({"verbose": True}, kwargs)
    path = Path(filePath)
    if not path.exists():
        raise ValueError("Invalid file path")
    if parameters.verbose:
        print(f"Loading {filePath}")
    data = _read_pickle(path)
    if parameters.verbose:
        print(".... finished")
    return data, parameters
