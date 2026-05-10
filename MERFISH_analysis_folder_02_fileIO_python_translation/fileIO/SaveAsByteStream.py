from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

from ._utils import _write_pickle, increment_save_name, make_parameters


def SaveAsByteStream(filePath: str, variable: Any, **kwargs: Any):
    """Serialize a Python object to a .matb-style file using pickle."""
    parameters = make_parameters({"verbose": True, "overwrite": True}, kwargs)
    if not isinstance(filePath, str):
        raise ValueError("The first argument must be a file path")
    path = Path(filePath)
    if path.suffix != ".matb":
        warnings.warn("It is strongly recommended to save bytestreams as .matb files", RuntimeWarning)
    if not parameters.overwrite:
        while path.exists():
            path = Path(increment_save_name(path))
            if parameters.verbose:
                warnings.warn("Detected existing file. Increased counter on save file.", RuntimeWarning)
    path.parent.mkdir(parents=True, exist_ok=True)
    if parameters.verbose:
        print(f"Saving {path}")
    _write_pickle(path, variable)
    if parameters.verbose:
        print(".... finished")
    return str(path), parameters
