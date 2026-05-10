from __future__ import annotations

from pathlib import Path
from typing import Any
import warnings

from . import _utils
from ._utils import make_parameters


def SetFigureSavePath(newPath: str, **kwargs: Any):
    """Set or create the global figure-save path used by SaveFigure."""
    if newPath is None:
        raise ValueError("A data path is required.")
    parameters = make_parameters({"makeDir": False, "incrementDir": False, "verbose": True, "incrementDelimiter": "_"}, kwargs)
    path = Path(newPath)
    if parameters.makeDir:
        if path.exists() and parameters.incrementDir:
            count = 1
            base = path
            while path.exists():
                path = Path(f"{base}{parameters.incrementDelimiter}{count}")
                count += 1
        path.mkdir(parents=True, exist_ok=True)
        _utils.FIGURE_SAVE_PATH = str(path) + "/"
        save_path = _utils.FIGURE_SAVE_PATH
    else:
        if not path.exists():
            warnings.warn("Provided path does not exist", RuntimeWarning)
            save_path = ""
        else:
            _utils.FIGURE_SAVE_PATH = str(path) + "/"
            save_path = _utils.FIGURE_SAVE_PATH
    if parameters.verbose:
        print(f"Current Save Figure Path: {_utils.FIGURE_SAVE_PATH}")
    return save_path, parameters
