from __future__ import annotations

from typing import Any
import numpy as np

from .._utils import _copy_value, _default_molecule_list


def CreateImageDataStructure(numElements: int):
    """Create a list of empty imageData dictionaries."""
    if numElements < 0:
        raise ValueError("Invalid values for numElements and numHybs")
    default = {
        "name": "",
        "filePath": "",
        "infFilePath": "",
        "uID": "",
        "movieType": "",
        "hybNum": -1,
        "cellNum": -1,
        "isFiducial": False,
        "binType": "",
        "delimiters": [],
        "imageH": 0,
        "imageW": 0,
        "Stage_X": 0,
        "Stage_Y": 0,
        "focusLockQuality": 0,
        "mList": _default_molecule_list(0),
        "tform": np.eye(3),
        "warpErrors": np.zeros(5),
        "hasFiducialError": False,
        "fiducialErrorMessage": None,
        "fidUID": "",
    }
    return [_copy_value(default) for _ in range(int(numElements))]
