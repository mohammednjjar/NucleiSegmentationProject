from __future__ import annotations

import numpy as np

from .._utils import _copy_value, _default_molecule_list


def CreateWordsStructure(numElements: int, numHybs: int):
    """Create a list of empty word dictionaries."""
    if numElements < 0 or numHybs < 0:
        raise ValueError("Invalid values for numElements and numHybs")
    default = {
        "uID": "",
        "codeword": np.zeros(numHybs, dtype=bool),
        "measuredCodeword": np.zeros(numHybs, dtype=bool),
        "intCodeword": np.nan,
        "geneName": "",
        "isExactMatch": False,
        "isCorrectedMatch": False,
        "numOnBits": 0,
        "onBits": np.full(numHybs, np.nan),
        "measuredOnBits": np.full(numHybs, np.nan),
        "bitOrder": np.arange(1, numHybs + 1),
        "numHyb": np.nan,
        "imageNames": [None] * numHybs,
        "imagePaths": [None] * numHybs,
        "imageUIDs": [None] * numHybs,
        "imageX": np.nan,
        "imageY": np.nan,
        "wordCentroidX": np.nan,
        "wordCentroidY": np.nan,
        "cellID": np.nan,
        "wordNumInCell": np.nan,
        "fiducialUIDs": [None] * numHybs,
        "hasFiducialError": np.zeros(numHybs, dtype=bool),
        "paddedCellID": np.full(numHybs, np.nan),
        "focusLockQuality": np.full(numHybs, np.nan),
        "mListInds": np.full(numHybs, np.nan),
    }
    default.update({key: value.copy() for key, value in _default_molecule_list(numHybs).items()})
    return [_copy_value(default) for _ in range(int(numElements))]
