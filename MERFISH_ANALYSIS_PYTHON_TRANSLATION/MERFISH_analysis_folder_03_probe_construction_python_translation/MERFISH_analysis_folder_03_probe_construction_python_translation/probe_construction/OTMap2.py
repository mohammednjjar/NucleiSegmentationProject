"""Python translation of MERFISH_analysis/probe_construction/OTMap2.m."""
from __future__ import annotations

from typing import Sequence
import numpy as np


class OTMap2:
    """Dictionary-backed key/value map used by OTTable.

    MATLAB purpose: containers.Map implementation for faster lookup after an
    OTTable is built.  Duplicate keys are accumulated.
    """

    def __init__(self, initialData: Sequence[Sequence[float]] | np.ndarray | None = None):
        self.data: dict[float, float] = {}
        if initialData is not None:
            arr = np.asarray(initialData, dtype=float)
            if arr.size == 0:
                return
            if arr.ndim != 2 or arr.shape[0] != 2:
                raise ValueError("initialData must be a numeric 2xN array")
            self.AddToMap(arr)

    def AddToMap(self, newData: Sequence[Sequence[float]] | np.ndarray) -> None:
        arr = np.asarray(newData, dtype=float)
        if arr.size == 0:
            return None
        if arr.ndim != 2 or arr.shape[0] != 2:
            raise ValueError("newData must be a numeric 2xN array")
        for key, value in zip(arr[0, :], arr[1, :]):
            k = float(key)
            self.data[k] = self.data.get(k, 0.0) + float(value)
        return None

    def GetValues(self, keys: Sequence[float] | np.ndarray) -> np.ndarray:
        query = np.asarray(keys, dtype=float).ravel()
        return np.asarray([self.data.get(float(k), 0.0) for k in query], dtype=float)

    def GetTable(self) -> np.ndarray:
        if not self.data:
            return np.zeros((2, 0), dtype=float)
        keys = np.asarray(sorted(self.data.keys()), dtype=float)
        values = np.asarray([self.data[float(k)] for k in keys], dtype=float)
        return np.vstack([keys, values])

    def keys(self) -> np.ndarray:
        return self.GetTable()[0, :]

    def values(self) -> np.ndarray:
        return self.GetTable()[1, :]

    def length(self) -> int:
        return len(self.data)

    def __len__(self) -> int:
        return self.length()
