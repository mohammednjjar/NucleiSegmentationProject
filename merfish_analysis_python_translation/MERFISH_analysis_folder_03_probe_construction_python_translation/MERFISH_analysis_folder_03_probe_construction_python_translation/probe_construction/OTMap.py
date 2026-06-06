"""Python translation of MERFISH_analysis/probe_construction/OTMap.m."""
from __future__ import annotations

from typing import Iterable, Sequence
import numpy as np


class OTMap:
    """Array-backed key/value map used by OTTable.

    MATLAB purpose: store key/value pairs as a 2xN numeric table and accumulate
    duplicate key values.  This Python class keeps the same public methods.
    """

    def __init__(self, initialData: Sequence[Sequence[float]] | np.ndarray | None = None):
        if initialData is None:
            self.data = np.zeros((2, 0), dtype=float)
        else:
            arr = np.asarray(initialData, dtype=float)
            if arr.size == 0:
                arr = np.zeros((2, 0), dtype=float)
            if arr.ndim != 2 or arr.shape[0] != 2:
                raise ValueError("initialData must be a numeric 2xN array")
            self.data = self._accumulate(arr)

    @staticmethod
    def _accumulate(arr: np.ndarray) -> np.ndarray:
        if arr.size == 0:
            return np.zeros((2, 0), dtype=float)
        keys, inverse = np.unique(arr[0, :], return_inverse=True)
        values = np.zeros(keys.size, dtype=float)
        np.add.at(values, inverse, arr[1, :].astype(float))
        return np.vstack([keys.astype(float), values])

    def AddToMap(self, newData: Sequence[Sequence[float]] | np.ndarray) -> None:
        arr = np.asarray(newData, dtype=float)
        if arr.size == 0:
            return None
        if arr.ndim != 2 or arr.shape[0] != 2:
            raise ValueError("newData must be a numeric 2xN array")
        self.data = self._accumulate(np.hstack([self.data, arr]))
        return None

    def GetValues(self, keys: Sequence[float] | np.ndarray) -> np.ndarray:
        query = np.asarray(keys, dtype=float).ravel()
        lookup = {float(k): float(v) for k, v in zip(self.data[0, :], self.data[1, :])}
        return np.asarray([lookup.get(float(k), 0.0) for k in query], dtype=float)

    def GetTable(self) -> np.ndarray:
        return self.data.copy()

    def keys(self) -> np.ndarray:
        return self.data[0, :].copy()

    def values(self) -> np.ndarray:
        return self.data[1, :].copy()

    def length(self) -> int:
        return int(self.data.shape[1])

    def __len__(self) -> int:
        return self.length()
