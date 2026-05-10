"""Python translation of MERFISH_analysis/probe_construction/OTTable.m."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterable, List, Sequence
import numpy as np

try:
    from .OTMap import OTMap
    from .OTMap2 import OTMap2
    from .utils import from_possible_fasta_struct, rolling_hash, save_pickle, load_pickle, ensure_dir
except ImportError:
    from OTMap import OTMap
    from OTMap2 import OTMap2
    from utils import from_possible_fasta_struct, rolling_hash, save_pickle, load_pickle, ensure_dir


class OTTable:
    """Off-target table for seed-sequence penalties.

    MATLAB purpose: count every valid n-mer in target sequences, optionally using
    per-sequence weights, and use those counts as off-target penalties.
    """

    def __init__(
        self,
        targetSequences: Any = None,
        seedLength: int | None = None,
        weights: Sequence[float] | None = None,
        verbose: bool = False,
        mapType: type | None = None,
        convertMapType: bool = True,
        parallel: Any = None,
        transferAbund: bool = False,
        name: str = "",
    ):
        self.name = name
        self.verbose = bool(verbose)
        self.seedLength = int(seedLength) if seedLength is not None else 0
        self.numPar = 0
        self.parallel = parallel
        self.mapType = mapType if mapType is not None else OTMap
        self.numEntries = 0
        self.uniformWeight = False
        self.data = self.mapType()
        self.hashBase = np.asarray([], dtype=np.int64)

        if targetSequences is None:
            return
        if seedLength is None or int(seedLength) <= 0:
            raise ValueError("seedLength must be a positive integer")

        self.seedLength = int(seedLength)
        self.hashBase = 4 ** np.arange(self.seedLength - 1, -1, -1, dtype=np.int64)
        self.data = self.mapType()

        if hasattr(targetSequences, "intSequences"):
            seqs = [np.asarray(s, dtype=np.int16) for s in targetSequences.intSequences]
            if transferAbund and getattr(targetSequences, "abundLoaded", False):
                weights = getattr(targetSequences, "abundance")
        else:
            seqs = from_possible_fasta_struct(targetSequences)

        if len(seqs) == 0:
            return
        if weights is None or len(weights) == 0:
            weights_arr = np.ones(len(seqs), dtype=float)
            self.uniformWeight = True
        else:
            weights_arr = np.asarray(weights, dtype=float).ravel()
            if weights_arr.size != len(seqs):
                raise ValueError("weights must be equal in length to targetSequences")
            self.uniformWeight = False

        key_values: List[np.ndarray] = []
        for seq, weight in zip(seqs, weights_arr):
            hashes, valid = rolling_hash(seq, self.seedLength)
            valid_hashes = hashes[valid]
            if valid_hashes.size:
                key_values.append(np.vstack([valid_hashes.astype(float), np.full(valid_hashes.size, float(weight))]))
        if key_values:
            all_data = np.hstack(key_values)
            self.data.AddToMap(all_data)
        if convertMapType and isinstance(self.data, OTMap):
            self.data = OTMap2(self.data.GetTable())
            self.mapType = OTMap2
        self.numEntries = len(self.data)

    def IsValidSequence(self, localSeq: str | Sequence[int] | np.ndarray) -> np.ndarray | int:
        hashes, valid = rolling_hash(localSeq, self.seedLength)
        if hashes.size == 0:
            return 0
        return valid

    def CalculatePenalty(self, seq: str | Sequence[int] | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        hashes, valid = rolling_hash(seq, self.seedLength)
        if hashes.size == 0:
            return np.asarray([], dtype=float), np.asarray([], dtype=float)
        penalty = np.full(hashes.size, np.nan, dtype=float)
        if np.any(valid):
            penalty[valid] = self.data.GetValues(hashes[valid].astype(float))
        return penalty, hashes

    def calculate_penalty(self, seq: str | Sequence[int] | np.ndarray) -> np.ndarray:
        return self.CalculatePenalty(seq)[0]

    def keys(self) -> np.ndarray:
        return self.data.keys()

    def values(self) -> np.ndarray:
        return self.data.values()

    def SetParallel(self, p: Any) -> None:
        self.parallel = p
        self.numPar = int(getattr(p, "NumWorkers", 0) or 0)
        return None

    def __add__(self, other: "OTTable") -> "OTTable":
        if not isinstance(other, OTTable):
            raise TypeError("Cannot add OTTable with an object of a different class")
        if self.seedLength != other.seedLength:
            raise ValueError("Cannot add two OTTables with different seed lengths")
        out = OTTable([], self.seedLength, mapType=self.mapType, verbose=self.verbose)
        data = np.hstack([self.data.GetTable(), other.data.GetTable()])
        out.data = self.mapType()
        out.data.AddToMap(data)
        out.numEntries = len(out.data)
        out.uniformWeight = bool(self.uniformWeight and other.uniformWeight)
        return out

    @staticmethod
    def sum(otTables: Sequence["OTTable"]) -> "OTTable":
        tables = list(otTables)
        if not tables:
            return OTTable()
        seed_lengths = {t.seedLength for t in tables}
        if len(seed_lengths) != 1:
            raise ValueError("otTable arrays must have the same seed length")
        out = OTTable([], tables[0].seedLength, mapType=tables[0].mapType, verbose=tables[0].verbose)
        all_data = [t.data.GetTable() for t in tables if t.numEntries]
        if all_data:
            out.data.AddToMap(np.hstack(all_data))
        out.numEntries = len(out.data)
        out.uniformWeight = all(t.uniformWeight for t in tables)
        return out

    def Save(self, dirPath: str | Path) -> None:
        path = ensure_dir(dirPath)
        payload = {
            "name": self.name,
            "verbose": self.verbose,
            "seedLength": self.seedLength,
            "numEntries": self.numEntries,
            "uniformWeight": self.uniformWeight,
            "mapType": self.mapType.__name__,
            "hashBase": self.hashBase,
            "data": self.data.GetTable(),
        }
        save_pickle(path / "OTTable.pkl", payload)
        np.asarray(payload["data"], dtype=float).tofile(path / "data.bin")
        save_pickle(path / "numEntries.matb", self.numEntries)
        save_pickle(path / "seedLength.matb", self.seedLength)
        save_pickle(path / "uniformWeight.matb", self.uniformWeight)
        save_pickle(path / "name.matb", self.name)
        return None

    @staticmethod
    def Load(dirPath: str | Path, verbose: bool = True, mapType: type | None = None) -> "OTTable":
        path = Path(dirPath)
        if not path.is_dir():
            raise ValueError("The provided path is not valid")
        payload_path = path / "OTTable.pkl"
        if payload_path.exists():
            payload = load_pickle(payload_path)
            cls = mapType if mapType is not None else (OTMap2 if payload.get("mapType") == "OTMap2" else OTMap)
            obj = OTTable([], int(payload["seedLength"]), mapType=cls, verbose=verbose)
            obj.name = payload.get("name", "")
            obj.uniformWeight = bool(payload.get("uniformWeight", False))
            obj.hashBase = np.asarray(payload.get("hashBase", 4 ** np.arange(obj.seedLength - 1, -1, -1)), dtype=np.int64)
            obj.data = cls(np.asarray(payload["data"], dtype=float))
            obj.numEntries = len(obj.data)
            return obj
        seed_length = int(load_pickle(path / "seedLength.matb"))
        num_entries = int(load_pickle(path / "numEntries.matb"))
        cls = mapType if mapType is not None else OTMap2
        obj = OTTable([], seed_length, mapType=cls, verbose=verbose)
        if (path / "data.bin").exists():
            arr = np.fromfile(path / "data.bin", dtype=float)
            if arr.size:
                obj.data = cls(arr.reshape(2, num_entries))
        obj.numEntries = len(obj.data)
        obj.uniformWeight = bool(load_pickle(path / "uniformWeight.matb")) if (path / "uniformWeight.matb").exists() else False
        obj.name = load_pickle(path / "name.matb") if (path / "name.matb").exists() else ""
        return obj
