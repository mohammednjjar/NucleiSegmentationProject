"""Python translation of MERFISH_analysis/probe_construction/TargetRegions.m."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, List, Sequence
import numpy as np

try:
    from .utils import ensure_dir, fasta_write, int_to_sequence, load_pickle, save_pickle, sequence_to_int
except ImportError:
    from utils import ensure_dir, fasta_write, int_to_sequence, load_pickle, save_pickle, sequence_to_int


class TargetRegions:
    """Container for designed target-region records."""

    map = "*ACGTRYKMSWBDHVN-"

    def __init__(
        self,
        geneName: str = "",
        id: str = "",
        geneSequence: str | Sequence[int] | np.ndarray | None = None,
        startPos: Sequence[int] | np.ndarray | None = None,
        regionLength: Sequence[int] | np.ndarray | None = None,
        GC: Sequence[float] | np.ndarray | None = None,
        Tm: Sequence[float] | np.ndarray | None = None,
        specificity: Sequence[float] | np.ndarray | None = None,
        isoSpecificity: Sequence[float] | np.ndarray | None = None,
        penalties: Sequence[Sequence[float]] | np.ndarray | None = None,
        penaltyNames: Sequence[str] | None = None,
        sequence: Sequence[str] | None = None,
    ):
        self.geneName = geneName
        self.id = id
        self.sequence = list(sequence) if sequence is not None else []
        self.startPos = np.asarray(startPos if startPos is not None else [], dtype=int)
        self.regionLength = np.asarray(regionLength if regionLength is not None else [], dtype=int)
        self.GC = np.asarray(GC if GC is not None else [], dtype=float)
        self.Tm = np.asarray(Tm if Tm is not None else [], dtype=float)
        self.specificity = np.asarray(specificity if specificity is not None else [], dtype=float)
        self.isoSpecificity = np.asarray(isoSpecificity if isoSpecificity is not None else [], dtype=float)
        self.penalties = np.asarray(penalties if penalties is not None else [], dtype=float)
        self.penaltyNames = list(penaltyNames) if penaltyNames is not None else []
        self.numRegions = int(self.startPos.size)

        if geneSequence is not None and not self.sequence and self.numRegions:
            if isinstance(geneSequence, str):
                seq_text = geneSequence
            else:
                seq_text = int_to_sequence(geneSequence)
            self.sequence = []
            for start, length in zip(self.startPos, self.regionLength):
                start0 = int(start) - 1 if int(start) >= 1 else int(start)
                self.sequence.append(seq_text[start0 : start0 + int(length)])

    def fastawrite(self, filePath: str | Path, overwrite: bool = False) -> None:
        path = Path(filePath)
        if overwrite and path.exists():
            path.unlink()
        headers: List[str] = []
        for i in range(self.numRegions):
            parts = [
                f"id={self.id}",
                f"geneName={self.geneName}",
                f"startPos={int(self.startPos[i])}",
                f"regionLength={int(self.regionLength[i])}",
                f"GC={float(self.GC[i]) if self.GC.size else np.nan}",
                f"Tm={float(self.Tm[i]) if self.Tm.size else np.nan}",
                f"specificity={float(self.specificity[i]) if self.specificity.size else np.nan}",
                f"isoSpecificity={float(self.isoSpecificity[i]) if self.isoSpecificity.size else np.nan}",
            ]
            if self.penalties.size and self.penaltyNames:
                pen = np.asarray(self.penalties)
                if pen.ndim == 1:
                    pen = pen.reshape(1, -1)
                for row, name in enumerate(self.penaltyNames):
                    parts.append(f"p_{name}={float(pen[row, i])}")
            headers.append(" ".join(parts))
        fasta_write(path, headers, self.sequence, append=not overwrite)
        return None

    def Save(self, dirPath: str | Path) -> None:
        path = ensure_dir(dirPath)
        payload = self.to_dict()
        save_pickle(path / "TargetRegions.pkl", payload)
        for key, value in payload.items():
            should_save = value is not None
            if isinstance(value, str):
                should_save = value != ""
            elif isinstance(value, list):
                should_save = len(value) > 0
            elif isinstance(value, np.ndarray):
                should_save = value.size > 0
            if should_save:
                save_pickle(path / f"{key}.matb", value)
        return None

    def to_dict(self) -> dict[str, Any]:
        return {
            "geneName": self.geneName,
            "id": self.id,
            "sequence": self.sequence,
            "startPos": self.startPos,
            "regionLength": self.regionLength,
            "GC": self.GC,
            "Tm": self.Tm,
            "specificity": self.specificity,
            "isoSpecificity": self.isoSpecificity,
            "penalties": self.penalties,
            "penaltyNames": self.penaltyNames,
            "numRegions": self.numRegions,
        }

    @staticmethod
    def Load(dirPath: str | Path, verbose: bool = False) -> "TargetRegions":
        path = Path(dirPath)
        if not path.is_dir():
            raise ValueError("The provided path is not valid")
        payload_path = path / "TargetRegions.pkl"
        if payload_path.exists():
            payload = load_pickle(payload_path)
            payload.pop("numRegions", None)
            return TargetRegions(**payload)
        payload: dict[str, Any] = {}
        for field in ["geneName", "id", "sequence", "startPos", "regionLength", "GC", "Tm", "specificity", "isoSpecificity", "penalties", "penaltyNames"]:
            p = path / f"{field}.matb"
            if p.exists():
                payload[field] = load_pickle(p)
        return TargetRegions(**payload)
