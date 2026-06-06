"""Python translation of MERFISH_analysis/probe_construction/TRDesigner.m."""
from __future__ import annotations

from pathlib import Path
from typing import Any, List, Sequence
import numpy as np

try:
    from .OTTable import OTTable
    from .TargetRegions import TargetRegions
    from .Transcriptome import Transcriptome
    from .utils import ensure_dir, int_to_sequence, load_pickle, rolling_hash, save_pickle, sequence_to_int, sliding_mean, sliding_sum
except ImportError:
    from OTTable import OTTable
    from TargetRegions import TargetRegions
    from Transcriptome import Transcriptome
    from utils import ensure_dir, int_to_sequence, load_pickle, rolling_hash, save_pickle, sequence_to_int, sliding_mean, sliding_sum


class TRDesigner:
    """Target-region designer translated from MATLAB TRDesigner.m."""

    def __init__(
        self,
        transcriptome: Transcriptome | None = None,
        OTTables: Sequence[OTTable] | None = None,
        OTTableNames: Sequence[str] | None = None,
        parallel: Any = None,
        verbose: bool = True,
        forbiddenSeqs: Sequence[str | Sequence[int]] | None = None,
        specificityTable: OTTable | None = None,
        isoSpecificityTables: Sequence[OTTable] | None = None,
        alwaysDefinedSpecificity: bool = False,
    ):
        self.verbose = bool(verbose)
        self.transcriptome = transcriptome
        self.parallel = parallel
        self.numPar = int(getattr(parallel, "NumWorkers", 0) or 0)
        self.dG: List[np.ndarray] = []
        self.gc: List[np.ndarray] = []
        self.isValid: List[np.ndarray] = []
        self.OTTableNames: List[str] = []
        self.OTTables: List[OTTable] = []
        self.penalties: List[List[np.ndarray]] = []
        self.forbiddenSeqs: List[np.ndarray] = []
        self.isForbiddenSeq: List[List[np.ndarray]] = []
        self.specificity: List[np.ndarray] = []
        self.specificityTable = specificityTable
        self.isoSpecificity: List[np.ndarray] = []
        self.isoSpecificityTables: List[OTTable] = list(isoSpecificityTables or [])
        self.alwaysDefinedSpecificity = bool(alwaysDefinedSpecificity)

        if transcriptome is None:
            return
        if not isinstance(transcriptome, Transcriptome):
            raise TypeError("transcriptome must be a Transcriptome object")

        self.gc = [(seq == 1) | (seq == 2) for seq in transcriptome.intSequences]
        self.isValid = [(seq >= 0) & (seq <= 3) for seq in transcriptome.intSequences]
        self.dG = [TRDesigner.SantaLuciaNearestNeighbor(seq) for seq in transcriptome.intSequences]

        tables = list(OTTables or [])
        names = list(OTTableNames or [])
        if len(tables) != len(names):
            raise ValueError("An equal number of names and OT tables must be provided")
        for table, name in zip(tables, names):
            self.AddOTTable(table, name)
        if specificityTable is not None:
            self.AddSpecificityTable(specificityTable, self.isoSpecificityTables)
        for seq in forbiddenSeqs or []:
            self.AddForbiddenSeq(seq)

    def AddForbiddenSeq(self, seq: str | Sequence[int] | np.ndarray, replace: bool = False) -> None:
        if replace:
            self.forbiddenSeqs = []
            self.isForbiddenSeq = []
        int_seq = sequence_to_int(seq, acgt_only=True)
        if int_seq.size == 0 or np.any((int_seq < 0) | (int_seq > 3)):
            raise ValueError("The provided sequence is invalid")
        self.forbiddenSeqs.append(int_seq)
        hits: List[np.ndarray] = []
        for local_seq in self.transcriptome.intSequences:
            hashes, valid = rolling_hash(local_seq, int_seq.size)
            f_hash, f_valid = rolling_hash(int_seq, int_seq.size)
            forbidden_hash = f_hash[0] if f_hash.size and f_valid[0] else -2
            hits.append((hashes == forbidden_hash) & valid)
        self.isForbiddenSeq.append(hits)
        return None

    def AddOTTable(self, otTable: OTTable, tableName: str, replace: bool = False) -> None:
        if not isinstance(otTable, OTTable) or not isinstance(tableName, str):
            raise TypeError("Both an OTTable and a name must be provided")
        if replace:
            self.OTTables = []
            self.OTTableNames = []
            self.penalties = []
        self.OTTables.append(otTable)
        self.OTTableNames.append(tableName)
        penalty = [otTable.CalculatePenalty(seq)[0] for seq in self.transcriptome.intSequences]
        self.penalties.append(penalty)
        return None

    def AddSpecificityTable(self, specificityTable: OTTable, isoSpecificityTables: Sequence[OTTable] | None = None) -> None:
        if not isinstance(specificityTable, OTTable):
            raise TypeError("The provided specificity table is not a valid OTTable")
        self.specificityTable = specificityTable
        self.isoSpecificityTables = list(isoSpecificityTables or [])
        if self.specificityTable.uniformWeight:
            normalization = np.ones(self.transcriptome.numTranscripts, dtype=float)
        else:
            normalization = self.transcriptome.abundance
        self.specificity = []
        self.isoSpecificity = []
        if not self.isoSpecificityTables:
            for s, seq in enumerate(self.transcriptome.intSequences):
                penalty = specificityTable.CalculatePenalty(seq)[0]
                with np.errstate(divide="ignore", invalid="ignore"):
                    spec = normalization[s] / penalty
                self.specificity.append(spec)
                self.isoSpecificity.append(np.ones_like(spec, dtype=float))
        else:
            table_by_name = {table.name: table for table in self.isoSpecificityTables}
            for s, seq in enumerate(self.transcriptome.intSequences):
                name = self.transcriptome.geneNames[s]
                local_table = table_by_name.get(name)
                if local_table is None:
                    self.specificity.append(np.full(max(0, len(seq) - specificityTable.seedLength + 1), np.nan))
                    self.isoSpecificity.append(np.full(max(0, len(seq) - specificityTable.seedLength + 1), np.nan))
                    continue
                iso_counts = local_table.CalculatePenalty(seq)[0]
                total_counts = specificityTable.CalculatePenalty(seq)[0]
                with np.errstate(divide="ignore", invalid="ignore"):
                    self.isoSpecificity.append(normalization[s] / iso_counts)
                    self.specificity.append(iso_counts / total_counts)
        if self.alwaysDefinedSpecificity:
            self.specificity = [np.where(np.isnan(x), 1.0, x) for x in self.specificity]
        return None

    def SetParallel(self, p: Any) -> None:
        self.parallel = p
        self.numPar = int(getattr(p, "NumWorkers", 0) or 0)
        return None

    def GetRegionForbiddenSeqs(self, len: int, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, inds: Sequence[int] | None = None) -> tuple[list[np.ndarray], list[str], list[str]]:
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID, inds=inds)
        out: List[np.ndarray] = []
        for idx in indices:
            n = max(0, self.transcriptome.intSequences[idx].size - len + 1)
            keep = np.ones(n, dtype=bool)
            for fseq, hits_all in zip(self.forbiddenSeqs, self.isForbiddenSeq):
                window_len = len - fseq.size + 1
                if window_len <= 0:
                    continue
                hits = hits_all[idx]
                has = sliding_sum(hits.astype(float), window_len) > 0
                keep = keep & ~has[:n]
            out.append(keep)
        return out, ids, names

    def GetRegionValidity(self, len: int, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, inds: Sequence[int] | None = None) -> tuple[list[np.ndarray], list[str], list[str]]:
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID, inds=inds)
        out = [sliding_sum(((self.transcriptome.intSequences[i] > 3) | (self.transcriptome.intSequences[i] < 0)).astype(float), len) == 0 for i in indices]
        return out, ids, names

    def GetRegionGC(self, len: int, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, inds: Sequence[int] | None = None) -> tuple[list[np.ndarray], list[str], list[str]]:
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID, inds=inds)
        out = [sliding_mean(self.gc[i].astype(float), len) for i in indices]
        return out, ids, names

    def GetRegionPenalty(self, len: int, OTtableName: str, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, inds: Sequence[int] | None = None) -> tuple[list[np.ndarray], list[str], list[str]]:
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID, inds=inds)
        if OTtableName not in self.OTTableNames:
            raise ValueError("The specified table does not exist")
        table_id = self.OTTableNames.index(OTtableName)
        window_len = len - self.OTTables[table_id].seedLength + 1
        out = [sliding_sum(np.nan_to_num(self.penalties[table_id][i], nan=0.0), window_len) for i in indices]
        return out, ids, names

    def GetRegionSpecificity(self, len: int, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, inds: Sequence[int] | None = None) -> tuple[list[np.ndarray], list[str], list[str]]:
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID, inds=inds)
        window_len = len - self.specificityTable.seedLength + 1
        out = [sliding_mean(self.specificity[i], window_len) for i in indices]
        return out, ids, names

    def GetRegionIsoSpecificity(self, len: int, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, inds: Sequence[int] | None = None) -> tuple[list[np.ndarray], list[str], list[str]]:
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID, inds=inds)
        if self.isoSpecificityTables:
            seed_len = self.isoSpecificityTables[0].seedLength
        elif self.specificityTable is not None:
            seed_len = self.specificityTable.seedLength
        else:
            seed_len = len
        window_len = len - seed_len + 1
        out = [sliding_mean(self.isoSpecificity[i], window_len) if i < self.isoSpecificity.__len__() else np.ones(max(0, self.transcriptome.intSequences[i].size - len + 1)) for i in indices]
        return out, ids, names

    def GetRegionTm(self, len: int, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, monovalentSalt: float = 0.3, probeConc: float = 5e-9, inds: Sequence[int] | None = None) -> tuple[list[np.ndarray], list[str], list[str]]:
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID, inds=inds)
        out: List[np.ndarray] = []
        for idx in indices:
            int_seq = self.transcriptome.intSequences[idx]
            dg = self.dG[idx]
            if int_seq.size < len:
                out.append(np.asarray([], dtype=float))
                continue
            h = sliding_sum(dg[0, :], len - 1)
            s = sliding_sum(dg[1, :], len - 1)
            five_prime_at = (int_seq[: int_seq.size - len + 1] == 0) | (int_seq[: int_seq.size - len + 1] == 3)
            three_prime_at = (int_seq[len - 1 :] == 0) | (int_seq[len - 1 :] == 3)
            h = h + 0.2 + 2.2 * five_prime_at + 2.2 * three_prime_at
            s = s - 5.7 + 6.9 * five_prime_at + 6.9 * three_prime_at
            s = s + 0.368 * (len - 1) * np.log(monovalentSalt)
            tm = h * 1000.0 / (s + 1.9872 * np.log(probeConc)) - 273.15
            out.append(tm)
        return out, ids, names

    def DesignTargetRegions(
        self,
        geneName: Sequence[str] | None = None,
        geneID: Sequence[str] | None = None,
        regionLength: int | Sequence[int] = 30,
        Tm: Sequence[float] | None = None,
        GC: Sequence[float] | None = None,
        specificity: Sequence[float] | None = None,
        isoSpecificity: Sequence[float] | None = None,
        monovalentSalt: float = 0.3,
        probeConc: float = 5e-9,
        OTTables: Sequence[Any] | None = None,
        includeSequence: bool = True,
        threePrimeSpace: int = 0,
        removeForbiddenSeqs: bool = False,
    ) -> list[TargetRegions]:
        lengths = [int(regionLength)] if isinstance(regionLength, int) else [int(x) for x in regionLength]
        indices, ids, names, _ = self.transcriptome.GetInternalInds(geneName=geneName, geneID=geneID)
        table_filters = list(OTTables or [])
        results: List[TargetRegions] = []
        for idx, id_value, name in zip(indices, ids, names):
            region_props = []
            for length in lengths:
                valid = self.GetRegionValidity(length, inds=[idx])[0][0]
                gc = self.GetRegionGC(length, inds=[idx])[0][0]
                tm = self.GetRegionTm(length, monovalentSalt=monovalentSalt, probeConc=probeConc, inds=[idx])[0][0]
                spec = self.GetRegionSpecificity(length, inds=[idx])[0][0] if self.specificityTable is not None else np.ones_like(gc)
                iso = self.GetRegionIsoSpecificity(length, inds=[idx])[0][0] if (self.specificityTable is not None or self.isoSpecificity) else np.ones_like(gc)
                keep = valid.copy()
                if removeForbiddenSeqs and self.forbiddenSeqs:
                    keep = keep & self.GetRegionForbiddenSeqs(length, inds=[idx])[0][0]
                pad = 10 * np.finfo(float).eps
                if GC is not None:
                    keep = keep & (gc >= GC[0] - pad) & (gc <= GC[1] + pad)
                if Tm is not None:
                    keep = keep & (tm >= Tm[0] - pad) & (tm <= Tm[1] + pad)
                if specificity is not None:
                    keep = keep & (spec >= specificity[0] - pad) & (spec <= specificity[1] + pad)
                if isoSpecificity is not None:
                    keep = keep & (iso >= isoSpecificity[0] - pad) & (iso <= isoSpecificity[1] + pad)
                if table_filters:
                    for j in range(0, len(table_filters), 2):
                        table_name = table_filters[j]
                        bounds = table_filters[j + 1]
                        pen = self.GetRegionPenalty(length, str(table_name), inds=[idx])[0][0]
                        keep = keep & (pen >= bounds[0] - pad) & (pen <= bounds[1] + pad)
                starts = np.where(keep)[0] + 1
                for start in starts:
                    start0 = int(start) - 1
                    region_props.append([start, length, tm[start0], gc[start0], spec[start0], iso[start0]])
            selected = TRDesigner.TileRegions(np.asarray(region_props, dtype=float).T if region_props else np.zeros((6, 0)), threePrimeSpace)
            results.append(TargetRegions(
                id=id_value,
                geneName=name,
                geneSequence=self.transcriptome.intSequences[idx] if includeSequence else None,
                startPos=selected[0, :].astype(int) if selected.size else [],
                regionLength=selected[1, :].astype(int) if selected.size else [],
                Tm=selected[2, :] if selected.size else [],
                GC=selected[3, :] if selected.size else [],
                specificity=selected[4, :] if selected.size else [],
                isoSpecificity=selected[5, :] if selected.size else [],
            ))
        return results

    def Save(self, dirPath: str | Path) -> None:
        path = ensure_dir(dirPath)
        payload = {key: value for key, value in self.__dict__.items() if key not in {"parallel"}}
        save_pickle(path / "TRDesigner.pkl", payload)
        if self.transcriptome is not None:
            self.transcriptome.Save(path / "transcriptome")
        for i, table in enumerate(self.OTTables, start=1):
            table.Save(path / f"OTTable_{i}")
        if self.specificityTable is not None:
            self.specificityTable.Save(path / "specificityTable")
        return None

    @staticmethod
    def Load(dirPath: str | Path, lightweight: bool = False) -> "TRDesigner":
        path = Path(dirPath)
        if not path.is_dir():
            raise ValueError("The provided path is not valid")
        if (path / "TRDesigner.pkl").exists():
            payload = load_pickle(path / "TRDesigner.pkl")
            obj = TRDesigner()
            obj.__dict__.update(payload)
            return obj
        transcriptome = Transcriptome.Load(path / "transcriptome") if (path / "transcriptome").is_dir() else None
        return TRDesigner(transcriptome=transcriptome)

    @staticmethod
    def SantaLuciaNearestNeighbor(intSeq: str | Sequence[int] | np.ndarray) -> np.ndarray:
        seq = sequence_to_int(intSeq, acgt_only=False)
        if seq.size < 2:
            return np.empty((2, 0), dtype=float)
        dg = np.full((2, seq.size - 1), np.nan, dtype=float)
        nn_id = 4 * seq[:-1] + seq[1:]
        valid = (seq[:-1] >= 0) & (seq[:-1] <= 3) & (seq[1:] >= 0) & (seq[1:] <= 3)
        h = np.asarray([-7.6, -8.4, -7.8, -7.2, -8.5, -8.0, -10.6, -7.8, -8.2, -9.8, -8.0, -8.4, -7.2, -8.2, -8.5, -7.6], dtype=float)
        s = np.asarray([-21.3, -22.4, -21.0, -20.4, -22.7, -19.9, -27.2, -21.0, -22.2, -24.4, -19.9, -22.4, -21.3, -22.2, -22.7, -21.3], dtype=float)
        dg[0, valid] = h[nn_id[valid].astype(int)]
        dg[1, valid] = s[nn_id[valid].astype(int)]
        return dg

    @staticmethod
    def TileRegions(regionProps: Sequence[Sequence[float]] | np.ndarray, padLength: int) -> np.ndarray:
        arr = np.asarray(regionProps, dtype=float)
        if arr.size == 0:
            return np.zeros((6, 0), dtype=float)
        if arr.ndim != 2 or arr.shape[0] != 6:
            raise ValueError("regionProps must be a 6xN array")
        order = np.argsort(arr[0, :])
        arr = arr[:, order]
        start_pos = arr[0, :]
        next_available = start_pos + arr[1, :] + int(padLength)
        keep = [0]
        while True:
            current_next = next_available[keep[-1]]
            candidates = np.where(start_pos >= current_next)[0]
            if candidates.size == 0:
                break
            keep.append(int(candidates[0]))
        return arr[:, keep]
