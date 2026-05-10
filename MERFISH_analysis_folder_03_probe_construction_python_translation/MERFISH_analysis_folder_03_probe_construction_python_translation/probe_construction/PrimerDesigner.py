"""Python translation of MERFISH_analysis/probe_construction/PrimerDesigner.m."""
from __future__ import annotations

from pathlib import Path
from typing import Any, List, Sequence
import uuid
import numpy as np

try:
    from .OTTable import OTTable
    from .TRDesigner import TRDesigner
    from .utils import ensure_dir, fasta_write, int_to_sequence, load_pickle, reverse_complement_int, rolling_hash, save_pickle, sequence_to_int
except ImportError:
    from OTTable import OTTable
    from TRDesigner import TRDesigner
    from utils import ensure_dir, fasta_write, int_to_sequence, load_pickle, reverse_complement_int, rolling_hash, save_pickle, sequence_to_int


class PrimerDesigner:
    """Orthogonal primer designer translated from MATLAB PrimerDesigner.m."""

    def __init__(
        self,
        verbose: bool = True,
        ntComposition: Sequence[float] | None = None,
        OTTables: Sequence[OTTable] | None = None,
        OTTableNames: Sequence[str] | None = None,
        parallel: Any = None,
        seqs: Sequence[str | Sequence[int]] | np.ndarray | None = None,
        primerLength: int = 20,
        numPrimersToGenerate: int = 0,
        homologyMax: int = 8,
        monovalentSalt: float = 0.3,
        primerConc: float = 0.5e-6,
        seqsToRemove: Sequence[str] | None = None,
    ):
        self.verbose = bool(verbose)
        self.ntComposition = np.asarray(ntComposition if ntComposition is not None else [0.25, 0.25, 0.25, 0.25], dtype=float)
        self.ntComposition = self.ntComposition / self.ntComposition.sum()
        self.OTTables = list(OTTables or [])
        self.OTTableNames = list(OTTableNames or [])
        self.monovalentSalt = float(monovalentSalt)
        self.primerConc = float(primerConc)
        self.seqsToRemove = list(seqsToRemove or ["AAAA", "TTTT", "GGGG", "CCCC"])
        self.numPrimers = 0
        self.primerLength = int(primerLength)
        self.seqs = np.empty((0, self.primerLength), dtype=np.int16)
        self.parallel = parallel
        self.numPar = int(getattr(parallel, "NumWorkers", 0) or 0)
        self.homologyMax = int(homologyMax)
        self.gc = np.asarray([], dtype=float)
        self.Tm = np.asarray([], dtype=float)
        self.penalties = np.empty((0, len(self.OTTables)), dtype=float)
        self.homologyMat = np.empty((0, 0), dtype=np.int8)
        self.seqHash = np.empty((0, 0), dtype=np.int64)
        self.seqRCHash = np.empty((0, 0), dtype=np.int64)

        if seqs is not None and len(seqs) > 0:
            if isinstance(seqs, np.ndarray) and seqs.ndim == 2 and np.issubdtype(seqs.dtype, np.integer):
                arr = seqs.astype(np.int16)
            else:
                arr = np.vstack([sequence_to_int(s, acgt_only=True) for s in seqs]).astype(np.int16)
            if np.any((arr < 0) | (arr > 3)):
                raise ValueError("The sequences must contain only A/C/G/T or integer values 0..3")
            self.seqs = arr
            self.primerLength = arr.shape[1]
            self.numPrimers = arr.shape[0]
            self.CalculatePrimerProperties()
        elif numPrimersToGenerate > 0:
            self.AddRandomSequences(ntComposition=self.ntComposition, primerLength=self.primerLength, numPrimersToGenerate=int(numPrimersToGenerate))

    def AddRandomSequences(self, ntComposition: Sequence[float] | None = None, primerLength: int | None = None, numPrimersToGenerate: int = 1_000_000) -> None:
        comp = np.asarray(ntComposition if ntComposition is not None else self.ntComposition, dtype=float)
        if comp.size != 4:
            raise ValueError("ntComposition must have four entries")
        comp = comp / comp.sum()
        primer_len = int(primerLength or self.primerLength)
        count = int(numPrimersToGenerate)
        rng = np.random.default_rng()
        seq = rng.choice(np.arange(4, dtype=np.int16), size=(count, primer_len), p=comp).astype(np.int16)
        if self.seqs.size and self.seqs.shape[1] != primer_len:
            raise ValueError("New primer length must match existing primer length")
        self.seqs = np.vstack([self.seqs, seq]) if self.seqs.size else seq
        self.primerLength = primer_len
        self.numPrimers = self.seqs.shape[0]
        self.CalculatePrimerProperties()
        return None

    def AddPrimer(self, seq: str | Sequence[int] | np.ndarray) -> None:
        int_seq = sequence_to_int(seq, acgt_only=True)
        if int_seq.size != self.primerLength:
            raise ValueError("The provided sequence must match primerLength")
        if np.any((int_seq < 0) | (int_seq > 3)):
            raise ValueError("All elements in integer sequences must be 0, 1, 2, or 3")
        self.seqs = np.vstack([self.seqs, int_seq.reshape(1, -1)]) if self.seqs.size else int_seq.reshape(1, -1)
        self.numPrimers = self.seqs.shape[0]
        self.CalculatePrimerProperties()
        return None

    def CalculatePrimerProperties(self, monovalentSalt: float | None = None, primerConc: float | None = None) -> None:
        if monovalentSalt is not None:
            self.monovalentSalt = float(monovalentSalt)
        if primerConc is not None:
            self.primerConc = float(primerConc)
        if self.numPrimers == 0:
            self.gc = np.asarray([], dtype=float)
            self.Tm = np.asarray([], dtype=float)
            self.penalties = np.empty((0, len(self.OTTables)), dtype=float)
            return None
        self.gc = np.mean((self.seqs == 1) | (self.seqs == 2), axis=1).astype(float)
        tm_values = []
        for row in self.seqs:
            dg = TRDesigner.SantaLuciaNearestNeighbor(row)
            totals = np.nansum(dg, axis=1)
            five_at = row[0] in (0, 3)
            three_at = row[-1] in (0, 3)
            h = totals[0] + 0.2 + 2.2 * five_at + 2.2 * three_at
            s = totals[1] - 5.7 + 6.9 * five_at + 6.9 * three_at
            s = s + 0.368 * (self.primerLength - 1) * np.log(self.monovalentSalt)
            tm_values.append(h * 1000.0 / (s + 1.9872 * np.log(self.primerConc)) - 273.15)
        self.Tm = np.asarray(tm_values, dtype=float)
        self.penalties = np.full((self.numPrimers, len(self.OTTables)), np.nan, dtype=float)
        for col, table in enumerate(self.OTTables):
            for row_id, row in enumerate(self.seqs):
                fwd = np.nansum(table.CalculatePenalty(row)[0])
                rev = np.nansum(table.CalculatePenalty(reverse_complement_int(row))[0])
                self.penalties[row_id, col] = fwd + rev
        return None

    def CutPrimers(self, Tm: Sequence[float] | None = None, GC: Sequence[float] | None = None, OTTables: Sequence[Any] | None = None) -> np.ndarray:
        keep = np.ones(self.numPrimers, dtype=bool)
        if GC is not None:
            keep = keep & (self.gc >= GC[0]) & (self.gc <= GC[1])
        if Tm is not None:
            keep = keep & (self.Tm >= Tm[0]) & (self.Tm <= Tm[1])
        filters = list(OTTables or [])
        for i in range(0, len(filters), 2):
            name = filters[i]
            bounds = filters[i + 1]
            if name not in self.OTTableNames:
                raise ValueError("Unrecognized OTTable name")
            col = self.OTTableNames.index(name)
            keep = keep & (self.penalties[:, col] >= bounds[0]) & (self.penalties[:, col] <= bounds[1])
        self._apply_keep(keep)
        return keep

    def RemoveForbiddenSeqs(self, seqsToRemove: Sequence[str] | None = None) -> np.ndarray:
        if seqsToRemove:
            self.seqsToRemove = list(seqsToRemove)
        keep = np.ones(self.numPrimers, dtype=bool)
        for forbidden in self.seqsToRemove:
            fseq = sequence_to_int(forbidden, acgt_only=True)
            if np.any((fseq < 0) | (fseq > 3)):
                raise ValueError("Invalid forbidden sequence")
            f_hash, f_valid = rolling_hash(fseq, fseq.size)
            target_hash = f_hash[0] if f_hash.size and f_valid[0] else -2
            for i, seq in enumerate(self.seqs):
                hashes, valid = rolling_hash(seq, fseq.size)
                keep[i] = keep[i] and not bool(np.any((hashes == target_hash) & valid))
        self._apply_keep(keep)
        return keep

    def RemoveSelfCompPrimers(self, homologyMax: int | None = None) -> np.ndarray:
        if homologyMax is not None:
            self.homologyMax = int(homologyMax)
        keep = np.ones(self.numPrimers, dtype=bool)
        seq_hash = []
        rc_hash = []
        for i, seq in enumerate(self.seqs):
            h, v = rolling_hash(seq, self.homologyMax)
            rh, rv = rolling_hash(reverse_complement_int(seq), self.homologyMax)
            valid_h = set(h[v].tolist())
            valid_rh = set(rh[rv].tolist())
            keep[i] = len(valid_h.intersection(valid_rh)) == 0
            seq_hash.append(h)
            rc_hash.append(rh)
        max_len = max((x.size for x in seq_hash), default=0)
        self.seqHash = np.full((len(seq_hash), max_len), -1, dtype=np.int64)
        self.seqRCHash = np.full((len(rc_hash), max_len), -1, dtype=np.int64)
        for i, h in enumerate(seq_hash):
            self.seqHash[i, : h.size] = h
        for i, h in enumerate(rc_hash):
            self.seqRCHash[i, : h.size] = h
        self._apply_keep(keep)
        return keep

    def RemoveHomologousPrimers(self, homologyMax: int | None = None) -> np.ndarray:
        if homologyMax is not None:
            self.homologyMax = int(homologyMax)
        hash_sets: List[set[int]] = []
        for seq in self.seqs:
            h, v = rolling_hash(seq, self.homologyMax)
            rh, rv = rolling_hash(reverse_complement_int(seq), self.homologyMax)
            hash_sets.append(set(h[v].tolist()).union(set(rh[rv].tolist())))
        n = self.numPrimers
        homology = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i + 1, n):
                if hash_sets[i].intersection(hash_sets[j]):
                    homology[i, j] = 1
                    homology[j, i] = 1
        self.homologyMat = homology
        active = np.ones(n, dtype=bool)
        while True:
            row_sum = homology[active][:, active].sum(axis=1)
            if row_sum.size == 0 or np.max(row_sum) == 0:
                break
            active_indices = np.where(active)[0]
            remove_id = active_indices[int(np.argmax(row_sum))]
            active[remove_id] = False
        self._apply_keep(active)
        return active

    def WriteFasta(self, filePath: str | Path, namePrefix: str = "", fieldPad: int | None = None) -> None:
        prefix = namePrefix or str(uuid.uuid4())[:8]
        pad = int(fieldPad if fieldPad is not None else max(1, int(np.ceil(np.log10(max(self.numPrimers, 1) + 1)))))
        headers = []
        seqs = []
        for i, seq in enumerate(self.seqs, start=1):
            headers.append(f"{prefix}-{i:0{pad}d} Tm={self.Tm[i-1]} GC={self.gc[i-1]}")
            seqs.append(int_to_sequence(seq))
        fasta_write(filePath, headers, seqs, append=False)
        return None

    def SetParallel(self, p: Any) -> None:
        self.parallel = p
        self.numPar = int(getattr(p, "NumWorkers", 0) or 0)
        return None

    def Save(self, dirPath: str | Path) -> None:
        path = ensure_dir(dirPath)
        payload = {key: value for key, value in self.__dict__.items() if key not in {"parallel"}}
        save_pickle(path / "PrimerDesigner.pkl", payload)
        for i, table in enumerate(self.OTTables, start=1):
            table.Save(path / f"OTTable_{i}")
        return None

    @staticmethod
    def Load(dirPath: str | Path) -> "PrimerDesigner":
        path = Path(dirPath)
        if not path.is_dir():
            raise ValueError("The provided path is not valid")
        payload = load_pickle(path / "PrimerDesigner.pkl")
        obj = PrimerDesigner(numPrimersToGenerate=0)
        obj.__dict__.update(payload)
        return obj

    def _apply_keep(self, keep: np.ndarray) -> None:
        self.seqs = self.seqs[keep, :]
        self.numPrimers = self.seqs.shape[0]
        self.gc = self.gc[keep] if self.gc.size else self.gc
        self.Tm = self.Tm[keep] if self.Tm.size else self.Tm
        self.penalties = self.penalties[keep, :] if self.penalties.size else self.penalties
        return None
