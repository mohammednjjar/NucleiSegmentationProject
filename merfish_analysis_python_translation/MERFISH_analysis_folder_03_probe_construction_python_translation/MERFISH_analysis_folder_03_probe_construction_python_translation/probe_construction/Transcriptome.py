"""Python translation of MERFISH_analysis/probe_construction/Transcriptome.m."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterable, List, Sequence
import re
import numpy as np

try:
    from .utils import as_list, ensure_dir, fasta_read, int_to_sequence, load_pickle, save_pickle, sequence_to_int
except ImportError:
    from utils import as_list, ensure_dir, fasta_read, int_to_sequence, load_pickle, save_pickle, sequence_to_int


class Transcriptome:
    """Transcriptome container with sequence, gene-name, ID, abundance, and CDS indexing."""

    def __init__(
        self,
        transcriptome: Any = None,
        verbose: bool = False,
        headerType: str = "cufflinks",
        IDType: str = "",
        abundPath: str | None = None,
    ):
        self.verbose = bool(verbose)
        self.headerType = headerType
        self.IDType = IDType
        self.abundPath = abundPath
        self.transPath = ""
        self.abundLoaded = False
        self.numTranscripts = 0
        self.numGenes = 0
        self.ids: List[str] = []
        self.geneNames: List[str] = []
        self.intSequences: List[np.ndarray] = []
        self.abundance = np.asarray([], dtype=float)
        self.id2Ind: dict[str, int] = {}
        self.name2Ind: dict[str, List[int]] = {}
        self.cds = np.empty((0, 2), dtype=int)
        self.idVersion: List[str] = []

        if transcriptome is None:
            return

        if isinstance(transcriptome, str):
            path = Path(transcriptome)
            if not path.is_file():
                raise ValueError("Invalid path to target sequences")
            self.transPath = str(path)
            records = fasta_read(path)
            self._parse_fasta_records(records)
        elif isinstance(transcriptome, list) and len(transcriptome) >= 3 and not isinstance(transcriptome[0], dict):
            ids, names, seqs = transcriptome[0], transcriptome[1], transcriptome[2]
            abund = transcriptome[3] if len(transcriptome) > 3 else []
            cds = transcriptome[4] if len(transcriptome) > 4 else None
            versions = transcriptome[5] if len(transcriptome) > 5 else None
            self.ids = [str(x) for x in ids]
            self.geneNames = [str(x) for x in names]
            self.intSequences = [sequence_to_int(s) for s in seqs]
            if abund is not None and len(abund) > 0:
                self.abundance = np.asarray(abund, dtype=float).ravel()
                self.abundLoaded = True
            if cds is None:
                self.cds = -np.ones((len(self.ids), 2), dtype=int)
            else:
                self.cds = np.asarray(cds, dtype=int).reshape(len(self.ids), 2)
            if versions is None:
                self.idVersion = [""] * len(self.ids)
            else:
                self.idVersion = [str(v) for v in versions]
        elif isinstance(transcriptome, list) and (len(transcriptome) == 0 or isinstance(transcriptome[0], dict)):
            self._parse_fasta_records(transcriptome)
        elif isinstance(transcriptome, dict) and "Header" in transcriptome and "Sequence" in transcriptome:
            self._parse_fasta_records([transcriptome])
        else:
            raise ValueError("Invalid transcriptome input")

        if self.abundance.size == 0:
            if abundPath:
                self.AddAbundances(abundPath)
            else:
                self.abundance = np.ones(len(self.ids), dtype=float)
        self.UpdateInternalIndexing()

    def _parse_fasta_records(self, records: Sequence[dict[str, str]]) -> None:
        self.ids = []
        self.geneNames = []
        self.intSequences = []
        self.cds = -np.ones((len(records), 2), dtype=int)
        self.idVersion = [""] * len(records)
        for i, rec in enumerate(records):
            header = rec.get("Header", "")
            seq = rec.get("Sequence", "")
            if self.headerType == "cufflinks":
                id_match = re.match(r"(?P<id>\S+)", header)
                gene_match = re.search(r"gene=(\S+)", header)
                version_match = re.search(r"transcript_version=(\S+)", header)
                cds_match = re.search(r"CDS=(\d+)-(\d+)", header)
                self.ids.append(id_match.group("id") if id_match else str(i + 1))
                self.geneNames.append(gene_match.group(1) if gene_match else "")
                self.idVersion[i] = version_match.group(1) if version_match else ""
                if cds_match:
                    self.cds[i, :] = [int(cds_match.group(1)), int(cds_match.group(2))]
            elif self.headerType == "ensembl":
                gene_match = re.search(r"gene:(\S+)", header)
                self.ids.append(gene_match.group(1) if gene_match else header.split()[0])
                self.geneNames.append("")
            elif self.headerType == "custom":
                pieces = header.split()
                self.ids.append(pieces[0] if pieces else str(i + 1))
                self.geneNames.append(pieces[1] if len(pieces) > 1 else pieces[0] if pieces else "")
            else:
                raise ValueError("headerType must be cufflinks, ensembl, or custom")
            self.intSequences.append(sequence_to_int(seq))

    def AddEntries(self, names: Sequence[str], ids: Sequence[str], seqs: Sequence[str], abunds: Sequence[float], cds: Sequence[Sequence[int]] | None = None, idVersions: Sequence[str] | None = None) -> None:
        n = len(names)
        self.ids.extend([str(x) for x in ids])
        self.geneNames.extend([str(x) for x in names])
        self.intSequences.extend([sequence_to_int(s) for s in seqs])
        self.abundance = np.concatenate([self.abundance, np.asarray(abunds, dtype=float).ravel()])
        new_cds = -np.ones((n, 2), dtype=int) if cds is None else np.asarray(cds, dtype=int).reshape(n, 2)
        self.cds = np.vstack([self.cds, new_cds]) if self.cds.size else new_cds
        self.idVersion.extend([""] * n if idVersions is None else [str(v) for v in idVersions])
        self.UpdateInternalIndexing()
        return None

    def UpdateInternalIndexing(self) -> None:
        self.numTranscripts = len(self.ids)
        self.id2Ind = {id_value: i for i, id_value in enumerate(self.ids)}
        self.name2Ind = {}
        for i, name in enumerate(self.geneNames):
            self.name2Ind.setdefault(name, []).append(i)
        self.numGenes = len(self.name2Ind)
        if self.abundance.size == 0 and self.numTranscripts:
            self.abundance = np.ones(self.numTranscripts, dtype=float)
        if self.cds.size == 0 and self.numTranscripts:
            self.cds = -np.ones((self.numTranscripts, 2), dtype=int)
        if len(self.idVersion) < self.numTranscripts:
            self.idVersion.extend([""] * (self.numTranscripts - len(self.idVersion)))
        return None

    def AddAbundances(self, *args: Any) -> tuple[list[str], list[str]]:
        if len(args) == 1 and isinstance(args[0], str):
            path = Path(args[0])
            if not path.is_file():
                raise ValueError("Invalid path to abundance data")
            found_ids: List[str] = []
            fpkm: List[float] = []
            with path.open("r", encoding="utf-8") as handle:
                header = handle.readline().rstrip("\n").split("\t")
                id_col = 0
                fpkm_col = 9 if len(header) > 9 else len(header) - 1
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) <= max(id_col, fpkm_col):
                        continue
                    found_ids.append(parts[id_col])
                    try:
                        fpkm.append(float(parts[fpkm_col]))
                    except ValueError:
                        fpkm.append(0.0)
        elif len(args) == 2:
            found_ids = [str(x) for x in args[0]]
            fpkm = [float(x) for x in args[1]]
        else:
            raise ValueError("AddAbundances expects a file path or ids plus abundance vector")
        self.abundance = np.zeros(self.numTranscripts, dtype=float)
        found_map = {id_value: float(value) for id_value, value in zip(found_ids, fpkm)}
        common = []
        for i, id_value in enumerate(self.ids):
            if id_value in found_map:
                self.abundance[i] = found_map[id_value]
                common.append(id_value)
        self.abundLoaded = True
        common_set = set(common)
        return [x for x in found_ids if x not in common_set], [x for x in self.ids if x not in common_set]

    def GetAbundanceByID(self, ids: str | Sequence[str]) -> np.ndarray:
        ids_list = [ids] if isinstance(ids, str) else list(ids)
        values = np.full(len(ids_list), np.nan, dtype=float)
        for i, id_value in enumerate(ids_list):
            idx = self.id2Ind.get(str(id_value))
            if idx is not None:
                values[i] = self.abundance[idx]
        return values

    def GetAbundanceByName(self, names: str | Sequence[str], returnType: str = "isoformFunction", isoformFunc: Callable[[Sequence[float]], float] = sum) -> Any:
        names_list = [names] if isinstance(names, str) else list(names)
        all_values: List[np.ndarray] = []
        for name in names_list:
            inds = self.name2Ind.get(str(name), [])
            all_values.append(self.abundance[inds])
        if returnType == "allIsoforms":
            return all_values
        return np.asarray([float(isoformFunc(vals)) if len(vals) else np.nan for vals in all_values], dtype=float)

    def GetSequencesByName(self, names: str | Sequence[str]) -> Any:
        single = isinstance(names, str)
        names_list = [names] if single else list(names)
        out: List[Any] = []
        for name in names_list:
            inds = self.name2Ind.get(str(name), [])
            seqs = [int_to_sequence(self.intSequences[i]) for i in inds]
            out.append(seqs[0] if len(seqs) == 1 else seqs)
        return out[0] if single else out

    def GetSequenceByID(self, ids: str | Sequence[str]) -> Any:
        single = isinstance(ids, str)
        ids_list = [ids] if single else list(ids)
        out: List[str] = []
        for id_value in ids_list:
            idx = self.id2Ind.get(str(id_value))
            out.append(int_to_sequence(self.intSequences[idx]) if idx is not None else "")
        return out[0] if single else out

    def CDSByID(self, ids: str | Sequence[str]) -> np.ndarray:
        ids_list = [ids] if isinstance(ids, str) else list(ids)
        out = -np.ones((len(ids_list), 2), dtype=int)
        for i, id_value in enumerate(ids_list):
            idx = self.id2Ind.get(str(id_value))
            if idx is not None:
                out[i, :] = self.cds[idx, :]
        return out

    def GetIDVersion(self, ids: str | Sequence[str]) -> List[str]:
        ids_list = [ids] if isinstance(ids, str) else list(ids)
        out: List[str] = []
        for id_value in ids_list:
            idx = self.id2Ind.get(str(id_value))
            out.append(self.idVersion[idx] if idx is not None else "")
        return out

    def GetIDsByName(self, names: str | Sequence[str]) -> Any:
        single = isinstance(names, str)
        names_list = [names] if single else list(names)
        out = [[self.ids[i] for i in self.name2Ind.get(str(name), [])] for name in names_list]
        return out[0] if single else out

    def GetNameById(self, id: str) -> str:
        idx = self.id2Ind.get(str(id))
        return self.geneNames[idx] if idx is not None else ""

    def GetNames(self) -> List[str]:
        return list(self.name2Ind.keys())

    def Slice(self, geneID: Sequence[str] | None = None, geneName: Sequence[str] | None = None) -> "Transcriptome":
        if geneID:
            inds = self.GetInternalInds(geneID=list(geneID))[0]
        elif geneName:
            inds = self.GetInternalInds(geneName=list(geneName))[0]
        else:
            inds = list(range(self.numTranscripts))
        return Transcriptome([
            [self.ids[i] for i in inds],
            [self.geneNames[i] for i in inds],
            [self.intSequences[i] for i in inds],
            self.abundance[inds],
            self.cds[inds, :],
            [self.idVersion[i] for i in inds],
        ], verbose=self.verbose, headerType=self.headerType, IDType=self.IDType, abundPath=self.abundPath)

    def Save(self, dirPath: str | Path) -> None:
        path = ensure_dir(dirPath)
        payload = {
            "verbose": self.verbose,
            "headerType": self.headerType,
            "IDType": self.IDType,
            "abundPath": self.abundPath,
            "transPath": self.transPath,
            "abundLoaded": self.abundLoaded,
            "numTranscripts": self.numTranscripts,
            "numGenes": self.numGenes,
            "ids": self.ids,
            "geneNames": self.geneNames,
            "intSequences": self.intSequences,
            "abundance": self.abundance,
            "cds": self.cds,
            "idVersion": self.idVersion,
        }
        save_pickle(path / "Transcriptome.pkl", payload)
        for key, value in payload.items():
            save_pickle(path / f"{key}.matb", value)
        return None

    def GetInternalInds(self, geneName: Sequence[str] | None = None, geneID: Sequence[str] | None = None, inds: Sequence[int] | None = None) -> tuple[list[int], list[str], list[str], list[bool]]:
        if inds is not None and len(inds) > 0:
            out_inds = [int(i) for i in inds]
            return out_inds, [self.ids[i] for i in out_inds], [self.geneNames[i] for i in out_inds], [True] * len(out_inds)
        if geneName is not None and len(geneName) > 0:
            out_inds: List[int] = []
            valid: List[bool] = []
            for name in geneName:
                local = self.name2Ind.get(str(name), [])
                valid.append(bool(local))
                out_inds.extend(local)
            return out_inds, [self.ids[i] for i in out_inds], [self.geneNames[i] for i in out_inds], valid
        if geneID is not None and len(geneID) > 0:
            out_inds = []
            valid = []
            for id_value in geneID:
                idx = self.id2Ind.get(str(id_value))
                valid.append(idx is not None)
                if idx is not None:
                    out_inds.append(idx)
            return out_inds, [self.ids[i] for i in out_inds], [self.geneNames[i] for i in out_inds], valid
        all_inds = list(range(self.numTranscripts))
        return all_inds, self.ids.copy(), self.geneNames.copy(), [True] * self.numTranscripts

    @staticmethod
    def Load(dirPath: str | Path, verbose: bool = False) -> "Transcriptome":
        path = Path(dirPath)
        if not path.is_dir():
            raise ValueError("The provided path is not valid")
        payload_path = path / "Transcriptome.pkl"
        if payload_path.exists():
            payload = load_pickle(payload_path)
        else:
            fields = ["ids", "geneNames", "intSequences", "abundance", "cds", "idVersion", "headerType", "IDType", "abundPath"]
            payload = {field: load_pickle(path / f"{field}.matb") for field in fields if (path / f"{field}.matb").exists()}
        obj = Transcriptome()
        for key, value in payload.items():
            setattr(obj, key, value)
        obj.verbose = verbose
        obj.UpdateInternalIndexing()
        return obj
