"""Python translation of `example_scripts/library_design_example.m`.

Purpose: demonstrate MERFISH probe-library design by constructing target
regions, specificity/off-target tables, readout assignments, encoding probes,
primers, and final oligo FASTA files.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable
import csv
import random
import re
import shutil
import time

import numpy as np

try:
    from fileIO.LoadCodebook import LoadCodebook
except Exception:
    LoadCodebook = None

try:
    from probe_construction.Transcriptome import Transcriptome
    from probe_construction.OTTable import OTTable
    from probe_construction.TRDesigner import TRDesigner
    from probe_construction.TargetRegions import TargetRegions
    from probe_construction.PrimerDesigner import PrimerDesigner
except Exception:
    Transcriptome = None
    OTTable = None
    TRDesigner = None
    TargetRegions = None
    PrimerDesigner = None

try:
    from .script_utils import (
        FastaRecord,
        ensure_dir,
        import_symbol,
        page_break,
        read_fasta,
        reverse_complement,
        set_figure_save_path,
        strip_spaces,
        tic,
        toc,
        write_fasta,
    )
except ImportError:
    from script_utils import (
        FastaRecord,
        ensure_dir,
        import_symbol,
        page_break,
        read_fasta,
        reverse_complement,
        set_figure_save_path,
        strip_spaces,
        tic,
        toc,
        write_fasta,
    )


@dataclass
class LibraryDesignPaths:
    MERFISHAnalysisPath: str = "."

    @property
    def basePath(self) -> str:
        return str(Path(self.MERFISHAnalysisPath) / "MERFISH_Examples2") + "/"

    @property
    def rawTranscriptomeFasta(self) -> str:
        return str(Path(self.basePath) / "transcripts.fasta")

    @property
    def fpkmPath(self) -> str:
        return str(Path(self.basePath) / "isoforms.fpkm_tracking")

    @property
    def ncRNAPath(self) -> str:
        return str(Path(self.basePath) / "Homo_sapiens.GRCh38.ncrna.fa")

    @property
    def readoutPath(self) -> str:
        return str(Path(self.basePath) / "readouts.fasta")

    @property
    def codebookPath(self) -> str:
        return str(Path(self.basePath) / "codebook.csv")

    @property
    def analysisSavePath(self) -> str:
        return set_figure_save_path(Path(self.basePath) / "libraryDesign", make_dir=True)

    def derived(self) -> dict[str, str]:
        base = Path(self.analysisSavePath)
        return {
            "rRNAtRNAPath": str(base / "rRNAtRNA.fa"),
            "transcriptomePath": str(base / "transcriptomeObj"),
            "specificityTablePath": str(base / "specificityTable"),
            "isoSpecificityTablePath": str(base / "isoSpecificityTables"),
            "trDesignerPath": str(base / "trDesigner"),
            "trRegionsPath": str(base / "tr_GC_43_63_Tm_66_76_Len_30_30_IsoSpec_0.75_1_Spec_0.75_1"),
        }


def _require_probe_classes() -> tuple[Any, Any, Any, Any, Any]:
    transcriptome_cls = Transcriptome or import_symbol(["probe_construction.Transcriptome.Transcriptome"])
    ot_table_cls = OTTable or import_symbol(["probe_construction.OTTable.OTTable"])
    tr_designer_cls = TRDesigner or import_symbol(["probe_construction.TRDesigner.TRDesigner"])
    target_regions_cls = TargetRegions or import_symbol(["probe_construction.TargetRegions.TargetRegions"])
    primer_designer_cls = PrimerDesigner or import_symbol(["probe_construction.PrimerDesigner.PrimerDesigner"])
    return transcriptome_cls, ot_table_cls, tr_designer_cls, target_regions_cls, primer_designer_cls


def load_codebook_csv(path: str | Path) -> list[dict[str, str]]:
    if LoadCodebook is not None:
        loaded = LoadCodebook(str(path))
        if isinstance(loaded, list):
            return [dict(x) if isinstance(x, dict) else vars(x) for x in loaded]
    with Path(path).expanduser().open("r", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        rows = list(csv.DictReader(handle, dialect=dialect))
    normalized: list[dict[str, str]] = []
    for row in rows:
        local = {str(k).strip(): "" if v is None else str(v).strip() for k, v in row.items()}
        lowered = {k.lower(): v for k, v in local.items()}
        normalized.append(
            {
                "id": lowered.get("id", lowered.get("geneid", lowered.get("gene_id", ""))),
                "name": lowered.get("name", lowered.get("genename", lowered.get("gene_name", ""))),
                "barcode": lowered.get("barcode", lowered.get("codeword", lowered.get("code", ""))),
            }
        )
    return normalized


def extract_biotype(header: str) -> str:
    match = re.search(r"gene_biotype:(\S+)", header)
    return match.group(1) if match else ""


def load_or_create_rrna_trna(nc_rna_path: str, rrna_trna_path: str) -> list[FastaRecord]:
    out_path = Path(rrna_trna_path)
    if out_path.exists():
        print(f"Found and loading: {out_path}")
        timer = tic()
        records = read_fasta(out_path)
        print(f".... completed in {toc(timer)}s")
        print(f"Loaded {len(records)} sequences")
        return records

    print(f"Loading: {nc_rna_path}")
    timer = tic()
    nc_rnas = read_fasta(nc_rna_path)
    print(f"... completed in {toc(timer)} s")
    print(f"Found {len(nc_rnas)} sequences")
    biotypes_to_keep = {"rRNA", "tRNA", "Mt_rRNA", "Mt_tRNA"}
    page_break()
    print("Keeping the following types:")
    for biotype in sorted(biotypes_to_keep):
        print(f"   {biotype}")
    rrna_trna = [rec for rec in nc_rnas if extract_biotype(rec.Header) in biotypes_to_keep]
    print(f"Keeping {len(rrna_trna)} ncRNAs")
    if out_path.exists():
        out_path.unlink()
    write_fasta(out_path, rrna_trna)
    print(f"Wrote: {out_path}")
    return rrna_trna


def load_or_create_transcriptome(raw_fasta: str, fpkm_path: str, transcriptome_path: str, transcriptome_cls: Any) -> Any:
    path = Path(transcriptome_path)
    if not path.exists():
        transcriptome = transcriptome_cls(raw_fasta, abundPath=fpkm_path, verbose=True)
        transcriptome.Save(path)
        return transcriptome
    return transcriptome_cls.Load(path)


def build_isoform_specificity_tables(transcriptome: Any, save_path: str, ot_table_cls: Any) -> list[Any]:
    path = Path(save_path)
    if path.exists():
        loaded = ot_table_cls.Load(path, verbose=True)
        return loaded if isinstance(loaded, list) else [loaded]

    names = transcriptome.GetNames()
    ids_by_name = transcriptome.GetIDsByName(names)
    timer = tic()
    page_break()
    print("Starting construction of isoform specificity tables")
    print(f"Started on {time.ctime()}")
    tables: list[Any] = []
    for i, (name, gene_ids) in enumerate(zip(names, ids_by_name), start=1):
        local_transcriptome = transcriptome.Slice(geneID=gene_ids)
        table = ot_table_cls(local_transcriptome, 17, verbose=False, transferAbund=True)
        table.name = name
        tables.append(table)
        if i % 500 == 0:
            print(f"... completed {i} of {len(names)} genes")
    print(f".... completed in {toc(timer)} s")
    if hasattr(ot_table_cls, "Save"):
        # Object arrays are stored by each translated OTTable in its own directory.
        ensure_dir(path)
        for i, table in enumerate(tables, start=1):
            table.Save(path / f"table_{i:05d}")
    return tables


def build_total_specificity_table(iso_tables: list[Any], save_path: str, ot_table_cls: Any) -> Any:
    path = Path(save_path)
    if path.exists():
        return ot_table_cls.Load(path, verbose=True)
    first = iso_tables[0]
    first.verbose = True
    specificity_table = ot_table_cls.sum(iso_tables) if hasattr(ot_table_cls, "sum") else sum(iso_tables)
    specificity_table.name = "Transcriptome Specificity"
    first.verbose = False
    specificity_table.verbose = True
    specificity_table.Save(path)
    return specificity_table


def load_or_create_tr_designer(
    tr_designer_path: str,
    transcriptome: Any,
    rrna_table: Any,
    specificity_table: Any,
    iso_specificity_tables: list[Any],
    tr_designer_cls: Any,
) -> Any:
    path = Path(tr_designer_path)
    if not path.exists():
        designer = tr_designer_cls(
            transcriptome=transcriptome,
            OTTables=[rrna_table],
            OTTableNames=["rRNA"],
            specificityTable=specificity_table,
            isoSpecificityTables=iso_specificity_tables,
            parallel=None,
        )
        designer.Save(path)
        return designer
    return tr_designer_cls.Load(path)


def load_or_design_target_regions(tr_regions_path: str, tr_designer: Any, target_regions_cls: Any) -> list[Any]:
    path = Path(tr_regions_path)
    if not path.exists():
        target_regions = tr_designer.DesignTargetRegions(
            regionLength=30,
            GC=[0.43, 0.63],
            Tm=[66, 76],
            isoSpecificity=[0.75, 1],
            specificity=[0.75, 1],
            OTTables=["rRNA", [0, 0]],
        )
        print(f"... completed: {time.ctime()}")
        ensure_dir(path)
        for i, region in enumerate(target_regions, start=1):
            region.Save(path / f"region_{i:05d}")
        return target_regions
    print(f"Found: {path}")
    loaded = target_regions_cls.Load(path)
    return loaded if isinstance(loaded, list) else [loaded]


def filter_target_regions(target_regions: list[Any], final_ids: list[str]) -> list[Any]:
    final_id_set = set(final_ids)
    return [r for r in target_regions if getattr(r, "id", "") in final_id_set]


def build_possible_oligos(
    analysis_save_path: str,
    library_name: str,
    final_ids: list[str],
    final_genes: list[str],
    barcodes: np.ndarray,
    readouts: list[FastaRecord],
    final_target_regions: list[Any],
    rrna_table: Any,
    num_probes_per_gene: int,
) -> list[FastaRecord]:
    oligos_path = Path(analysis_save_path) / f"{library_name}_possible_oligos.fasta"
    if oligos_path.exists():
        raise FileExistsError("Found existing possible oligos file")

    oligos: list[FastaRecord] = []
    for i, local_gene_name in enumerate(final_genes):
        page_break()
        print(f"Designing probes for {library_name}: {local_gene_name}")
        possible_readouts = [readouts[j] for j in np.where(barcodes[i, :])[0]]
        matched_regions = [r for r in final_target_regions if getattr(r, "geneName", "") == local_gene_name]
        if not matched_regions:
            continue
        t_region = matched_regions[0]
        seqs: list[str] = []
        headers: list[str] = []
        for p in range(int(getattr(t_region, "numRegions", 0))):
            local_readouts = random.sample(possible_readouts, min(3, len(possible_readouts)))
            if len(local_readouts) < 3:
                continue
            start_pos = int(t_region.startPos[p])
            region_seq = t_region.sequence[p]
            gc = float(t_region.GC[p]) if len(t_region.GC) else float("nan")
            tm = float(t_region.Tm[p]) if len(t_region.Tm) else float("nan")
            specificity = float(t_region.specificity[p]) if len(t_region.specificity) else float("nan")
            core_name = f"{t_region.geneName}__{t_region.id}__{start_pos}__{len(region_seq)}__{gc}__{tm}__{specificity}"
            if random.random() > 0.5:
                headers.append(f"{library_name} {local_readouts[0].Header} {core_name} {local_readouts[1].Header} {local_readouts[2].Header}")
                seqs.append(
                    "A "
                    + reverse_complement(local_readouts[0].Sequence)
                    + " "
                    + reverse_complement(region_seq)
                    + " A "
                    + reverse_complement(local_readouts[1].Sequence)
                    + " "
                    + reverse_complement(local_readouts[2].Sequence)
                )
            else:
                headers.append(f"{library_name} {local_readouts[0].Header} {local_readouts[1].Header} {core_name} {local_readouts[2].Header}")
                seqs.append(
                    "A "
                    + reverse_complement(local_readouts[0].Sequence)
                    + " "
                    + reverse_complement(local_readouts[1].Sequence)
                    + " "
                    + reverse_complement(region_seq)
                    + " A "
                    + reverse_complement(local_readouts[2].Sequence)
                )
        print(f"... constructed {len(seqs)} possible probes")
        penalties = []
        for seq in seqs:
            penalty, _ = rrna_table.CalculatePenalty(reverse_complement(strip_spaces(seq)))
            penalties.append(np.nansum(penalty) > 0)
        keep_indices = [idx for idx, has_penalty in enumerate(penalties) if not has_penalty]
        remove_indices = [idx for idx in range(len(seqs)) if idx not in keep_indices]
        print(f"... removing {len(remove_indices)} probes")
        for idx in remove_indices:
            print(f"...     {headers[idx]}")
        random.shuffle(keep_indices)
        keep_indices = keep_indices[: min(len(keep_indices), num_probes_per_gene)]
        print(f"... keeping {len(keep_indices)} probes")
        if len(keep_indices) < num_probes_per_gene:
            print(f"Warning: not enough probes for {i + 1}: {local_gene_name}")
        for idx in keep_indices:
            oligos.append(FastaRecord(headers[idx], seqs[idx]))

    page_break()
    print(f"Writing: {oligos_path}")
    timer = tic()
    write_fasta(oligos_path, oligos)
    print(f"... completed in {toc(timer)}")
    return oligos


def design_primers(analysis_save_path: str, library_name: str, oligos: list[FastaRecord], ot_table_cls: Any, primer_designer_cls: Any) -> list[FastaRecord]:
    primers_path = Path(analysis_save_path) / f"{library_name}_possible_primers.fasta"
    if primers_path.exists():
        raise FileExistsError("Found existing primers")
    page_break()
    print(f"Designing primers for {library_name}")
    seq_rcomplement = [reverse_complement(strip_spaces(x.Sequence)) for x in oligos]
    all_seqs = [strip_spaces(x.Sequence) for x in oligos] + seq_rcomplement
    encoding_probe_table = ot_table_cls([{"Sequence": s} for s in all_seqs], 15, verbose=True, parallel=None)
    primer_designer = primer_designer_cls(
        numPrimersToGenerate=1000,
        primerLength=20,
        OTTables=encoding_probe_table,
        OTTableNames=["encoding"],
        parallel=None,
    )
    primer_designer.CutPrimers(Tm=[70, 72], GC=[0.5, 0.65], OTTables=["encoding", [0, 0]])
    primer_designer.RemoveForbiddenSeqs()
    primer_designer.RemoveSelfCompPrimers(homologyMax=6)
    primer_designer.RemoveHomologousPrimers(homologyMax=8)
    primer_designer.WriteFasta(primers_path)
    return read_fasta(primers_path)


def add_primers(analysis_save_path: str, library_name: str, oligos: list[FastaRecord], primers: list[FastaRecord]) -> list[FastaRecord]:
    used_primers = primers[:2]
    if len(used_primers) < 2:
        raise ValueError("At least two valid primers are required")
    page_break()
    print("Adding primers")
    final_primers_path = Path(analysis_save_path) / f"{library_name}_primers.fasta"
    if final_primers_path.exists():
        raise FileExistsError("Found existing final primers path")
    write_fasta(final_primers_path, used_primers)
    print(f"Wrote: {final_primers_path}")
    name1 = used_primers[0].Header.split()[0]
    name2 = used_primers[1].Header.split()[0]
    final_oligos: list[FastaRecord] = []
    for oligo in oligos:
        string_parts = oligo.Header.split()
        header = " ".join([string_parts[0], name1] + string_parts[1:] + [name2])
        sequence = used_primers[0].Sequence + " " + oligo.Sequence + " " + reverse_complement(used_primers[1].Sequence)
        final_oligos.append(FastaRecord(header, sequence))
    return final_oligos


def final_crosscheck_and_write(analysis_save_path: str, library_name: str, final_oligos: list[FastaRecord], rrna_table: Any) -> list[FastaRecord]:
    page_break()
    print("Running final cross checks and building final fasta file")
    oligos_path = Path(analysis_save_path) / f"{library_name}_oligos.fasta"
    if oligos_path.exists():
        print("Warning: found existing oligos")
        return read_fasta(oligos_path)
    timer = tic()
    print("Searching oligos for homology")
    has_rrna_penalty = []
    for oligo in final_oligos:
        penalty, _ = rrna_table.CalculatePenalty(reverse_complement(strip_spaces(oligo.Sequence)))
        has_rrna_penalty.append(np.nansum(penalty) > 0)
    print(f"... completed in {toc(timer)} s")
    kept = [oligo for oligo, has_penalty in zip(final_oligos, has_rrna_penalty) if not has_penalty]
    removed = [oligo for oligo, has_penalty in zip(final_oligos, has_rrna_penalty) if has_penalty]
    print(f"... found {len(removed)} oligos to remove")
    for oligo in removed:
        print(f"...     {oligo.Header}")
    write_fasta(oligos_path, kept)
    print(f"Wrote: {oligos_path}")
    return kept


def run_library_design_example(paths: LibraryDesignPaths | None = None, num_workers: int = 8) -> dict[str, Any]:
    paths = paths or LibraryDesignPaths()
    derived = paths.derived()
    transcriptome_cls, ot_table_cls, tr_designer_cls, target_regions_cls, primer_designer_cls = _require_probe_classes()

    page_break()
    rrna_trna = load_or_create_rrna_trna(paths.ncRNAPath, derived["rRNAtRNAPath"])
    transcriptome = load_or_create_transcriptome(paths.rawTranscriptomeFasta, paths.fpkmPath, derived["transcriptomePath"], transcriptome_cls)
    iso_specificity_tables = build_isoform_specificity_tables(transcriptome, derived["isoSpecificityTablePath"], ot_table_cls)
    specificity_table = build_total_specificity_table(iso_specificity_tables, derived["specificityTablePath"], ot_table_cls)

    print(f"Configured parallel pool equivalent with {num_workers} workers")
    rrna_table = ot_table_cls(rrna_trna, 15, verbose=True, parallel=None)

    page_break()
    print("Slicing transcriptome based on expression level: >= 1e-2 FPKM")
    ids = list(transcriptome.ids)
    abundance = np.asarray(transcriptome.GetAbundanceByID(ids), dtype=float)
    good_ids = [id_value for id_value, value in zip(ids, abundance) if value >= 1e-2]
    sliced_transcriptome = transcriptome.Slice(geneID=good_ids)

    tr_designer = load_or_create_tr_designer(
        derived["trDesignerPath"],
        sliced_transcriptome,
        rrna_table,
        specificity_table,
        iso_specificity_tables,
        tr_designer_cls,
    )
    target_regions = load_or_design_target_regions(derived["trRegionsPath"], tr_designer, target_regions_cls)

    page_break()
    timer = tic()
    print(f"Loading: {paths.readoutPath}")
    readouts = read_fasta(paths.readoutPath)
    print(f"Found {len(readouts)} oligos in {toc(timer)} s")
    codebook = load_codebook_csv(paths.codebookPath)
    final_ids = [row["id"] for row in codebook]
    final_genes = [row["name"] for row in codebook]
    barcodes = np.asarray([[char == "1" for char in row["barcode"]] for row in codebook], dtype=bool)
    final_target_regions = filter_target_regions(target_regions, final_ids)

    num_probes_per_gene = 92
    library_name = "L1E1"
    page_break()
    print(f"Designing oligos for {library_name}")
    print(f"... {num_probes_per_gene} probes per gene")
    used_readout_path = Path(paths.analysisSavePath) / f"{library_name}_used_readouts.fasta"
    if used_readout_path.exists():
        print(f"Warning: found {used_readout_path}")
        used_readout_path.unlink()
    write_fasta(used_readout_path, readouts)
    page_break()
    print(f"Wrote {len(readouts)} readouts to {used_readout_path}")

    oligos = build_possible_oligos(
        paths.analysisSavePath,
        library_name,
        final_ids,
        final_genes,
        barcodes,
        readouts,
        final_target_regions,
        rrna_table,
        num_probes_per_gene,
    )
    primers = design_primers(paths.analysisSavePath, library_name, oligos, ot_table_cls, primer_designer_cls)
    final_oligos = add_primers(paths.analysisSavePath, library_name, oligos, primers)
    kept_oligos = final_crosscheck_and_write(paths.analysisSavePath, library_name, final_oligos, rrna_table)

    source_file = Path(__file__)
    if source_file.exists():
        shutil.copy2(source_file, Path(paths.analysisSavePath) / source_file.name)
        page_break()
        print(f"Copied analysis script to {paths.analysisSavePath}{source_file.name}")
    return {
        "paths": asdict(paths),
        "derived": derived,
        "numReadouts": len(readouts),
        "numCodebookRows": len(codebook),
        "numTargetRegions": len(target_regions),
        "numOligos": len(oligos),
        "numFinalOligos": len(kept_oligos),
    }


def main() -> dict[str, Any]:
    return run_library_design_example()


if __name__ == "__main__":
    main()
