"""Python translation of `example_scripts/deprecated/ExampleLibraryConstruction_140genes_script.m`.

Purpose: demonstrate the older OligoArray-based construction of a MERFISH
probe library for 140 target genes.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

try:
    from ..script_utils import ensure_dir, import_symbol, read_fasta, tic, toc
except ImportError:
    from script_utils import ensure_dir, import_symbol, read_fasta, tic, toc


@dataclass
class DeprecatedLibraryConstructionPaths:
    MERFISHAnalysisPath: str = "."
    legacyBLASTPath: str = ""
    oligoArrayAuxPath: str = ""
    oligoArrayExe: str = ""

    @property
    def exampleDataPath(self) -> str:
        return str(Path(self.MERFISHAnalysisPath) / "MERFISH_Examples") + "/"

    @property
    def dataFolder(self) -> str:
        return str(Path(self.exampleDataPath) / "probe_construction_data") + "/"

    @property
    def geneFasta(self) -> str:
        return str(Path(self.dataFolder) / "TargetGeneSeqs.fasta")

    @property
    def blastLib(self) -> str:
        return str(Path(self.dataFolder) / "TargetGeneSeqs.fasta")

    @property
    def saveFolder(self) -> str:
        return str(Path(self.exampleDataPath) / "probe_construction_output") + "/"


def run_deprecated_library_construction(paths: DeprecatedLibraryConstructionPaths | None = None) -> dict[str, Any]:
    paths = paths or DeprecatedLibraryConstructionPaths()
    ensure_dir(paths.saveFolder)
    time_run = tic()

    build_blast_lib = import_symbol([
        "deprecated.probe_construction.BuildBLASTlib.BuildBLASTlib",
        "BuildBLASTlib.BuildBLASTlib",
    ])
    oligo_array_cmd = import_symbol([
        "deprecated.probe_construction.OligoArrayCmd.OligoArrayCmd",
        "OligoArrayCmd.OligoArrayCmd",
    ])
    batch_launch_oligo_array = import_symbol([
        "deprecated.probe_construction.BatchLaunchOligoArray.BatchLaunchOligoArray",
        "BatchLaunchOligoArray.BatchLaunchOligoArray",
    ])
    compile_oligo_array_output = import_symbol([
        "deprecated.probe_construction.CompileOligoArrayOutput.CompileOligoArrayOutput",
        "CompileOligoArrayOutput.CompileOligoArrayOutput",
    ])
    assemble_probes = import_symbol([
        "deprecated.probe_construction.AssembleProbes.AssembleProbes",
        "AssembleProbes.AssembleProbes",
    ])

    build_blast_lib(paths.blastLib, legacy=True, blastPath=str(Path(paths.legacyBLASTPath) / ""))
    oligo_array_command, save_path = oligo_array_cmd(
        savePath=paths.saveFolder,
        saveDir=paths.saveFolder,
        minTm=70,
        probeLength=30,
        crosshybeT=72,
        secstructT=76,
        blastLib=paths.blastLib,
        maskedSeq='"GGGG;CCCC;TTTTTT;AAAA"',
        numParallel=5,
        GbMem=1,
        oligoArrayExe=paths.oligoArrayExe,
    )
    gene_list = read_fasta(paths.geneFasta)
    set_path_command = f"PATH={paths.legacyBLASTPath};{paths.oligoArrayAuxPath};%PATH%; & "
    oligo_array_command = set_path_command + oligo_array_command

    batch_launch_oligo_array(
        oligo_array_command,
        gene_list,
        batchsize=20,
        maxTime=30,
        savePath=save_path,
        maxCPU=95,
        runExternal=True,
    )
    compile_oligo_array_output(save_path, blastLib=paths.blastLib, zeroOffTarget=True)

    probe_data_path = Path(paths.saveFolder) / "ProbeData.mat"
    try:
        from scipy.io import loadmat
        probe_payload = loadmat(str(probe_data_path), squeeze_me=True, struct_as_record=False)
        probe_data = probe_payload.get("ProbeData")
    except Exception as exc:
        raise RuntimeError(f"Could not load ProbeData from {probe_data_path}") from exc

    assemble_probes(
        probe_data,
        numCntrls=5,
        numBlanks=5,
        numOligos=200,
        subLibNum=1,
        numTotalBits=16,
        onBits=4,
        bitsPerProbe=2,
        universal=False,
        primerFasta=str(Path(paths.dataFolder) / "PrimerSeqs.fasta"),
        readoutFasta=str(Path(paths.dataFolder) / "ReadoutSeqs.fasta"),
        cntrlFasta=str(Path(paths.dataFolder) / "BLASTedRandomSeqs.fasta"),
        libSaveFolder=str(Path(paths.saveFolder) / "AssembledProbes") + "/",
    )
    total_minutes = toc(time_run) / 60.0
    print(f"Probe library assembly completed in {total_minutes:.5g} min.")
    return {"paths": asdict(paths), "oligoArrayCommand": oligo_array_command, "savePath": save_path, "totalMinutes": total_minutes}


def main() -> dict[str, Any]:
    return run_deprecated_library_construction()


if __name__ == "__main__":
    main()
