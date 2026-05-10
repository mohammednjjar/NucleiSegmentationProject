"""Python translation of `deprecated/probe_construction/OligoArrayCmd.m`."""
from __future__ import annotations

from datetime import date
from pathlib import Path
import math


def OligoArrayCmd(blastLib: str = r'D:\Data\Genomics\DmelGenome\legacyBLASTlib\Dmel_Genome.fasta',
                  savePath: str = '', saveDir: str = '', inputPath: str = '',
                  oligoArrayExe: str = r'C:\Users\Alistair\Documents\Research\Software\OligoArray\OligoArray2.jar',
                  probeLength: int = 42, probeLengthMax: int | None = None, minTm: float = 80,
                  crosshybeT: float = 72, secstructT: float = 72, maxFragment: float = 1e3,
                  GbMem: float = 3, maskedSeq: str = '"GGGG;CCCC"', minPercentGC: float = 10,
                  maxPercentGC: float = 90, maxTm: float = 100, minProbeDist: float | None = None,
                  fastaName: str = '', numParallel: int = 1, verbose: bool = True):
    if probeLengthMax is None:
        probeLengthMax = probeLength
    if not saveDir and not savePath:
        raise ValueError('Specify saveDir or savePath')
    if saveDir:
        saveFolder = f'{date.today():%y-%m-%d}_Probes_{int(probeLength)}mers_t{minTm:g}Cx{crosshybeT:g}Cs{secstructT:g}C/'
        savePath = str(Path(saveDir) / saveFolder)
        Path(savePath).mkdir(parents=True, exist_ok=True)
        if verbose:
            print(f'Creating Folder {savePath}')
    if not inputPath:
        inputPath = savePath
    if minProbeDist is None:
        minProbeDist = probeLengthMax + 1
    n_oligos = round(maxFragment / probeLength)
    mem = math.ceil(GbMem)
    command = (
        f'java -Xmx{mem}g -jar {oligoArrayExe}'
        f' -i {inputPath}genename.fasta'
        f' -d {blastLib}'
        f' -o {savePath}genename_oligos.txt'
        f' -r {savePath}genename_failed.txt'
        f' -R {savePath}genename_log.txt'
        f' -n {n_oligos}'
        f' -l {probeLength}'
        f' -L {probeLengthMax}'
        f' -D {int(maxFragment)}'
        f' -t {minTm:g}'
        f' -T {maxTm:g}'
        f' -s {secstructT:g}'
        f' -x {crosshybeT:g}'
        f' -p {minPercentGC:g}'
        f' -P {maxPercentGC:g}'
        f' -m {maskedSeq}'
        f' -g {minProbeDist:g}'
        f' -N {int(numParallel)}'
    )
    if fastaName:
        command = command.replace('genename', Path(fastaName).stem)
    return command, savePath


def oligo_array_cmd(*args, **kwargs):
    return OligoArrayCmd(*args, **kwargs)
