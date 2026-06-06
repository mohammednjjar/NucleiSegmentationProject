"""Python translation of `deprecated/probe_construction/AssembleProbes.m`."""
from __future__ import annotations

from pathlib import Path
import math

try:
    from ..helpers import read_fasta, gen_secded
    from .AddCntrlSeqs import AddCntrlSeqs
    from .GenProbe import GenProbe
except ImportError:
    from deprecated.helpers import read_fasta, gen_secded
    from deprecated.probe_construction.AddCntrlSeqs import AddCntrlSeqs
    from deprecated.probe_construction.GenProbe import GenProbe


def AssembleProbes(ProbeData, primerFasta: str = '', readoutFasta: str = '', cntrlFasta: str = '',
                   dataFolder: str = '', libSaveFolder: str = '', numCntrls: int = 0,
                   numBlanks: int = 0, numOligos: int = 200, subLibNum: int = 1,
                   universal: bool = True, numTotalBits: int = 16, numDataBits: int = 11,
                   onBits: int = 4, bitsPerProbe: int = 2, seed: int | None = None):
    if dataFolder:
        base = Path(dataFolder)
        primerFasta = str(base / 'PrimerSeqs.fasta')
        readoutFasta = str(base / 'ReadoutSeqs.fasta')
        cntrlFasta = str(base / 'BLASTedRandomSeqs.fasta')
    if not readoutFasta:
        raise ValueError('dataFolder or readoutFasta is required')
    if not libSaveFolder:
        libSaveFolder = str(Path(readoutFasta).parent / 'AssembledProbes')
    Path(libSaveFolder).mkdir(parents=True, exist_ok=True)
    index_primers = read_fasta(primerFasta)
    universal_seq = index_primers[0]['Sequence'] if universal and index_primers else ''
    cntrl = read_fasta(cntrlFasta) if cntrlFasta else []
    readout_seqs = read_fasta(readoutFasta)
    numOligos = math.floor(numOligos / math.comb(onBits, bitsPerProbe)) * math.comb(onBits, bitsPerProbe)
    codebook_table = gen_secded(numTotalBits, numDataBits, onBits)
    num_genes = len(ProbeData.get('GeneName', ProbeData.get('CommonName', []))) if isinstance(ProbeData, dict) else len(ProbeData)
    num_words = num_genes + numCntrls + numBlanks
    if codebook_table.shape[0] < num_words:
        raise ValueError('codebook is not large enough to encode all genes and controls')
    if codebook_table.shape[0] > num_words:
        numBlanks = int(codebook_table.shape[0] - numCntrls - num_genes)
    picked = AddCntrlSeqs(ProbeData, cntrl, numOligos, numCntrls, numBlanks, seed=seed)
    oligo_seq, oligo_name, _, _ = GenProbe(codebook_table, picked, readout_seqs, index_primers,
                                           numOligos, num_genes + numCntrls, libSaveFolder,
                                           bitsPerProbe=bitsPerProbe, subLibNum=subLibNum,
                                           universal=universal_seq, seed=seed)
    return oligo_seq, oligo_name, codebook_table


def assemble_probes(*args, **kwargs):
    return AssembleProbes(*args, **kwargs)
