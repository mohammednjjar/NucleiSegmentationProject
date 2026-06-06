"""Python translation of `deprecated/probe_construction/BuildBLASTlib.m`."""
from __future__ import annotations

from pathlib import Path
import shutil

try:
    from ..helpers import run_command
except ImportError:
    from deprecated.helpers import run_command


def BuildBLASTlib(fastaFile: str, blastPath: str = '', legacy: bool = False, execute: bool = True):
    fasta = Path(fastaFile)
    if not fasta.exists():
        raise FileNotFoundError(f'A valid FASTA path is required: {fasta}')
    executable = 'formatdb' if legacy else 'makeblastdb'
    exe_name = executable + ('.exe' if blastPath else '')
    exe_path = str(Path(blastPath) / exe_name) if blastPath else executable
    if execute and shutil.which(exe_path) is None and not Path(exe_path).exists():
        raise FileNotFoundError(f'Could not find {executable}; set blastPath or use execute=False')
    database_out = str(fasta.with_suffix(''))
    if legacy:
        command = f'{exe_path} -i {fasta} -o T -p F'
    else:
        command = f'{exe_path} -dbtype nucl -in {fasta} -parse_seqids -out {database_out}'
    result = run_command(command, execute=execute)
    result['database'] = database_out
    return result


def build_blast_lib(*args, **kwargs):
    return BuildBLASTlib(*args, **kwargs)
