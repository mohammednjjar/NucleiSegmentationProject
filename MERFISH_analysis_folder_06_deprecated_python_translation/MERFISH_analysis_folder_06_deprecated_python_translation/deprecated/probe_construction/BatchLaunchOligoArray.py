"""Python translation of `deprecated/probe_construction/BatchLaunchOligoArray.m`."""
from __future__ import annotations

from pathlib import Path
import re
import time
from typing import Sequence, Mapping, Any

try:
    from ..helpers import read_fasta, sanitize_locus_name, write_fasta_records, run_command
except ImportError:
    from deprecated.helpers import read_fasta, sanitize_locus_name, write_fasta_records, run_command


def _extract_save_path(command: str) -> str:
    match = re.search(r'\s-o\s+(.+?)genename_oligos\.txt', command)
    if match:
        return match.group(1)
    return './'


def BatchLaunchOligoArray(oligoArrayCommandtemp: str, geneFasta: str | Sequence[Mapping[str, Any]],
                          savePath: str = '', headerSeps: str = '', maxFragment: int = 1000,
                          batchsize: int = 4, maxTime: float = 30, maxBlastTime: float = 15,
                          maxCPU: float = 75, refreshTime: float = 20, runExternal: bool = True,
                          verbose: bool = True, runSbatch: bool = False,
                          remotePath: str = '/n/home05/boettiger/', execute: bool | None = None):
    if execute is None:
        execute = bool(runExternal or runSbatch)
    records = read_fasta(geneFasta) if isinstance(geneFasta, (str, Path)) else [dict(x) for x in geneFasta]
    if not savePath:
        savePath = _extract_save_path(oligoArrayCommandtemp)
    Path(savePath).mkdir(parents=True, exist_ok=True)
    launched = []
    for g, rec in enumerate(records, start=1):
        header = str(rec.get('Header', ''))
        sequence = str(rec.get('Sequence', ''))
        locus_name = sanitize_locus_name(header)
        if headerSeps:
            parts = header.split(headerSeps)
            if len(parts) >= 3:
                locus_name = parts[1]
        fragments = []
        for start in range(0, len(sequence), int(maxFragment)):
            fragment = sequence[start:start + int(maxFragment)]
            if fragment:
                fragments.append({'Header': f'{locus_name}_pt{len(fragments) + 1}', 'Sequence': fragment})
        fasta_out = Path(savePath) / f'{locus_name}.fasta'
        write_fasta_records(fasta_out, fragments, append=False)
        command = oligoArrayCommandtemp.replace('genename', locus_name)
        if verbose:
            print('-----------------------------------------------------------------')
            print(f'Running file {g} of {len(records)}:')
            print(f'     {command}')
            print('-----------------------------------------------------------------')
        if runSbatch:
            bat_file = Path(remotePath) / f'Run{locus_name}.sh'
            bat_file.parent.mkdir(parents=True, exist_ok=True)
            bat_file.write_text('\n'.join([
                '#!/bin/bash', '#SBATCH -n 1', f'#SBATCH -t {int(2 * maxTime)}',
                '#SBATCH -p serial_requeue', '#SBATCH --mem=3000', command
            ]), encoding='utf-8')
            result = run_command(f'sbatch {bat_file}', execute=execute)
        else:
            result = run_command(command, execute=execute)
        result['fasta'] = str(fasta_out)
        launched.append(result)
        if execute and batchsize > 0 and (g % batchsize) == 0:
            time.sleep(0.05)
    return launched


def batch_launch_oligo_array(*args, **kwargs):
    return BatchLaunchOligoArray(*args, **kwargs)
