"""Python translation of `deprecated/probe_construction/CompileOligoArrayOutput.m`."""
from __future__ import annotations

from pathlib import Path
import csv
import math
from typing import Sequence

try:
    from ..helpers import read_fasta, save_pickle_and_json
    from ..misc.StringFind import StringFind
except ImportError:
    from deprecated.helpers import read_fasta, save_pickle_and_json
    from deprecated.misc.StringFind import StringFind


def _to_float(value: str):
    try:
        return float(value)
    except Exception:
        return float('nan')


def CompileOligoArrayOutput(oligo_folder: str, savePath: str = '', blastLib: str = '', zeroOffTarget: bool = True):
    folder = Path(oligo_folder)
    if not folder.exists():
        raise FileNotFoundError(folder)
    output_dir = Path(savePath) if savePath else folder
    output_dir.mkdir(parents=True, exist_ok=True)
    all_names = [rec['Header'] for rec in read_fasta(blastLib)] if blastLib else ['']
    oligofiles = sorted(folder.glob('*_oligos.txt'))
    probe_data = {
        'Nprobes': [], 'GeneName': [], 'CommonName': [], 'FivePrimeEnd': [], 'FreeEnergy': [],
        'Enthalpy': [], 'Entropy': [], 'MeltingTemp': [], 'Sequence': [], 'IsoformName': [],
        'FPKM': [], 'GeneNumber': [], 'Locus': []
    }
    for file in oligofiles:
        genename = file.name.replace('_oligos.txt', '')
        rows = []
        with file.open('r', encoding='utf-8', errors='ignore') as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) >= 10:
                    rows.append(row)
        kept = []
        if rows:
            for row in rows:
                off_target = False
                if zeroOffTarget:
                    target_mask, _ = StringFind(all_names, genename, boolean=True)
                    off_target_names = [n for n, is_target in zip(all_names, target_mask) if not is_target]
                    hits, _ = StringFind([row[7]], off_target_names, cellOutput=True)
                    off_target = any(len(h) > 0 for h in hits)
                if not off_target:
                    kept.append(row)
        probe_data['Nprobes'].append(len(kept))
        probe_data['GeneName'].append(genename)
        probe_data['CommonName'].append(genename)
        probe_data['IsoformName'].append('')
        probe_data['FPKM'].append(float('nan'))
        probe_data['GeneNumber'].append(None)
        probe_data['Locus'].append('')
        if kept:
            probe_data['FivePrimeEnd'].append([_to_float(r[1]) for r in kept])
            probe_data['FreeEnergy'].append([_to_float(r[3]) for r in kept])
            probe_data['Enthalpy'].append([_to_float(r[4]) for r in kept])
            probe_data['Entropy'].append([_to_float(r[5]) for r in kept])
            probe_data['MeltingTemp'].append([_to_float(r[6]) for r in kept])
            probe_data['Sequence'].append([r[8] for r in kept])
        else:
            probe_data['FivePrimeEnd'].append([float('nan')])
            probe_data['FreeEnergy'].append([float('nan')])
            probe_data['Enthalpy'].append([float('nan')])
            probe_data['Entropy'].append([float('nan')])
            probe_data['MeltingTemp'].append([float('nan')])
            probe_data['Sequence'].append([])
    pkl_path, json_path = save_pickle_and_json(probe_data, output_dir / 'ProbeData')
    try:
        from scipy.io import savemat
        savemat(output_dir / 'ProbeData.mat', {'ProbeData': probe_data})
    except Exception:
        (output_dir / 'ProbeData.mat.txt').write_text('Saved ProbeData.pkl and ProbeData.json as Python equivalents.\n', encoding='utf-8')
    return probe_data


def compile_oligo_array_output(*args, **kwargs):
    return CompileOligoArrayOutput(*args, **kwargs)
