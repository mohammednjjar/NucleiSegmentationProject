"""Python translation of `deprecated/probe_construction/AddCntrlSeqs.m`."""
from __future__ import annotations

import random
from typing import Any, Mapping, Sequence

try:
    from ..helpers import table_to_records, records_to_table
except ImportError:
    from deprecated.helpers import table_to_records, records_to_table


def _insert(records, idx, rec):
    records[idx:idx] = [rec]


def AddCntrlSeqs(pickedGenes: Mapping[str, Any] | Sequence[Mapping[str, Any]], cntrl: Sequence[Mapping[str, Any]],
                 numOligos: int, numCntrls: int, numBlanks: int, distribute: bool = True, seed: int | None = None):
    rng = random.Random(seed)
    records = table_to_records(pickedGenes)
    rand_seqs = [str(c.get('Sequence', '')) for c in cntrl]
    rng.shuffle(rand_seqs)
    j = 0
    if distribute:
        if numCntrls > 0:
            insert_positions = [round((len(records) + numCntrls) / numCntrls * (n + 1)) - 1 for n in range(numCntrls)]
            for n, i in enumerate(insert_positions, start=1):
                seqs = rand_seqs[j:j + 2 * numOligos]
                rec = {'CommonName': f'notarget{n:03d}', 'Sequence': seqs,
                       'FivePrimeEnd': [40 * (k + 1) for k in range(j, j + len(seqs))],
                       'IsoformName': '', 'GeneName': '', 'Nprobes': 2 * numOligos, 'FPKM': 0,
                       'MeltingTemp': [], 'FreeEnergy': [], 'Enthalpy': [], 'Entropy': [], 'GeneNumber': [], 'Locus': ''}
                _insert(records, min(max(i, 0), len(records)), rec)
                j += max(2 * numOligos - 1, 1)
        if numBlanks > 0:
            insert_positions = [round((len(records) + numBlanks) / numBlanks * (n + 1)) - 1 for n in range(numBlanks)]
            for n, i in enumerate(insert_positions, start=1):
                rec = {'CommonName': f'blank{n:03d}', 'Sequence': [], 'FivePrimeEnd': [j + 1],
                       'IsoformName': '', 'GeneName': '', 'Nprobes': 0, 'FPKM': 0,
                       'MeltingTemp': [], 'FreeEnergy': [], 'Enthalpy': [], 'Entropy': [], 'GeneNumber': [], 'Locus': ''}
                _insert(records, min(max(i, 0), len(records)), rec)
                j += max(2 * numOligos - 1, 1)
    else:
        for n in range(1, numCntrls + 1):
            seqs = rand_seqs[j:j + 2 * numOligos]
            records.append({'CommonName': f'notarget{n:03d}', 'Sequence': seqs,
                            'FivePrimeEnd': [40 * (k + 1) for k in range(j, j + len(seqs))],
                            'IsoformName': '', 'GeneName': '', 'Nprobes': 2 * numOligos, 'FPKM': 0,
                            'MeltingTemp': [], 'FreeEnergy': [], 'Enthalpy': [], 'Entropy': [], 'GeneNumber': [], 'Locus': ''})
            j += max(2 * numOligos - 1, 1)
        for n in range(1, numBlanks + 1):
            records.append({'CommonName': f'blank{n:03d}', 'Sequence': [], 'FivePrimeEnd': [j + 1],
                            'IsoformName': '', 'GeneName': '', 'Nprobes': 0, 'FPKM': 0,
                            'MeltingTemp': [], 'FreeEnergy': [], 'Enthalpy': [], 'Entropy': [], 'GeneNumber': [], 'Locus': ''})
            j += max(2 * numOligos - 1, 1)
    return records_to_table(records)


def add_cntrl_seqs(*args, **kwargs):
    return AddCntrlSeqs(*args, **kwargs)
