"""Python translation of `deprecated/probe_construction/GenProbe.m`."""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any, Mapping, Sequence
import numpy as np
import random

try:
    from ..helpers import table_to_records, reverse_complement, write_fasta_records
except ImportError:
    from deprecated.helpers import table_to_records, reverse_complement, write_fasta_records


def _as_bool_list(value, length):
    if isinstance(value, (list, tuple, np.ndarray)):
        vals = [bool(x) for x in value]
    else:
        vals = [bool(value)]
    if len(vals) < length:
        vals.extend([vals[-1]] * (length - len(vals)))
    return vals


def GenProbe(libCodebook, pickedGenes, secondaries, indexPrimers, numOligos: int, numGenesAndCntrls: int,
             saveFolder: str, bitsPerProbe: int = 2, subLibNum: int = 1, primerNum: int = 1,
             universal: str = '', distanceSort: bool = True, shuffleCodewords=1, shufflePriSec=1,
             shufflePri=1, changeSecondaries=0, saveIsoName: bool = True, seed: int | None = None):
    rng = random.Random(seed)
    codebook = np.asarray(libCodebook, dtype=np.uint8)
    genes = table_to_records(pickedGenes)
    seconds = table_to_records(secondaries)
    primers = table_to_records(indexPrimers)
    if len(primers) < primerNum + 2:
        raise ValueError('not enough index primers')
    num_words, num_letters = codebook.shape
    num_shuffles = max(len(np.atleast_1d(shuffleCodewords)), len(np.atleast_1d(shufflePriSec)), len(np.atleast_1d(shufflePri)), len(np.atleast_1d(changeSecondaries)))
    scw = _as_bool_list(shuffleCodewords, num_shuffles)
    sps = _as_bool_list(shufflePriSec, num_shuffles)
    sp = _as_bool_list(shufflePri, num_shuffles)
    cs = _as_bool_list(changeSecondaries, num_shuffles)
    Path(saveFolder).mkdir(parents=True, exist_ok=True)
    all_oligo_seq: list[list[str]] = []
    all_oligo_name: list[list[str]] = []
    first_codebook = codebook.copy()
    pri_sec_shuffle: dict[tuple[int, int], list[int]] = {}
    pri_draw_shuffle: dict[tuple[int, int], list[int]] = {}
    common_flag = 'universal ' if universal else ''
    for r in range(num_shuffles):
        primerNum += 1
        fwd = primers[primerNum - 1]['Sequence']
        fwd_name = primers[primerNum - 1]['Header']
        primerNum += 1
        rev = reverse_complement(primers[primerNum - 1]['Sequence'])
        rev_name = primers[primerNum - 1]['Header']
        current_seconds = seconds[r * num_letters:(r + 1) * num_letters] if cs[r] else seconds[:num_letters]
        if scw[r]:
            order = list(range(num_words))
            rng.shuffle(order)
            local_codebook = codebook[order, :]
            if distanceSort and num_words > 1:
                ref = local_codebook[0]
                dists = np.mean(local_codebook != ref, axis=1)
                local_codebook = local_codebook[np.argsort(dists), :]
            first_codebook = local_codebook.copy()
        else:
            local_codebook = first_codebook.copy()
        codebook_records = []
        for gene, word in zip(genes, local_codebook):
            common = str(gene.get('CommonName', ''))
            iso = str(gene.get('IsoformName', ''))
            name = f'{common}    {iso}' if saveIsoName else common
            secondary_names = ' '.join(str(s['Header']) for s in current_seconds[:num_letters])
            codebook_records.append({'Header': ' '.join(str(int(x)) for x in word), 'Sequence': f'{name}    {secondary_names}'})
        write_fasta_records(Path(saveFolder) / f'E{subLibNum}_codebook.fasta', codebook_records, append=False)
        oligo_file = Path(saveFolder) / f'E{subLibNum}_oligos.fasta'
        if oligo_file.exists():
            oligo_file.unlink()
        local_seq: list[str] = []
        local_name: list[str] = []
        j = 1
        for n, gene in enumerate(genes[:num_words]):
            common = str(gene.get('CommonName', ''))
            if 'blank' in common:
                continue
            gene_on_bits = [idx for idx, val in enumerate(local_codebook[n]) if val]
            sec_combos = list(combinations(gene_on_bits, int(bitsPerProbe)))
            if not sec_combos:
                continue
            seqs = list(gene.get('Sequence', []) or [])
            fiveprime = list(gene.get('FivePrimeEnd', [0] * len(seqs)) or [])
            num_oligos_used = min(len(seqs), int(numOligos))
            combo_count = len(sec_combos)
            num_oligos_used = (num_oligos_used // combo_count) * combo_count
            if num_oligos_used == 0:
                continue
            if sps[r]:
                pri_order = list(range(num_oligos_used))
                rng.shuffle(pri_order)
                pri_sec_shuffle[(r, n)] = pri_order
            else:
                pri_order = pri_sec_shuffle.get((0, n), list(range(num_oligos_used)))
            if sp[r]:
                draw = list(range(len(seqs)))
                rng.shuffle(draw)
                pri_draw_shuffle[(r, n)] = draw
            else:
                draw = pri_draw_shuffle.get((0, n), list(range(len(seqs))))
            selected_seqs = [seqs[i] for i in draw[:num_oligos_used]]
            selected_pos = [fiveprime[i] if i < len(fiveprime) else 0 for i in draw[:num_oligos_used]]
            per_combo = num_oligos_used // combo_count
            for k, sec_idx_tuple in enumerate(sec_combos):
                for i in range(per_combo):
                    local_index = pri_order[k * per_combo + i]
                    sec_parts = []
                    sec_names = []
                    for sec_idx in sec_idx_tuple:
                        sec_parts.append(reverse_complement(current_seconds[sec_idx]['Sequence']))
                        sec_names.append(str(current_seconds[sec_idx]['Header']))
                    while len(sec_parts) < 4:
                        sec_parts.append('')
                        sec_names.append('')
                    oligo = ' '.join(x for x in [fwd, universal, sec_parts[2], sec_parts[0], reverse_complement(str(selected_seqs[local_index])), sec_parts[1], sec_parts[3], rev] if x != '')
                    iso = str(gene.get('IsoformName', ''))
                    pos = selected_pos[local_index]
                    name = ' '.join(x for x in [f'E{subLibNum}', fwd_name, common_flag.strip(), sec_names[2], sec_names[0], f'{common}__{iso}__{pos}', sec_names[1], sec_names[3], rev_name, f'probe{j:03d}'] if x != '')
                    local_seq.append(oligo)
                    local_name.append(name)
                    write_fasta_records(oligo_file, [{'Header': name, 'Sequence': oligo}], append=True)
                    j += 1
        all_oligo_seq.append(local_seq)
        all_oligo_name.append(local_name)
        subLibNum += 1
    return all_oligo_seq, all_oligo_name, primerNum, {'bitsPerProbe': bitsPerProbe, 'subLibNum': subLibNum, 'primerNum': primerNum}


def gen_probe(*args, **kwargs):
    return GenProbe(*args, **kwargs)
