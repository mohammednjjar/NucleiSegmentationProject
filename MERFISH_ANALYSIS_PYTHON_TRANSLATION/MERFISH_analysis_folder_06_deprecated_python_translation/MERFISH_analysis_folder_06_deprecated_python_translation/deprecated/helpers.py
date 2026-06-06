"""Shared Python utilities for the translated root `deprecated/` folder.

The functions in this file replace small MATLAB helpers used by the old MERFISH
probe-construction/report scripts: FASTA parsing/writing, structure-field access,
reverse-complement sequence handling, Hamming-neighbor generation, and plotting.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import json
import math
import os
import pickle
import re
import subprocess
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


Record = dict[str, Any]


def get_field(obj: Any, name: str, default: Any = None) -> Any:
    if isinstance(obj, Mapping):
        return obj.get(name, default)
    return getattr(obj, name, default)


def set_field(obj: Any, name: str, value: Any) -> None:
    if isinstance(obj, dict):
        obj[name] = value
    else:
        setattr(obj, name, value)


def ensure_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, list):
        return value
    return [value]


def table_to_records(table: Any) -> list[Record]:
    """Convert MATLAB-style struct arrays or dict-of-lists into list-of-dicts."""
    if table is None:
        return []
    if isinstance(table, list):
        return [dict(x) if isinstance(x, Mapping) else vars(x).copy() for x in table]
    if isinstance(table, tuple):
        return table_to_records(list(table))
    if isinstance(table, Mapping):
        values = {k: ensure_list(v) for k, v in table.items()}
        lengths = [len(v) for v in values.values()]
        if not lengths:
            return []
        n = max(lengths)
        records: list[Record] = []
        for i in range(n):
            row: Record = {}
            for k, v in values.items():
                if len(v) == n:
                    row[k] = v[i]
                elif len(v) == 1:
                    row[k] = v[0]
                else:
                    row[k] = v
            records.append(row)
        return records
    return [vars(table).copy()]


def records_to_table(records: Sequence[Mapping[str, Any]]) -> Record:
    keys: list[str] = []
    for rec in records:
        for key in rec.keys():
            if key not in keys:
                keys.append(key)
    return {key: [rec.get(key) for rec in records] for key in keys}


def field_vector(records: Sequence[Any], name: str, default: Any = None) -> list[Any]:
    return [get_field(rec, name, default) for rec in records]


def flatten(values: Iterable[Any]) -> list[Any]:
    out: list[Any] = []
    for item in values:
        if isinstance(item, np.ndarray):
            out.extend(item.ravel().tolist())
        elif isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out


def read_fasta(path: str | os.PathLike[str]) -> list[Record]:
    path = Path(path)
    records: list[Record] = []
    header = None
    seq_chunks: list[str] = []
    with path.open('r', encoding='utf-8') as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append({'Header': header, 'Sequence': ''.join(seq_chunks)})
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(re.sub(r'\s+', '', line))
    if header is not None:
        records.append({'Header': header, 'Sequence': ''.join(seq_chunks)})
    return records


def write_fasta_records(path: str | os.PathLike[str], records: Sequence[Mapping[str, Any]], append: bool = False, line_width: int = 70) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = 'a' if append else 'w'
    with path.open(mode, encoding='utf-8') as handle:
        for rec in records:
            header = str(rec.get('Header', rec.get('header', '')))
            seq = str(rec.get('Sequence', rec.get('sequence', '')))
            if not header.startswith('>'):
                header = '>' + header
            handle.write(header + '\n')
            compact = re.sub(r'\s+', '', seq)
            for start in range(0, len(compact), line_width):
                handle.write(compact[start:start + line_width] + '\n')
            handle.write('\n')


_DNA_COMP = str.maketrans('ACGTUacgtuNn', 'TGCAAtgcaaNn')


def reverse_complement(seq: str) -> str:
    return str(seq).translate(_DNA_COMP)[::-1]


def sanitize_locus_name(name: str) -> str:
    return re.sub(r'[\s]', '_', re.sub(r'[/\\]', lambda m: '[' if m.group(0) == '/' else ']', re.sub(r'[:]', 'C', re.sub(r'[.]', 'P', str(name).replace(',', '')))))


def bits_to_int(bits: Sequence[int] | str) -> int:
    if isinstance(bits, str):
        clean = re.sub(r'\s+', '', bits)
        if not clean:
            return 0
        return int(clean, 2)
    out = 0
    for bit in bits:
        out = (out << 1) | int(bit)
    return int(out)


def int_to_bits(value: int, width: int) -> np.ndarray:
    return np.array([(int(value) >> i) & 1 for i in range(width - 1, -1, -1)], dtype=np.uint8)


def generate_surrounding_codewords(codeword: Sequence[int] | str, distance: int) -> list[np.ndarray]:
    bits = np.array([int(x) for x in (list(codeword.strip()) if isinstance(codeword, str) else codeword)], dtype=np.uint8)
    if distance < 0:
        raise ValueError('distance must be non-negative')
    words: list[np.ndarray] = []
    for idxs in combinations(range(len(bits)), distance):
        mutated = bits.copy()
        for idx in idxs:
            mutated[idx] = 1 - mutated[idx]
        words.append(mutated)
    return words


def _required_parity_bits(k: int) -> int:
    r = 1
    while 2 ** r < k + r + 1:
        r += 1
    return r


def _all_messages(k: int) -> np.ndarray:
    return np.vstack([int_to_bits(i, k) for i in range(2 ** k)])


def _hamming_encode(messages: np.ndarray, code_length: int) -> np.ndarray:
    msg = np.asarray(messages, dtype=np.uint8)
    if msg.ndim == 1:
        msg = msg.reshape(1, -1)
    k = msg.shape[1]
    r = _required_parity_bits(k)
    parity_positions = {2 ** i for i in range(r)}
    data_positions = [pos for pos in range(1, code_length + 1) if pos not in parity_positions]
    code = np.zeros((msg.shape[0], code_length), dtype=np.uint8)
    for j, pos in enumerate(data_positions[:k]):
        code[:, pos - 1] = msg[:, j]
    for p in parity_positions:
        covered = [idx for idx in range(1, code_length + 1) if idx & p and idx != p]
        code[:, p - 1] = np.mod(code[:, np.asarray(covered) - 1].sum(axis=1), 2)
    return code


def gen_secded(num_letters: int, num_data_bits: int, on_bits: int | None = None) -> np.ndarray:
    """Python equivalent of the old GenSECDED use in deprecated scripts."""
    num_letters = int(num_letters)
    num_data_bits = int(num_data_bits)
    if num_letters == 12:
        msg = _all_messages(11)
        words = _hamming_encode(msg, 15)
        parity = np.mod(words.sum(axis=1), 2).reshape(-1, 1).astype(np.uint8)
        secded = np.concatenate([parity, words], axis=1)[0::16, :12]
    else:
        msg = _all_messages(num_data_bits)
        words = _hamming_encode(msg, num_letters - 1)
        parity = np.mod(words.sum(axis=1), 2).reshape(-1, 1).astype(np.uint8)
        secded = np.concatenate([parity, words], axis=1)
    if on_bits is not None:
        secded = secded[secded.sum(axis=1) == int(on_bits)]
    else:
        secded = secded[~((secded.sum(axis=1) <= 4) | (secded.sum(axis=1) >= num_letters - 4))]
    return secded.astype(np.uint8)


def pearson_corr(x: Sequence[float], y: Sequence[float]) -> tuple[float, float]:
    xarr = np.asarray(x, dtype=float)
    yarr = np.asarray(y, dtype=float)
    mask = np.isfinite(xarr) & np.isfinite(yarr)
    if mask.sum() < 2:
        return (float('nan'), float('nan'))
    try:
        from scipy.stats import pearsonr
        r, p = pearsonr(xarr[mask], yarr[mask])
        return float(r), float(p)
    except Exception:
        r = np.corrcoef(xarr[mask], yarr[mask])[0, 1]
        return float(r), float('nan')


def histogram(values: Sequence[float], bins: int | Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        if isinstance(bins, int):
            edges = np.linspace(0, 1, bins + 1)
        else:
            edges = np.asarray(bins, dtype=float)
        centers = (edges[:-1] + edges[1:]) / 2
        return np.zeros_like(centers), centers
    counts, edges = np.histogram(arr, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    return counts, centers


def save_pickle_and_json(obj: Any, base_path: str | os.PathLike[str]) -> tuple[Path, Path]:
    base = Path(base_path)
    base.parent.mkdir(parents=True, exist_ok=True)
    pkl_path = base.with_suffix('.pkl')
    json_path = base.with_suffix('.json')
    with pkl_path.open('wb') as handle:
        pickle.dump(obj, handle)
    def safe(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.integer):
            return int(value)
        if isinstance(value, np.floating):
            return float(value)
        if isinstance(value, Mapping):
            return {str(k): safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [safe(v) for v in value]
        return value
    with json_path.open('w', encoding='utf-8') as handle:
        json.dump(safe(obj), handle, indent=2)
    return pkl_path, json_path


def maybe_make_figure(name: str | None = None, visible: str = 'on'):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    if name is not None:
        try:
            fig.canvas.manager.set_window_title(str(name))
        except Exception:
            setattr(fig, '_merfish_name', str(name))
    if str(visible).lower() == 'off':
        plt.close(fig)
    return fig


def save_figure(fig, output_dir: str | os.PathLike[str], name: str, formats: Sequence[str] = ('png',), close: bool = False) -> list[Path]:
    import matplotlib.pyplot as plt
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for fmt in formats:
        if str(fmt).lower() == 'fig':
            continue
        path = output / f'{name}.{fmt}'
        fig.savefig(path, bbox_inches='tight')
        paths.append(path)
    if close:
        plt.close(fig)
    return paths


def run_command(command: str | Sequence[str], execute: bool = True, shell: bool = True, cwd: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    if not execute:
        return {'command': command, 'returncode': None, 'stdout': '', 'stderr': ''}
    result = subprocess.run(command, shell=shell, cwd=cwd, text=True, capture_output=True)
    return {'command': command, 'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr}
