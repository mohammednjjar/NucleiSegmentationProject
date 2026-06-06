"""Python translation of `deprecated/probe_construction/WriteFasta.m`."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

try:
    from ..helpers import table_to_records, write_fasta_records
except ImportError:
    from deprecated.helpers import table_to_records, write_fasta_records


def WriteFasta(filename: str, Header: Any, seq: Any = None, Append: bool = False, Warnings: bool = True):
    path = Path(filename)
    if path.exists() and Warnings:
        action = 'appended to' if Append else 'overwritten'
        print(f'existing file {path} will be {action}')
    if isinstance(Header, Mapping) or isinstance(Header, list) and Header and isinstance(Header[0], Mapping):
        records = table_to_records(Header)
    elif isinstance(Header, (str, bytes)):
        records = [{'Header': str(Header), 'Sequence': '' if seq is None else str(seq)}]
    else:
        headers = list(Header)
        sequences = list(seq)
        records = [{'Header': h, 'Sequence': s} for h, s in zip(headers, sequences)]
    write_fasta_records(path, records, append=Append)
    return str(path)


def write_fasta(filename: str, header: Any, seq: Any = None, append: bool = False, warnings: bool = True):
    return WriteFasta(filename, header, seq, Append=append, Warnings=warnings)
