"""Python translation of `deprecated/misc/StringFind.m`."""
from __future__ import annotations

from typing import Iterable, Sequence, Any


def StringFind(cellString: Sequence[str], targetString: str | Sequence[str], exactly: bool = False,
               cellOutput: bool = False, boolean: bool = False):
    if cellString is None or targetString is None:
        raise ValueError('cellString and targetString are required')
    haystack = [str(x) for x in cellString]

    def search_one(target: str):
        if exactly:
            matches = [i for i, s in enumerate(haystack) if s == str(target)]
        else:
            matches = [i for i, s in enumerate(haystack) if str(target) in s]
        return matches

    if isinstance(targetString, (str, bytes)):
        idx = search_one(str(targetString))
        if boolean:
            mask = [False] * len(haystack)
            for i in idx:
                mask[i] = True
            return mask, []
        return idx, []

    target_list = [str(t) for t in targetString]
    grouped = [search_one(t) for t in target_list]
    not_idx = [i for i, matches in enumerate(grouped) if len(matches) == 0]
    if boolean:
        mask = [False] * len(haystack)
        for matches in grouped:
            for i in matches:
                mask[i] = True
        return mask, not_idx
    if cellOutput:
        return grouped, not_idx
    if grouped and max(len(x) for x in grouped) <= 1:
        return [x[0] for x in grouped if x], not_idx
    return grouped, not_idx


def string_find(cell_string, target_string, exactly: bool = False, cell_output: bool = False, as_boolean: bool = False):
    return StringFind(cell_string, target_string, exactly=exactly, cellOutput=cell_output, boolean=as_boolean)
