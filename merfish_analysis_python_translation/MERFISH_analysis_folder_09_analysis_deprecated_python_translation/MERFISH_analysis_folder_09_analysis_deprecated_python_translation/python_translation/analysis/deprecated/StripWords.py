
from __future__ import annotations

from .helpers import fields, parse_parameters, remove_fields


def strip_words(words, *, parameters=None, fieldsToKeep=None):
    """Python translation of `StripWords.m`.

    Keeps only selected fields in each word structure.
    """
    keep = fieldsToKeep or [
        "intCodeword", "geneName", "isExactMatch", "isCorrectedMatch",
        "imageX", "imageY", "wordCentroidX", "wordCentroidY", "cellID",
    ]
    params = parse_parameters({"fieldsToKeep": keep}, parameters)
    out = list(words)
    for word in out:
        existing = set(fields(word))
        fields_to_remove = existing.difference(params["fieldsToKeep"])
        remove_fields(word, fields_to_remove)
    return out


# MATLAB-style alias
StripWords = strip_words
