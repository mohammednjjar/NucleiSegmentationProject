
from __future__ import annotations

from .helpers import get_field, key_converter, parse_parameters, set_field


def decode_words(words, exactMap, correctableMap, *, parameters=None, keyType="binStr"):
    """Python translation of `DecodeWords.m`.

    Decodes each word using exact and correctable maps. Exact matches set `isExactMatch`; correctable matches set `isCorrectedMatch` and can overwrite `geneName`, matching the MATLAB order.
    """
    params = parse_parameters({"keyType": keyType}, parameters)
    conv = key_converter(params["keyType"])
    out = list(words)
    exact = exactMap or {}
    correctable = correctableMap or {}
    for word in out:
        key = conv(get_field(word, "codeword"))
        if key in exact:
            set_field(word, "geneName", exact[key])
            set_field(word, "isExactMatch", True)
        if key in correctable:
            set_field(word, "geneName", correctable[key])
            set_field(word, "isCorrectedMatch", True)
    return out, params


# MATLAB-style alias
DecodeWords = decode_words
