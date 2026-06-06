from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import Any, Callable, Sequence

from ._utils import _apply_conv, _split_preserve_delimiters, make_parameters


def BuildFileStructure(folderPath: str, **kwargs: Any):
    """Parse file names in a folder into a list of dictionaries.

    Python translation of BuildFileStructure.m. It supports extension filtering,
    delimiter-based parsing, regex named-group parsing, required/excluded flags,
    and field conversion functions.
    """
    defaults = {
        "fileExt": "*",
        "delimiters": ["_"],
        "fieldNames": ["field1"],
        "fieldConv": [str],
        "appendExtraFields": False,
        "excludeFlags": [],
        "requireFlag": "",
        "requireExactMatch": False,
        "containsDelimiters": None,
        "regExp": None,
    }
    parameters = make_parameters(defaults, kwargs)

    folder = Path(folderPath)
    if not folder.is_dir():
        raise ValueError("A valid folder is required.")

    if any(name in {"name", "filePath"} for name in parameters.fieldNames):
        raise ValueError("name and filePath are protected names.")

    if parameters.containsDelimiters is None:
        parameters.containsDelimiters = len(parameters.fieldConv)
    if isinstance(parameters.containsDelimiters, Sequence) and not isinstance(parameters.containsDelimiters, (str, bytes)):
        if len(parameters.containsDelimiters) > 1:
            raise ValueError("Only one field can have internal delimiters.")
        if len(parameters.containsDelimiters) == 1:
            parameters.containsDelimiters = int(parameters.containsDelimiters[0])
    parameters.containsDelimiters = int(parameters.containsDelimiters)

    while len(parameters.fieldConv) < len(parameters.fieldNames):
        parameters.fieldConv.append(str)

    ext = str(parameters.fileExt).lstrip(".")
    pattern = "*" if ext in {"*", ""} else f"*.{ext}"
    files = [p for p in sorted(folder.glob(pattern)) if p.is_file()]

    if parameters.requireFlag:
        files = [p for p in files if str(parameters.requireFlag) in p.name]
    if parameters.excludeFlags:
        files = [p for p in files if not any(flag in p.name for flag in parameters.excludeFlags)]

    parsed: list[dict[str, Any]] = []
    for path in files:
        if parameters.regExp:
            match = re.search(parameters.regExp, path.name)
            if not match:
                continue
            local = dict(match.groupdict())
            for idx, field in enumerate(parameters.fieldNames):
                if field not in local or local[field] is None:
                    local[field] = None
                else:
                    local[field] = _apply_conv(local[field], parameters.fieldConv[idx])
            local["name"] = path.name
            local["filePath"] = str(path)
            local["regExp"] = parameters.regExp
            parsed.append(local)
            continue

        split_text, found_delimiters = _split_preserve_delimiters(path.name, list(parameters.delimiters) + ["."])
        if split_text:
            split_text = split_text[:-1]
        if found_delimiters:
            found_delimiters = found_delimiters[:-1]

        if len(split_text) > len(parameters.fieldNames) and parameters.appendExtraFields:
            start = parameters.containsDelimiters - 1
            length_diff = len(split_text) - len(parameters.fieldNames) + 1
            finish = len(split_text) - (len(parameters.fieldNames) - parameters.containsDelimiters)
            combined: list[str] = []
            for i in range(start, finish):
                combined.append(split_text[i])
                if i < finish - 1 and i < len(found_delimiters):
                    combined.append(found_delimiters[i])
            old = split_text
            split_text = old[:start] + ["".join(combined)] + old[start + length_diff:]

        condition = len(split_text) == len(parameters.fieldNames) if parameters.requireExactMatch else len(split_text) <= len(parameters.fieldNames)
        if not condition:
            continue

        local = {"name": path.name, "filePath": str(path)}
        for j, text in enumerate(split_text):
            local[parameters.fieldNames[j]] = _apply_conv(text, parameters.fieldConv[j])
        local["delimiters"] = list(parameters.delimiters) + ["."]
        parsed.append(local)

    if len(parsed) < len(files):
        warnings.warn("Some files did not fit the specified pattern.", RuntimeWarning)
    return parsed, parameters
