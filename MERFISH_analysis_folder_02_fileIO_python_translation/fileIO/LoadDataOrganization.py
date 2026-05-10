from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from ._utils import make_parameters, page_break


def _coerce_cell(value: str):
    text = str(value).strip()
    if text == "":
        return ""
    for conv in (int, float):
        try:
            return conv(text)
        except ValueError:
            continue
    return text


def _parse_internal_numeric(value: Any, delimiter: str):
    if not isinstance(value, str):
        return value
    parts = [p.strip() for p in value.split(delimiter)]
    parsed = []
    for part in parts:
        if part == "":
            continue
        try:
            x = float(part)
            parsed.append(int(x) if x.is_integer() else x)
        except ValueError:
            return value
    return parsed


def LoadDataOrganization(dataOrgPath: str, **kwargs: Any):
    """Load a MERFISH data-organization CSV and build metadata."""
    parameters = make_parameters({"verbose": False, "internalDelimiter": ";"}, kwargs)
    path = Path(dataOrgPath)
    if not path.exists():
        raise ValueError("A valid path to a data organization file must be provided")
    if path.suffix.lower() != ".csv":
        raise ValueError("Only csv files are supported.")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        data_org = []
        for row in reader:
            data_org.append({k: _coerce_cell(v) for k, v in row.items()})

    for field in ["frame", "zPos"]:
        for record in data_org:
            if field in record:
                record[field] = _parse_internal_numeric(record[field], parameters.internalDelimiter)

    if "zPos" not in data_org[0] if data_org else True:
        for record in data_org:
            record["zPos"] = 0
        unique_z = [0]
    else:
        flat = []
        for record in data_org:
            z = record.get("zPos", 0)
            flat.extend(z if isinstance(z, list) else [z])
        unique_z = sorted(set(flat)) if flat else [0]

    meta_data = {
        "numDataChannels": len(data_org),
        "zPos": unique_z,
        "numZPos": len(unique_z),
    }

    if parameters.verbose:
        page_break()
        print(f"Loaded data organization file: {dataOrgPath}")
        print(f"Found {meta_data['numDataChannels']} data channels")
        print(f"Found {meta_data['numZPos']} z-stacks")
    return data_org, meta_data, parameters
