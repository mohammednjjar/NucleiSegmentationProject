from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from ._utils import _parse_codebook_header_rows, _read_csv_rows, make_parameters, page_break


def LoadCodebook(codebookPath: str, **kwargs: Any):
    """Load a MERFISH codebook CSV into (codebook, header, parameters)."""
    defaults = {"verbose": True, "barcodeConvFunc": str}
    parameters = make_parameters(defaults, kwargs)
    path = Path(codebookPath)
    if not path.exists():
        raise ValueError("A valid path to a codebook is required")

    if parameters.verbose:
        page_break()
        print(f"Loading codebook from: {codebookPath}")

    rows = _read_csv_rows(path)
    header, data_start = _parse_codebook_header_rows(rows)

    if parameters.verbose:
        for key, value in header.items():
            display = ", ".join(value) if isinstance(value, list) else value
            print(f"...{key}: {display}")

    version = str(header["version"])
    if version not in {"1.0", "1"}:
        raise ValueError(f"Unsupported codebook version: {version}")

    rmspace = lambda x: "".join(str(x).split())
    codebook = []
    for row in rows[data_start:]:
        if not row or all(str(x).strip() == "" for x in row):
            continue
        if len(row) < 3:
            raise ValueError(f"Invalid codebook data row: {row}")
        barcode_text = rmspace(row[2])
        codebook.append({
            "name": row[0].strip(),
            "id": row[1].strip(),
            "barcode": parameters.barcodeConvFunc(barcode_text),
        })

    if parameters.verbose:
        print(f"...loaded {len(codebook)} barcodes")
    return codebook, header, parameters
