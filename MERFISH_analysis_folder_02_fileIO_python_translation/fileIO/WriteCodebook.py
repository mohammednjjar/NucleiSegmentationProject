from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

from ._utils import _get_header_value, make_parameters, page_break


def WriteCodebook(codebookPath: str, barcodes, readouts: Sequence[Any], names: Sequence[str], ids: Sequence[str], **kwargs: Any):
    """Write a MERFISH codebook CSV.

    Translation note: the MATLAB source references undefined variables
    finalBarcodes/finalGenes/finalTargetRegions. This Python translation uses
    the documented function inputs: barcodes, readouts, names, and ids.
    """
    parameters = make_parameters({"verbose": True, "overwrite": False, "codebookName": "M3E1"}, kwargs)
    if any(x is None for x in [codebookPath, barcodes, readouts, names, ids]):
        raise ValueError("A valid path, barcodes, readouts, names, and ids are required")

    arr = np.asarray(barcodes)
    if arr.ndim != 2:
        raise ValueError("barcodes must be a 2-D matrix")
    if arr.shape[0] < len(names):
        raise ValueError("Insufficient barcodes provided for the provided names")
    if arr.shape[1] > len(readouts):
        raise ValueError("Insufficient readouts for the number of bits in the barcodes")
    if len(ids) < len(names):
        raise ValueError("Insufficient ids provided for the provided names")

    path = Path(codebookPath)
    if path.suffix.lower() != ".csv":
        import warnings
        warnings.warn("It is recommended to save codebooks as csv files", RuntimeWarning)
    if path.exists() and not parameters.overwrite:
        raise FileExistsError("Found existing codebook!")
    path.parent.mkdir(parents=True, exist_ok=True)

    if parameters.verbose:
        page_break()
        print(f"Writing codebook: {codebookPath}")

    bit_names = [_get_header_value(r) for r in readouts[: arr.shape[1]]]
    with path.open("w", newline="") as fh:
        fh.write("version, 1.0\n")
        fh.write(f"codebook_name, {parameters.codebookName}\n")
        fh.write("bit_names, " + ", ".join(bit_names) + "\n")
        fh.write("name, id, barcode\n")
        for i, name in enumerate(names):
            barcode = " ".join(str(int(x)) for x in arr[i, :])
            safe_id = str(ids[i]).replace(",", "_")
            fh.write(f"{name}, {safe_id}, {barcode}\n")

    if parameters.verbose:
        page_break()
        print(f"Wrote: {codebookPath}")
    return parameters
