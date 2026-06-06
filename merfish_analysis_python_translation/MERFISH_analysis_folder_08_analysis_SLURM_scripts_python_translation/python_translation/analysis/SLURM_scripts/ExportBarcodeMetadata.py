"""Translate analysis/SLURM_scripts/ExportBarcodeMetadata.m.

Purpose/use: export barcode metadata, then compute doublet-score metadata tables.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable, Optional

from .CalculateDoubletScore import calculate_doublet_score
from .slurm_utils import call_decoder_method, diary, ensure_dir, load_decoder, log_run_end, log_run_start, set_parallel_if_available, slurm_ntasks


def export_barcode_metadata(
    nDPath: str,
    decoder_loader: Optional[Callable[[str], Any]] = None,
    workers: Optional[int] = None,
) -> None:
    """Run BarcodesToCSV, then calculate doublet-score barcode metadata."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    log_file = Path(nDPath) / "log" / "calculate.log"
    ensure_dir(log_file.parent)
    with diary(log_file):
        start = log_run_start("ExportBarcodeMetadata")
        print(f"Loading MERFISH Decoder from {nDPath}")
        decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
        set_parallel_if_available(decoder, workers if workers is not None else slurm_ntasks())
        call_decoder_method(decoder, "BarcodesToCSV")
        calculate_doublet_score(nDPath, decoder_loader=decoder_loader)
        log_run_end("Completed", start)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export barcode metadata and doublet-score tables.")
    parser.add_argument("nDPath")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    export_barcode_metadata(args.nDPath, workers=args.workers)


if __name__ == "__main__":
    main()
