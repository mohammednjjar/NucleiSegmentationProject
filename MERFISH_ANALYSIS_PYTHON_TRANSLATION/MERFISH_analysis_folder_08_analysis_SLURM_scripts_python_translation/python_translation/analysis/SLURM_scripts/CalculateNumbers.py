"""Translate analysis/SLURM_scripts/CalculateNumbers.m.

Purpose/use: calculate barcode counts per feature and generate a feature-count report.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, diary, ensure_dir, load_decoder, log_run_end, log_run_start, set_parallel_if_available, slurm_ntasks


def calculate_numbers(
    nDPath: str,
    scratchPath: str,
    decoder_loader: Optional[Callable[[str], Any]] = None,
    workers: Optional[int] = None,
) -> None:
    """Run CalculateFeatureCounts and GenerateFeatureCountsReport on a decoder."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    if not scratchPath:
        raise ValueError("A path to scratch folder must be provided.")
    log_file = Path(nDPath) / "log" / "calculate.log"
    ensure_dir(log_file.parent)
    with diary(log_file):
        start = log_run_start("CalculateNumbers")
        print(f"Loading MERFISH Decoder from {nDPath}")
        decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
        set_parallel_if_available(decoder, workers if workers is not None else slurm_ntasks())
        call_decoder_method(decoder, "CalculateFeatureCounts")
        call_decoder_method(decoder, "GenerateFeatureCountsReport")
        log_run_end("Completed", start)


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate MERFISH barcode counts per feature.")
    parser.add_argument("nDPath")
    parser.add_argument("scratchPath")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    calculate_numbers(args.nDPath, args.scratchPath, workers=args.workers)


if __name__ == "__main__":
    main()
