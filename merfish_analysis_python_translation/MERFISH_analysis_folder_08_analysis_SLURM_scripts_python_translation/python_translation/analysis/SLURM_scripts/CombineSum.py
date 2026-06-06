"""Translate analysis/SLURM_scripts/CombineSum.m.

Purpose/use: combine raw-signal summation outputs and generate a summation report.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, diary, ensure_dir, load_decoder, log_run_end, log_run_start, set_decoder_overwrite


def combine_sum(
    nDPath: str,
    overwrite: bool = False,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> None:
    """Run CombineRawSum and GenerateSummationReport on a decoder."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    log_file = Path(nDPath) / "log" / "combine_sum.log"
    ensure_dir(log_file.parent)
    with diary(log_file):
        start = log_run_start("CombineSum")
        print(f"Loading MERFISH Decoder from {nDPath}")
        decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
        set_decoder_overwrite(decoder, overwrite)
        call_decoder_method(decoder, "CombineRawSum")
        call_decoder_method(decoder, "GenerateSummationReport")
        log_run_end("Completed", start)


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine raw MERFISH signal sums.")
    parser.add_argument("nDPath")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    combine_sum(args.nDPath, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
