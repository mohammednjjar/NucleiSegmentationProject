"""Translate analysis/SLURM_scripts/CombineBoundaries.m.

Purpose/use: combine segmentation boundaries/features after FOV-level segmentation.
"""

from __future__ import annotations

import argparse
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, load_decoder, log_run_end, log_run_start, set_decoder_overwrite


def combine_boundaries(
    nDPath: str,
    overwrite: bool = False,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> Any:
    """Load a decoder and run CombineFeatures()."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    start = log_run_start("CombineBoundaries")
    decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
    set_decoder_overwrite(decoder, overwrite)
    result = call_decoder_method(decoder, "CombineFeatures")
    log_run_end("Completed", start)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine MERFISH boundaries/features.")
    parser.add_argument("nDPath")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    combine_boundaries(args.nDPath, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
