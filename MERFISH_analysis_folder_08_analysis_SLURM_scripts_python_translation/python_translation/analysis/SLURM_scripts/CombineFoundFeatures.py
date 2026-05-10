"""Translate analysis/SLURM_scripts/CombineFoundFeatures.m.

Purpose/use: combine found features, generate the feature report, and export features to CSV.
"""

from __future__ import annotations

import argparse
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, load_decoder, log_run_end, log_run_start, set_decoder_overwrite


def combine_found_features(
    nDPath: str,
    overwrite: bool = False,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> None:
    """Run CombineFeatures, GenerateFoundFeaturesReport, and FoundFeaturesToCSV."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    start = log_run_start("CombineFoundFeatures")
    decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
    set_decoder_overwrite(decoder, overwrite)
    call_decoder_method(decoder, "CombineFeatures")
    call_decoder_method(decoder, "GenerateFoundFeaturesReport")
    call_decoder_method(decoder, "FoundFeaturesToCSV", downSampleFactor=10, zIndex=4)
    log_run_end("Completed", start)


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine found MERFISH features.")
    parser.add_argument("nDPath")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    combine_found_features(args.nDPath, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
