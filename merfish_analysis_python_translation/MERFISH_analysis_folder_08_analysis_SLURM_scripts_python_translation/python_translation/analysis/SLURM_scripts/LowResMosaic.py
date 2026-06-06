"""Translate analysis/SLURM_scripts/LowResMosaic.m.

Purpose/use: generate low-resolution MERFISH mosaic images from the decoder output.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, diary, ensure_dir, load_decoder, log_run_end, log_run_start


def low_res_mosaic(
    nDPath: str,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> Any:
    """Load a decoder and run GenerateLowResolutionMosaic()."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    log_file = Path(nDPath) / "log" / "low_res_mosaics.log"
    ensure_dir(log_file.parent)
    with diary(log_file):
        start = log_run_start("LowResMosaic")
        print(f"Loading MERFISH Decoder from {nDPath}")
        decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
        result = call_decoder_method(decoder, "GenerateLowResolutionMosaic")
        log_run_end("Completed", start)
        return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate MERFISH low-resolution mosaics.")
    parser.add_argument("nDPath")
    args = parser.parse_args()
    low_res_mosaic(args.nDPath)


if __name__ == "__main__":
    main()
