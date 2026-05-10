"""Translate analysis/SLURM_scripts/ParseFOV.m.

Purpose/use: parse decoded barcodes into segmentation feature boundaries for one FOV.
"""

from __future__ import annotations

import argparse
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, get_fov_id, load_decoder, log_run_end, log_run_start, set_decoder_overwrite


def parse_fov(
    base_path: str,
    array_id: Optional[int] = None,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> Any:
    """Load a decoder from base_path and run ParseFOV(fov_id)."""
    fov_id = get_fov_id(array_id)
    start = log_run_start("ParseFOV", fov_id)
    decoder = load_decoder(base_path, decoder_loader=decoder_loader)
    set_decoder_overwrite(decoder, False)
    result = call_decoder_method(decoder, "ParseFOV", fov_id)
    log_run_end(f"Completed parsing of fov {fov_id}", start)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse one MERFISH FOV.")
    parser.add_argument("base_path")
    parser.add_argument("--array-id", type=int, default=None)
    args = parser.parse_args()
    parse_fov(args.base_path, array_id=args.array_id)


if __name__ == "__main__":
    main()
