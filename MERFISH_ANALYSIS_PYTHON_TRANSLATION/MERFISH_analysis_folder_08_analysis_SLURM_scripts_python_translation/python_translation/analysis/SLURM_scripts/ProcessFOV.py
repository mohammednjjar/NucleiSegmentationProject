"""Translate analysis/SLURM_scripts/ProcessFOV.m.

Purpose/use: warp and preprocess one MERFISH field of view.
"""

from __future__ import annotations

import argparse
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, get_fov_id, load_decoder, log_run_end, log_run_start, set_decoder_overwrite


def process_fov(
    base_path: str,
    array_id: Optional[int] = None,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> None:
    """Load a decoder and run WarpFOV(fov_id), then PreprocessFOV(fov_id)."""
    fov_id = get_fov_id(array_id)
    start = log_run_start("ProcessFOV", fov_id)
    decoder = load_decoder(base_path, decoder_loader=decoder_loader)
    set_decoder_overwrite(decoder, False)
    call_decoder_method(decoder, "WarpFOV", fov_id)
    call_decoder_method(decoder, "PreprocessFOV", fov_id)
    log_run_end(f"Completed warping and processing of fov {fov_id}", start)


def main() -> None:
    parser = argparse.ArgumentParser(description="Warp and preprocess one MERFISH FOV.")
    parser.add_argument("base_path")
    parser.add_argument("--array-id", type=int, default=None)
    args = parser.parse_args()
    process_fov(args.base_path, array_id=args.array_id)


if __name__ == "__main__":
    main()
