"""Translate analysis/SLURM_scripts/Segment.m.

Purpose/use: segment all FOVs using a saved MERFISHDecoder.
"""

from __future__ import annotations

import argparse
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, load_decoder, log_run_end, log_run_start, set_parallel_if_available, slurm_ntasks


def segment(
    nDPath: str,
    decoder_loader: Optional[Callable[[str], Any]] = None,
    workers: Optional[int] = None,
) -> Any:
    """Run SegmentFOV(None) to segment all FOVs, matching MATLAB SegmentFOV([])."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    start = log_run_start("Segment")
    decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
    set_parallel_if_available(decoder, workers if workers is not None else slurm_ntasks())
    result = call_decoder_method(decoder, "SegmentFOV", None)
    log_run_end("Completed segmentation of all fovs", start)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Segment all MERFISH FOVs.")
    parser.add_argument("nDPath")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    segment(args.nDPath, workers=args.workers)


if __name__ == "__main__":
    main()
