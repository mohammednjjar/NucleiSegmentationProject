"""Translate analysis/SLURM_scripts/Optimize.m.

Purpose/use: generate warp report, initialize scale factors, and optimize decoding scale factors.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable, Optional

from .slurm_utils import call_decoder_method, diary, ensure_dir, load_decoder, log_run_end, log_run_start, set_parallel_if_available, slurm_ntasks


def optimize(
    nDPath: str,
    overwrite: bool = False,
    decoder_loader: Optional[Callable[[str], Any]] = None,
    workers: Optional[int] = None,
    iterations: int = 25,
    use_blanks: bool = False,
) -> None:
    """Run the MATLAB Optimize.m sequence on a translated decoder."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    log_file = Path(nDPath) / "log" / "optimize.log"
    ensure_dir(log_file.parent)
    with diary(log_file):
        start = log_run_start("Optimize")
        print(f"Loading MERFISH Decoder from {nDPath}")
        decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
        call_decoder_method(decoder, "GenerateWarpReport")
        call_decoder_method(decoder, "InitializeScaleFactors")
        set_parallel_if_available(decoder, workers if workers is not None else slurm_ntasks())
        call_decoder_method(decoder, "OptimizeScaleFactors", iterations, overwrite=overwrite, useBlanks=use_blanks)
        save_method = getattr(decoder, "Save", None)
        if callable(save_method):
            save_method()
        log_run_end("Completed", start)


def main() -> None:
    parser = argparse.ArgumentParser(description="Optimize MERFISH decoding scale factors.")
    parser.add_argument("nDPath")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--iterations", type=int, default=25)
    parser.add_argument("--use-blanks", action="store_true")
    args = parser.parse_args()
    optimize(args.nDPath, overwrite=args.overwrite, workers=args.workers, iterations=args.iterations, use_blanks=args.use_blanks)


if __name__ == "__main__":
    main()
