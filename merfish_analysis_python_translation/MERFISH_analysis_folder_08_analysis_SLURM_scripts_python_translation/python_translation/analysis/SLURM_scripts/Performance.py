"""Translate analysis/SLURM_scripts/Performance.m.

Purpose/use: validate barcode files and run MERFISH performance metrics against
abundance data.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from .slurm_utils import (
    build_fov_barcode_files,
    decoder_normalized_data_path,
    decoder_report_path,
    diary,
    ensure_dir,
    load_decoder,
    log_run_end,
    log_run_start,
    nested_get,
    read_binary_file_header,
    set_parallel_if_available,
    slurm_ntasks,
)


def _default_metrics_runner(nDPath: str, parameters: Dict[str, Any]) -> Any:
    """Import and call a translated MERFISHPerformanceMetrics function."""
    candidates = [
        "MERFISHPerformanceMetrics",
        "analysis.functions.MERFISHPerformanceMetrics",
        "python_translation.analysis.functions.MERFISHPerformanceMetrics",
    ]
    last_error: Optional[Exception] = None
    for module_name in candidates:
        try:
            module = __import__(module_name, fromlist=["*"])
        except ImportError as exc:
            last_error = exc
            continue
        for attr in ("MERFISHPerformanceMetrics", "merfish_performance_metrics", "run"):
            function = getattr(module, attr, None)
            if callable(function):
                return function(nDPath, parameters=parameters)
    raise RuntimeError(
        "Could not import a translated MERFISHPerformanceMetrics function. "
        "Pass metrics_runner=callable or add the translated analysis/functions package."
    ) from last_error


def performance(
    nDPath: str,
    aPath: str,
    scratchPath: str,
    decoder_loader: Optional[Callable[[str], Any]] = None,
    metrics_runner: Optional[Callable[[str, Dict[str, Any]], Any]] = None,
    workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Run the Python equivalent of the MATLAB Performance.m workflow."""
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    if not aPath:
        raise ValueError("A path to abundance data must be provided.")
    if not scratchPath:
        raise ValueError("A path to scratch folder must be provided.")

    log_file = Path(nDPath) / "log" / "performance.log"
    ensure_dir(log_file.parent)
    with diary(log_file):
        start = log_run_start("Performance")
        print(f"Loading MERFISH Decoder from {nDPath}")
        decoder = load_decoder(nDPath, decoder_loader=decoder_loader)

        barcode_files = build_fov_barcode_files(nDPath)
        print(f"Found {len(barcode_files)} barcode files")
        corrupt = []
        for record in barcode_files:
            header = read_binary_file_header(record["filePath"])
            if header["isCorrupt"]:
                print(f"Barcode file for fov {record['fov']} is corrupt")
                corrupt.append(record["fov"])
        if not corrupt:
            print("All barcode files appear to be complete and uncorrupted!")

        decoder_ids = set(int(v) for v in getattr(decoder, "fovIDs", getattr(decoder, "fov_ids", [])))
        file_ids = set(int(record["fov"]) for record in barcode_files)
        missing_fov_ids = sorted(decoder_ids.difference(file_ids))
        if missing_fov_ids:
            print(f"Discovered missing fov ids: {missing_fov_ids}")

        parallel_workers = workers if workers is not None else slurm_ntasks()
        set_parallel_if_available(decoder, parallel_workers)
        parameters = {
            "parallel": parallel_workers,
            "brightnessThreshold": nested_get(decoder, "parameters.quantification.minimumBarcodeBrightness"),
            "areaThreshold": nested_get(decoder, "parameters.quantification.minimumBarcodeArea"),
            "stageOrientation": nested_get(decoder, "parameters.decoding.stageOrientation"),
            "abundDataPath": aPath,
            "verbose": True,
            "outputPath": str(Path(decoder_normalized_data_path(decoder, nDPath)) / decoder_report_path(decoder) / "performance"),
        }
        print("Using the following parameters")
        for key, value in parameters.items():
            if key != "parallel":
                print(f" {key}: {value}")

        runner = metrics_runner or _default_metrics_runner
        result = runner(nDPath, parameters)
        log_run_end("Completed", start)
        return {
            "result": result,
            "corrupt_fov_ids": corrupt,
            "missing_fov_ids": missing_fov_ids,
            "parameters": parameters,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MERFISH performance metrics.")
    parser.add_argument("nDPath")
    parser.add_argument("aPath")
    parser.add_argument("scratchPath")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    performance(args.nDPath, args.aPath, args.scratchPath, workers=args.workers)


if __name__ == "__main__":
    main()
