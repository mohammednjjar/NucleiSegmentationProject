"""Translate analysis/SLURM_scripts/CalculateDoubletScore.m.

Purpose/use: calculate per-feature barcode counts, barcode centers of mass, and
coordinate variances used to identify potential segmentation/doublet errors.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

import numpy as np
import pandas as pd

from .slurm_utils import diary, ensure_dir, load_decoder, log_run_end, log_run_start


def _feature_attr(feature: Any, name: str) -> Any:
    if isinstance(feature, dict):
        return feature[name]
    return getattr(feature, name)


def _extract_features(decoder: Any) -> Sequence[Any]:
    method = getattr(decoder, "GetFoundFeatures", None)
    if callable(method):
        return method()
    value = getattr(decoder, "foundFeatures", None)
    if value is None:
        value = getattr(decoder, "found_features", None)
    if value is None:
        raise AttributeError("Decoder has no GetFoundFeatures(), foundFeatures, or found_features.")
    return value


def calculate_doublet_score(
    nDPath: str,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> dict[str, Path]:
    """Compute and export barcode count/center/variance tables.

    Output files match the MATLAB script names:
    barcode_counts.csv, barcode_center_of_mass_X.csv,
    barcode_center_of_mass_Y.csv, barcode_center_of_mass_var_X.csv,
    barcode_center_of_mass_var_Y.csv.
    """
    if not nDPath:
        raise ValueError("A normalized data path must be provided.")
    npath = Path(nDPath)
    log_file = npath / "log" / "calculate.log"
    ensure_dir(log_file.parent)
    with diary(log_file):
        start = log_run_start("CalculateDoubletScore")
        print(f"Loading MERFISH Decoder from {nDPath}")
        decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
        found_features = list(_extract_features(decoder))
        feature_ids = [int(_feature_attr(feature, "feature_id")) for feature in found_features]
        if not feature_ids:
            raise ValueError("No found features were returned by the decoder.")

        barcode_metadata_path = npath / "reports" / "barcode_metadata.csv"
        print("Loading barcode metadata")
        barcode_table = pd.read_csv(barcode_metadata_path)
        print(f"...found {len(barcode_table)} barcodes")
        barcode_table = barcode_table.loc[barcode_table["in_feature"] == 1].copy()
        print(f"...cutting to {len(barcode_table)} barcodes within features")

        num_barcodes = len(getattr(decoder, "codebook"))
        max_feature_id = max(feature_ids)
        counts = np.zeros((max_feature_id, num_barcodes), dtype=float)
        com_x = np.full((max_feature_id, num_barcodes), np.nan, dtype=float)
        com_y = np.full((max_feature_id, num_barcodes), np.nan, dtype=float)
        var_x = np.full((max_feature_id, num_barcodes), np.nan, dtype=float)
        var_y = np.full((max_feature_id, num_barcodes), np.nan, dtype=float)

        required = {"feature_id", "barcode_id", "abs_position_1", "abs_position_2"}
        missing = sorted(required.difference(barcode_table.columns))
        if missing:
            raise ValueError(f"barcode_metadata.csv is missing required columns: {missing}")

        for feature_id, local in barcode_table.groupby("feature_id", sort=True):
            row = int(feature_id) - 1
            if row < 0 or row >= max_feature_id:
                continue
            for barcode_id, sub in local.groupby("barcode_id", sort=True):
                col = int(barcode_id) - 1
                if col < 0 or col >= num_barcodes:
                    continue
                counts[row, col] = float(len(sub))
                com_x[row, col] = float(sub["abs_position_1"].mean())
                com_y[row, col] = float(sub["abs_position_2"].mean())
                var_x[row, col] = float(sub["abs_position_1"].var(ddof=1)) if len(sub) > 1 else 0.0
                var_y[row, col] = float(sub["abs_position_2"].var(ddof=1)) if len(sub) > 1 else 0.0

        save_path = ensure_dir(npath / "reports")
        outputs = {
            "counts": save_path / "barcode_counts.csv",
            "com_x": save_path / "barcode_center_of_mass_X.csv",
            "com_y": save_path / "barcode_center_of_mass_Y.csv",
            "var_x": save_path / "barcode_center_of_mass_var_X.csv",
            "var_y": save_path / "barcode_center_of_mass_var_Y.csv",
        }
        np.savetxt(outputs["counts"], counts, delimiter=",")
        np.savetxt(outputs["com_x"], com_x, delimiter=",")
        np.savetxt(outputs["com_y"], com_y, delimiter=",")
        np.savetxt(outputs["var_x"], var_x, delimiter=",")
        np.savetxt(outputs["var_y"], var_y, delimiter=",")
        log_run_end("Completed", start)
        return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate barcode doublet-score metadata.")
    parser.add_argument("nDPath")
    args = parser.parse_args()
    calculate_doublet_score(args.nDPath)


if __name__ == "__main__":
    main()
