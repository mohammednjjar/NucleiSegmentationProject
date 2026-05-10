from pathlib import Path
import tempfile
import numpy as np
import pandas as pd

from analysis.functions import MLists2Transform, GenerateGeoTransformReport, MERFISHPerformanceMetrics


def test_mlists_to_transform_identical_points():
    ref = pd.DataFrame({
        "x": [0.0, 1.0, 0.0, 1.0],
        "y": [0.0, 0.0, 1.0, 1.0],
        "frame": [1, 1, 1, 1],
    })
    mov = pd.DataFrame({
        "x": [0.0, 1.0, 0.0, 1.0],
        "y": [0.0, 0.0, 1.0, 1.0],
        "frame": [1, 1, 1, 1],
        "xc": [0.0, 0.0, 0.0, 0.0],
        "yc": [0.0, 0.0, 0.0, 0.0],
    })
    tforms, updated, residuals, inds, params = MLists2Transform(ref, mov)
    assert len(tforms) == 1
    assert np.allclose(updated["xc"].to_numpy(), mov["x"].to_numpy())
    assert np.allclose(updated["yc"].to_numpy(), mov["y"].to_numpy())
    assert residuals[0].shape[1] == 4
    assert inds[0].shape[1] == 2


def test_generate_geo_transform_report():
    ref = pd.DataFrame({"x": [0, 1], "y": [0, 0], "frame": [1, 1]})
    mov = pd.DataFrame({"x": [0, 1], "y": [0, 0], "frame": [1, 1], "xc": [0, 0], "yc": [0, 0]})
    tforms, _, residuals, _, _ = MLists2Transform(ref, mov)
    report, figs, params = GenerateGeoTransformReport(tforms, residuals, makeFigures=False)
    assert "muX" in report
    assert "errX" in report
    assert "offsetX" in report


def test_merfish_performance_metrics_csv_input():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        (root / "barcodes" / "barcode_fov").mkdir(parents=True)
        pd.DataFrame({
            "name": ["GeneA", "Blank-1", "GeneB"],
            "bit1": [1, 0, 1],
            "bit2": [0, 1, 1],
        }).to_csv(root / "test_codebook.csv", index=False)
        pd.DataFrame({
            "barcode_id": [1, 1, 2, 3],
            "is_exact": [1, 0, 1, 1],
            "error_dir": [0, 1, 0, 0],
            "error_bit": [0, 1, 0, 2],
            "total_magnitude": [20, 30, 40, 50],
            "area": [2, 3, 4, 5],
            "fov_id": [1, 1, 1, 1],
            "x": [10, 20, 30, 40],
            "y": [5, 15, 25, 35],
            "z": [1, 1, 1, 1],
        }).to_csv(root / "barcodes" / "barcode_fov" / "fov_1.csv", index=False)
        result = MERFISHPerformanceMetrics(
            root,
            verbose=False,
            logProgress=True,
            archive=True,
            visibleOption="off",
            formats=("png",),
            imageSize=(64, 64),
            numZPos=1,
        )
        assert result["countsPerCellExact"].shape == (3, 1)
        assert result["countsPerCellCorrected"].shape == (3, 1)
        assert (root / "barcodes" / "performance" / "countsPerCellExact.csv").exists()
        assert (root / "barcodes" / "performance" / "performance.log").exists()
