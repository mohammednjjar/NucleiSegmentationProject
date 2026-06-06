from __future__ import annotations

import csv
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "python_translation"
sys.path.insert(0, str(ROOT))

from analysis.SLURM_scripts.CalculateDoubletScore import calculate_doublet_score
from analysis.SLURM_scripts.CombineFoundFeatures import combine_found_features
from analysis.SLURM_scripts.DecodeFOV import decode_fov
from analysis.SLURM_scripts.MERFISHScheduler import merfish_scheduler
from analysis.SLURM_scripts.Optimize import optimize
from analysis.SLURM_scripts.ProcessFOV import process_fov


class Feature:
    def __init__(self, feature_id, uid):
        self.feature_id = feature_id
        self.uID = uid


class FakeDecoder:
    def __init__(self):
        self.calls = []
        self.codebook = ["a", "b", "c"]
        self.fovIDs = [1, 2]
        self.overwrite = None
        self.parameters = {
            "quantification": {"minimumBarcodeBrightness": 10, "minimumBarcodeArea": 2},
            "decoding": {"stageOrientation": "normal"},
        }
        self.normalizedDataPath = ""
        self.reportPath = "reports"

    def _call(self, name, *args, **kwargs):
        self.calls.append((name, args, kwargs))
        return name

    def DecodeFOV(self, fov):
        return self._call("DecodeFOV", fov)

    def WarpFOV(self, fov):
        return self._call("WarpFOV", fov)

    def PreprocessFOV(self, fov):
        return self._call("PreprocessFOV", fov)

    def CombineFeatures(self):
        return self._call("CombineFeatures")

    def GenerateFoundFeaturesReport(self):
        return self._call("GenerateFoundFeaturesReport")

    def FoundFeaturesToCSV(self, **kwargs):
        return self._call("FoundFeaturesToCSV", **kwargs)

    def GenerateWarpReport(self):
        return self._call("GenerateWarpReport")

    def InitializeScaleFactors(self):
        return self._call("InitializeScaleFactors")

    def SetParallel(self, workers):
        return self._call("SetParallel", workers)

    def OptimizeScaleFactors(self, iterations, **kwargs):
        return self._call("OptimizeScaleFactors", iterations, **kwargs)

    def Save(self):
        return self._call("Save")

    def GetFoundFeatures(self):
        return [Feature(1, "f1"), Feature(2, "f2")]


def test_fov_driver_calls_decoder_methods():
    dec = FakeDecoder()
    decode_fov("/tmp/fake", array_id=7, decoder_loader=lambda path: dec)
    assert dec.calls[-1][0] == "DecodeFOV"
    assert dec.calls[-1][1] == (7,)
    dec2 = FakeDecoder()
    process_fov("/tmp/fake", array_id=5, decoder_loader=lambda path: dec2)
    assert [c[0] for c in dec2.calls] == ["WarpFOV", "PreprocessFOV"]


def test_combine_found_features_sequence():
    dec = FakeDecoder()
    combine_found_features("/tmp/fake", overwrite=True, decoder_loader=lambda path: dec)
    assert dec.overwrite is True
    assert [c[0] for c in dec.calls] == ["CombineFeatures", "GenerateFoundFeaturesReport", "FoundFeaturesToCSV"]
    assert dec.calls[-1][2] == {"downSampleFactor": 10, "zIndex": 4}


def test_optimize_sequence():
    dec = FakeDecoder()
    optimize("/tmp/fake", overwrite=True, decoder_loader=lambda path: dec, workers=3)
    assert [c[0] for c in dec.calls] == [
        "GenerateWarpReport",
        "InitializeScaleFactors",
        "SetParallel",
        "OptimizeScaleFactors",
        "Save",
    ]
    assert dec.calls[3][1] == (25,)
    assert dec.calls[3][2] == {"overwrite": True, "useBlanks": False}


def test_calculate_doublet_score_outputs_csvs():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        reports = root / "reports"
        reports.mkdir()
        pd.DataFrame(
            {
                "in_feature": [1, 1, 1, 0],
                "feature_id": [1, 1, 2, 2],
                "barcode_id": [1, 2, 1, 3],
                "abs_position_1": [10.0, 20.0, 30.0, 99.0],
                "abs_position_2": [5.0, 6.0, 7.0, 99.0],
            }
        ).to_csv(reports / "barcode_metadata.csv", index=False)
        dec = FakeDecoder()
        outputs = calculate_doublet_score(str(root), decoder_loader=lambda path: dec)
        assert all(path.exists() for path in outputs.values())
        rows = list(csv.reader(open(outputs["counts"], newline="")))
        assert len(rows) == 2
        assert float(rows[0][0]) == 1.0
        assert float(rows[0][1]) == 1.0
        assert float(rows[1][0]) == 1.0


def test_scheduler_writes_expected_scripts():
    with tempfile.TemporaryDirectory() as tmp:
        dec = FakeDecoder()
        result = merfish_scheduler(
            nDPath=tmp,
            aControl="wdsi",
            skip_decoder_construction=True,
            decoder_loader=lambda path: dec,
            dry_run=True,
            load_python_command="",
            activate_environment_command="",
        )
        names = sorted(path.name for path in result.scripts)
        assert "w_f1.slurm" in names
        assert "d_f1.slurm" in names
        assert "s_f1.slurm" in names
        assert "combine_sum.slurm" in names
        text = (Path(tmp) / "scheduler" / "decoding" / "d_f1.slurm").read_text()
        assert "analysis.SLURM_scripts.DecodeFOV" in text
