
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "python_translation"
sys.path.insert(0, str(ROOT))

from analysis.deprecated.Warp2BestPair import warp2_best_pair
from analysis.deprecated.AlignFiducials import align_fiducials
from analysis.deprecated.TransformImageData import transform_image_data
from analysis.deprecated.CreateWords import create_words
from analysis.deprecated.DecodeWords import decode_words
from analysis.deprecated.StripWords import strip_words


def test_warp2_best_pair_aligns_translation():
    ref = np.array([[0.,0.],[10.,0.],[0.,10.],[10.,10.]])
    mov = ref + np.array([2., -3.])
    tform, errors = warp2_best_pair(ref, mov, maxD=20, showPlots=False)
    aligned = tform.inverse_points(mov)
    assert np.allclose(aligned, ref, atol=1e-6)
    assert np.nanmax(errors) < 1e-6


def test_align_and_transform_image_data():
    ref = {"mList": {"xc": [0,10,0,10], "yc": [0,0,10,10], "w": [80,80,80,80], "frame": [1,1,1,1]}, "cellNum": 1, "uID": "f1"}
    mov = {"mList": {"xc": [2,12,2,12], "yc": [-3,-3,7,7], "w": [80,80,80,80], "frame": [1,1,1,1]}, "cellNum": 1, "uID": "f2"}
    fids, params = align_fiducials([ref, mov], maxD=20, verbose=False, printedUpdates=False)
    img = [{"mList": {"x": [2], "y": [-3]}, "uID": "i1"}, {"mList": {"x": [12], "y": [7]}, "uID": "i2"}]
    out, _ = transform_image_data(img, fids, verbose=False, printedUpdates=False)
    assert np.isclose(out[1]["mList"]["xc"][0], 10.0)
    assert np.isclose(out[1]["mList"]["yc"][0], 10.0)
    assert out[1]["fidUID"] == "f2"


def test_create_decode_strip_words_per_localization():
    images = [
        {"name":"h1", "filePath":"h1.csv", "cellNum":1, "Stage_X":0, "Stage_Y":0, "uID":"u1", "fidUID":"f1", "hasFiducialError":False, "focusLockQuality":1, "mList":{"xc":[1, 20], "yc":[1, 20], "a":[5, 5]}},
        {"name":"h2", "filePath":"h2.csv", "cellNum":1, "Stage_X":0, "Stage_Y":0, "uID":"u2", "fidUID":"f2", "hasFiducialError":False, "focusLockQuality":1, "mList":{"xc":[1.2, 40], "yc":[1.1, 40], "a":[5, 5]}},
    ]
    words, params = create_words(images, numHybs=2, bitOrder=[1,2], maxDtoCentroid=0.5, verbose=False, printedUpdates=False)
    assert any(w["intCodeword"] == 3 for w in words)
    decoded, _ = decode_words(words, {"11":"GeneA"}, {}, keyType="binStr")
    assert any(w["geneName"] == "GeneA" and w["isExactMatch"] for w in decoded)
    stripped = strip_words(decoded, fieldsToKeep=["intCodeword", "geneName"])
    assert set(stripped[0].keys()).issubset({"intCodeword", "geneName"})


def test_create_words_common_centroid_empty_safe():
    images = [{"mList":{"xc": [], "yc": [], "a": []}, "cellNum": 1, "Stage_X": 0, "Stage_Y": 0}]
    words, _ = create_words(images, wordConstMethod="commonCentroid", numHybs=1, verbose=False, printedUpdates=False)
    assert words == []
