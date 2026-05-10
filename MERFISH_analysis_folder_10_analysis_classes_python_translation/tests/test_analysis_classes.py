from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
import tifffile

from analysis.classes import FoundFeature, MERFISHDecoder, SLURMJob, SLURMJobArray, Vrykolakas


def test_found_feature_geometry():
    labels = np.zeros((30, 30), dtype=int)
    labels[8:20, 9:22] = 1
    f = FoundFeature(labels, fov_id=1, fov_center_pos=[0, 0], pixel_size=1000, stage_orientation=[1, 1], bounding_box=[0, 0, 30, 30], z_pos=[0], feature_label=1)
    assert f.abs_volume > 0
    c = f.CalculateCentroid()
    assert c.shape == (3,)
    assert f.IsInFeature([0, 0, 0]) in (True, False)
    table = f.Feature2Table()
    assert "feature_uID" in table.columns


def test_decode_pixels():
    img = np.zeros((4, 4, 2), dtype=float)
    img[:2, :, 0] = 10
    img[2:, :, 1] = 10
    vectors = np.eye(2)
    decoded, mag, traces, dist = MERFISHDecoder.DecodePixels(img, [1, 1], vectors, 0.2)
    assert decoded.shape == (4, 4)
    assert decoded[:2].max() == 1
    assert decoded[2:].max() == 2


def test_decoder_file_flow(tmp_path: Path):
    raw = tmp_path / "raw"
    norm = tmp_path / "norm"
    raw.mkdir()
    arr = np.zeros((8, 8, 2), dtype=np.float32)
    arr[:, :, 0] = 1
    tifffile.imwrite(raw / "fov1.tif", arr)
    md = MERFISHDecoder(raw, norm, verbose=False)
    md.numBits = 2
    md.scaleFactors = np.ones(2)
    md.codebook = pd.DataFrame({"bit1": [1, 0], "bit2": [0, 1]})
    out = md.DecodeFOV([1])[0]
    assert isinstance(out, pd.DataFrame)
    assert (norm / "barcodes" / "decoded_fov_1.tif").exists()


def test_slurm_job_dry_run(tmp_path: Path):
    job = SLURMJob(["echo hello"], jobName="x", scriptPath=str(tmp_path), dryRun=True)
    assert job.Submit() is True
    assert job.script_file.exists()
    arr = SLURMJobArray([job], name="arr")
    short, long = arr.Status()
    assert "arr" in short
