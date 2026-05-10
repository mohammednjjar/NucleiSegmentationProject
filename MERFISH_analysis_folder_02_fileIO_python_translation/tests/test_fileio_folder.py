from __future__ import annotations

import csv
import tempfile
from pathlib import Path

import numpy as np

from fileIO import (
    BuildFileStructure,
    CatBinaryFiles,
    LoadByteStream,
    LoadCodebook,
    LoadDataOrganization,
    LoadSplitByteStream,
    ReadBinaryFile,
    ReadBinaryFileHeader,
    SaveAsByteStream,
    SaveSplitByteStream,
    WriteBinaryFile,
    WriteCodebook,
)
from fileIO.deprecated import BuildImageDataStructures, CreateImageDataStructure, CreateWordsStructure


def test_build_file_structure():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "STORM_01_03.dax").write_text("x")
        out, params = BuildFileStructure(
            d,
            fileExt="dax",
            delimiters=["_"],
            fieldNames=["movieType", "hybNum", "cellNum"],
            fieldConv=[str, int, int],
            requireExactMatch=True,
        )
        assert len(out) == 1
        assert out[0]["movieType"] == "STORM"
        assert out[0]["hybNum"] == 1
        assert out[0]["cellNum"] == 3


def test_codebook_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "codebook.csv"
        readouts = [{"Header": "bit1"}, {"Header": "bit2"}, {"Header": "bit3"}]
        WriteCodebook(str(path), np.array([[1, 0, 1], [0, 1, 1]]), readouts, ["GeneA", "GeneB"], ["idA", "idB"], overwrite=True, verbose=False)
        cb, header, _ = LoadCodebook(str(path), verbose=False)
        assert header["version"] == "1.0"
        assert cb[0]["name"] == "GeneA"
        assert cb[0]["barcode"] == "101"


def test_binary_roundtrip_and_cat():
    with tempfile.TemporaryDirectory() as d:
        p1 = Path(d) / "a.bin"
        p2 = Path(d) / "b.bin"
        records1 = [
            {"x": np.array([1.0, 2.0], dtype=np.float64), "id": np.array([7], dtype=np.uint32)},
            {"x": np.array([3.0, 4.0], dtype=np.float64), "id": np.array([8], dtype=np.uint32)},
        ]
        records2 = [{"x": np.array([5.0, 6.0], dtype=np.float64), "id": np.array([9], dtype=np.uint32)}]
        WriteBinaryFile(str(p1), records1, verbose=False)
        WriteBinaryFile(str(p2), records2, verbose=False)
        header, layout = ReadBinaryFileHeader(str(p1))
        assert header["numEntries"] == 2
        assert layout[0]["field"] == "x"
        loaded, header, _ = ReadBinaryFile(str(p1))
        assert len(loaded) == 2
        assert np.allclose(loaded[1]["x"], [3.0, 4.0])
        CatBinaryFiles(str(p1), str(p2))
        loaded2, header2, _ = ReadBinaryFile(str(p1))
        assert header2["numEntries"] == 3
        assert int(loaded2[2]["id"]) == 9


def test_bytestream_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "obj.matb"
        obj = {"a": [1, 2, 3], "b": "text"}
        SaveAsByteStream(str(path), obj, verbose=False)
        loaded, _ = LoadByteStream(str(path), verbose=False)
        assert loaded == obj


def test_split_bytestream_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        arr = np.arange(12).reshape(3, 4, order="F")
        header_path = Path(d) / "split.matb"
        SaveSplitByteStream(str(header_path), arr, blockSize=5, verbose=False)
        loaded, info, _ = LoadSplitByteStream(str(header_path), verbose=False)
        assert info["numBlocks"] == 3
        assert np.array_equal(loaded, arr)


def test_data_organization():
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "data.csv"
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["channel", "frame", "zPos"])
            writer.writeheader()
            writer.writerow({"channel": "DAPI", "frame": "1;2", "zPos": "0.0;1.0"})
        data, meta, _ = LoadDataOrganization(str(path))
        assert data[0]["frame"] == [1, 2]
        assert meta["numZPos"] == 2


def test_deprecated_structures():
    imgs = CreateImageDataStructure(2)
    words = CreateWordsStructure(1, 3)
    assert len(imgs) == 2
    assert imgs[0]["tform"].shape == (3, 3)
    assert len(words) == 1
    assert words[0]["codeword"].shape == (3,)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "movie_1.dax").write_text("x")
        image_data, _ = BuildImageDataStructures(d, fileExt="dax", delimiters=["_"], fieldNames=["movieType", "hybNum"], fieldConv=[str, int], requireExactMatch=True)
        assert image_data[0]["movieType"] == "movie"
        assert "mList" in image_data[0]
