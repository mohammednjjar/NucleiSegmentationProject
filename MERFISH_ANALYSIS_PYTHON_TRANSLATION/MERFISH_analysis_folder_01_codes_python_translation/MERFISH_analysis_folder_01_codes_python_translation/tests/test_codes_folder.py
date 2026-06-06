from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from codes import (
    CodebookToMap,
    FastaRecord,
    GenSECDED,
    GenerateExtendedHammingWords,
    GenerateSurroundingCodewords,
    SECDEDCorrectableWords,
)


def test_generate_surrounding():
    out = GenerateSurroundingCodewords(np.array([True, False, True]), 1)
    assert len(out) == 3
    assert all(x.dtype == bool for x in out)
    mat = GenerateSurroundingCodewords(np.array([True, False, True]), 2, logical=True)
    assert mat.shape == (3, 3)
    assert (mat.sum(axis=1) == np.array([2, 0, 2])).all()


def test_secded_correctable():
    out = SECDEDCorrectableWords(np.array([True, False, True, False]))
    assert len(out) == 4


def test_codebook_to_map_records():
    records = [
        FastaRecord("1010", "GeneA extra info"),
        FastaRecord("0101", "GeneB extra info"),
    ]
    m, genes, words, params = CodebookToMap(records)
    assert m[int("1010", 2)] == "GeneA"
    assert genes == ["GeneA", "GeneB"]
    assert words == [10, 5]


def test_codebook_to_map_fasta_and_error_correction():
    fasta = ">101\nGeneA info\n>010\nGeneB info\n"
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        f.write(fasta)
        path = f.name
    try:
        m, genes, words, params = CodebookToMap(path, errCorrFunc=SECDEDCorrectableWords, keyType="binStr")
        assert m["101"] == "GeneA"
        assert m["001"] in {"GeneA", "GeneB"}
    finally:
        Path(path).unlink(missing_ok=True)


def test_extended_hamming_words():
    words, gen, r = GenerateExtendedHammingWords(4)
    assert words.shape[0] == 16
    assert gen.shape[0] == 4
    assert words.shape[1] == gen.shape[1]
    filtered, _, _ = GenerateExtendedHammingWords(4, numOn=4)
    assert np.all(filtered.sum(axis=1) == 4)


def test_gen_secded():
    words = GenSECDED(16, 11, 4)
    assert words.shape[1] == 16
    assert np.all(words.sum(axis=1) == 4)
    special = GenSECDED(12, 11, 4)
    assert special.shape[1] == 12
    assert np.all(special.sum(axis=1) == 4)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
    print("All tests passed.")
