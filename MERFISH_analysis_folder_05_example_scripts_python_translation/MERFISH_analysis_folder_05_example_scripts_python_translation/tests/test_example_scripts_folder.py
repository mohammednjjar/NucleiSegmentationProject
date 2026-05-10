from pathlib import Path
import importlib
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def test_code_construction_mhd2():
    mod = importlib.import_module("example_scripts.code_construction_script")
    out = mod.construct_mhd2_code(num_bits=5, on_bits=2)
    assert out["words"].shape == (10, 5)
    assert out["unique_weights"].tolist() == [2]
    assert out["minimum_distances"].tolist() == [2.0]


def test_script_utils_fasta_and_reverse_complement():
    utils = importlib.import_module("example_scripts.script_utils")
    records = [utils.FastaRecord("r1", "ACGT"), utils.FastaRecord("r2", "TTAA")]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "x.fa"
        utils.write_fasta(path, records)
        loaded = utils.read_fasta(path)
    assert [r.Header for r in loaded] == ["r1", "r2"]
    assert utils.reverse_complement("ACGTN") == "NACGT"
