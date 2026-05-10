"""Python translation of `example_scripts/code_construction_script.m`.

Purpose: demonstrate construction of an MHD4 extended-Hamming MERFISH code and
an MHD2 constant-weight code.
"""
from __future__ import annotations

import numpy as np

try:
    from codes.GenerateExtendedHammingWords import GenerateExtendedHammingWords
except Exception:
    GenerateExtendedHammingWords = None

try:
    from .script_utils import combinations_to_binary, pairwise_hamming_matrix, page_break
except ImportError:
    from script_utils import combinations_to_binary, pairwise_hamming_matrix, page_break


def _generate_extended_hamming_words(num_data_bits: int) -> np.ndarray:
    if GenerateExtendedHammingWords is not None:
        words = GenerateExtendedHammingWords(num_data_bits)
        if isinstance(words, tuple):
            return np.asarray(words[0], dtype=np.uint8)
        return np.asarray(words, dtype=np.uint8)
    try:
        from MERFISH_analysis_folder_01_codes_python_translation.codes.GenerateExtendedHammingWords import GenerateExtendedHammingWords as local_func
    except Exception as exc:
        raise RuntimeError("GenerateExtendedHammingWords from folder 01 is required") from exc
    out = local_func(num_data_bits)
    return np.asarray(out[0] if isinstance(out, tuple) else out, dtype=np.uint8)


def construct_mhd4_code(num_data_bits: int = 11, hamming_weight: int = 4) -> dict[str, object]:
    extended_hamming_words = _generate_extended_hamming_words(num_data_bits)
    weights = extended_hamming_words.sum(axis=1)
    mhd4_words = extended_hamming_words[weights == int(hamming_weight), :]
    measured_distances = pairwise_hamming_matrix(mhd4_words)
    minimum_distance_by_word = measured_distances.min(axis=1)
    return {
        "words": mhd4_words,
        "unique_weights": np.unique(mhd4_words.sum(axis=1)),
        "minimum_distances": np.unique(minimum_distance_by_word),
    }


def construct_mhd2_code(num_bits: int = 14, on_bits: int = 4) -> dict[str, object]:
    mhd2_words = combinations_to_binary(num_bits, on_bits)
    measured_distances = pairwise_hamming_matrix(mhd2_words)
    minimum_distance_by_word = measured_distances.min(axis=1)
    return {
        "words": mhd2_words,
        "unique_weights": np.unique(mhd2_words.sum(axis=1)),
        "minimum_distances": np.unique(minimum_distance_by_word),
    }


def main() -> dict[str, dict[str, object]]:
    mhd4 = construct_mhd4_code()
    page_break()
    print(f"Constructed {mhd4['words'].shape[0]} barcodes/words")
    print("Found the following hamming weights")
    print(mhd4["unique_weights"])
    print("Found the following minimum HD")
    print(mhd4["minimum_distances"])

    mhd2 = construct_mhd2_code()
    page_break()
    print(f"Constructed {mhd2['words'].shape[0]} barcodes/words")
    print("Found the following hamming weights")
    print(mhd2["unique_weights"])
    print("Found the following minimum HD")
    print(mhd2["minimum_distances"])
    return {"MHD4": mhd4, "MHD2": mhd2}


if __name__ == "__main__":
    main()
