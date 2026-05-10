"""Python translation of `deprecated/reports/GenerateBitFlipReport.m`."""
from __future__ import annotations

import numpy as np

try:
    from ..helpers import bits_to_int, int_to_bits, generate_surrounding_codewords
except ImportError:
    from deprecated.helpers import bits_to_int, int_to_bits, generate_surrounding_codewords


def GenerateBitFlipReport(words, exactMap, reportsToGenerate=None, printedUpdates: bool = True,
                          saveAndClose: bool = False, useSubFolderForCellReport: bool = True,
                          overwrite: bool = True, figFormats=('png', 'fig'), probToUse: str = 'exact',
                          errCorrFunc=None, numHybs: int = 16):
    if errCorrFunc is None:
        errCorrFunc = lambda cw: generate_surrounding_codewords(cw, 1)
    gene_names = list(exactMap.values())
    codewords = list(exactMap.keys())
    normalized = []
    for cw in codewords:
        if isinstance(cw, str):
            normalized.append(np.array([int(c) for c in cw.strip()], dtype=np.uint8))
        else:
            normalized.append(int_to_bits(int(cw), numHybs))
    int_words = np.array([int(w['intCodeword']) if isinstance(w, dict) else int(getattr(w, 'intCodeword')) for w in words], dtype=int)
    counts = np.bincount(int_words, minlength=2 ** numHybs + 1)
    first = np.full((len(gene_names), numHybs, 2), np.nan)
    exact = np.full((len(gene_names), numHybs, 2), np.nan)
    num_counts = np.zeros(len(gene_names), dtype=float)
    for i, codeword in enumerate(normalized):
        correct = bits_to_int(codeword)
        surrounds = [bits_to_int(x) for x in errCorrFunc(codeword)]
        if not surrounds:
            continue
        total = counts[[correct] + surrounds].sum()
        local_first = counts[surrounds] / total if total > 0 else np.full(len(surrounds), np.nan)
        alpha = counts[surrounds] / counts[correct] if counts[correct] > 0 else np.full(len(surrounds), np.nan)
        local_exact = alpha / (1 + alpha)
        num_counts[i] = total
        for j, word_int in enumerate(surrounds):
            diff = np.flatnonzero(codeword != int_to_bits(word_int, numHybs))
            if diff.size:
                bit = int(diff[0])
                transition = 0 if word_int < correct else 1
                first[i, bit, transition] = local_first[j]
                exact[i, bit, transition] = local_exact[j]
    probs = exact if probToUse == 'exact' else first
    with np.errstate(invalid='ignore', divide='ignore'):
        weight_vec = num_counts / num_counts.sum() if num_counts.sum() else np.zeros_like(num_counts)
        hyb_prob = np.nanmean(probs, axis=0)
        hyb_err = np.nanstd(probs, axis=0)
        scaled = np.nansum(probs * weight_vec[:, None, None], axis=0)
        scaled_err = np.sqrt(np.nanvar(probs, axis=0)) / max(np.sqrt(len(weight_vec)), 1)
    report = {'geneNames': gene_names, 'counts': num_counts, 'firstOrderProbabilities': first,
              'exactProbabilities': exact, 'probabilities': probs, 'hybProb': hyb_prob,
              'hybProbErr': hyb_err, 'numCounts': num_counts, 'scaledHybProb': scaled,
              'scaledHybProbErr': scaled_err, 'figHandles': []}
    if reportsToGenerate:
        try:
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 2)
            axes[0].bar(np.arange(1, numHybs + 1), scaled[:, 0])
            axes[1].bar(np.arange(1, numHybs + 1), scaled[:, 1])
            report['figHandles'].append(fig)
        except Exception:
            report['figHandles'] = []
    return report, {'probToUse': probToUse, 'numHybs': numHybs}


def generate_bit_flip_report(*args, **kwargs):
    return GenerateBitFlipReport(*args, **kwargs)
