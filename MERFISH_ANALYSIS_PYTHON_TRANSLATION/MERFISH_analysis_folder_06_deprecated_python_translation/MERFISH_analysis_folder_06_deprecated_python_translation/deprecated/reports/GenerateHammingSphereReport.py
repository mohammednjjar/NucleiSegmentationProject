"""Python translation of `deprecated/reports/GenerateHammingSphereReport.m`."""
from __future__ import annotations

import re
import numpy as np

try:
    from ..helpers import bits_to_int, generate_surrounding_codewords
except ImportError:
    from deprecated.helpers import bits_to_int, generate_surrounding_codewords


def _field(obj, name):
    return obj.get(name) if isinstance(obj, dict) else getattr(obj, name)


def GenerateHammingSphereReport(words, exactMap, reportsToGenerate=None, printedUpdates: bool = True,
                                saveAndClose: bool = False, overwrite: bool = True, figFormats=('png', 'fig'),
                                subFolder: str = '', maxHammingSphere: int = 1,
                                blankWordIdentifiers=('blank', 'notarget'), colorMap='jet', numHistogramBins: int = 25):
    code_words = list(exactMap.keys())
    gene_names = list(exactMap.values())
    num_hybs = len(str(code_words[0]).strip())
    int_codewords = [bits_to_int(str(cw).strip()) for cw in code_words]
    word_ints = [_field(w, 'intCodeword') for w in words]
    counts = np.bincount(np.asarray(word_ints, dtype=int), minlength=2 ** num_hybs + 1)
    sphere_counts = np.zeros((maxHammingSphere + 1, len(gene_names)), dtype=float)
    sphere_counts[0, :] = counts[int_codewords]
    for dist in range(1, maxHammingSphere + 1):
        for j, cw in enumerate(code_words):
            surrounding = [bits_to_int(x) for x in generate_surrounding_codewords(str(cw).strip(), dist)]
            sphere_counts[dist, j] = counts[surrounding].sum() if surrounding else 0
    is_blank = np.zeros(len(gene_names), dtype=bool)
    for ident in blankWordIdentifiers:
        is_blank |= np.array([re.search(str(ident), str(g), flags=re.I) is not None for g in gene_names])
    denom = sphere_counts[[0, 1], :].sum(axis=0) if maxHammingSphere >= 1 else sphere_counts[0, :]
    ratio = np.divide(sphere_counts[0, :], denom, out=np.zeros_like(denom), where=denom != 0)
    hist_max = float(np.nanmax(ratio)) if ratio.size else 1.0
    bins = np.linspace(0, hist_max if hist_max > 0 else 1, int(numHistogramBins))
    n_non, edges = np.histogram(ratio[~is_blank], bins=bins)
    n_blank, _ = np.histogram(ratio[is_blank], bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2
    non_cdf = np.cumsum(n_non) / n_non.sum() if n_non.sum() else np.zeros_like(n_non, dtype=float)
    blank_cdf = np.cumsum(n_blank) / n_blank.sum() if n_blank.sum() else np.zeros_like(n_blank, dtype=float)
    sorted_ind = np.argsort(-ratio)
    report = {'geneNames': gene_names, 'hammingSphereCounts': sphere_counts,
              'blankNames': [g for g, b in zip(gene_names, is_blank) if b],
              'blankIDs': np.flatnonzero(is_blank).tolist(), 'nonBlankIDs': np.flatnonzero(~is_blank).tolist(),
              'hammingSphere01Ratio': ratio, 'nNonBlank': n_non, 'nBlank': n_blank,
              'histBins': centers, 'NonBlankCDF': non_cdf, 'blankCDF': blank_cdf,
              'sortedInd': sorted_ind.tolist(), 'figHandles': []}
    if reportsToGenerate:
        try:
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(2, 2)
            axes[0, 0].plot(ratio[sorted_ind])
            axes[1, 0].bar(centers, n_non, width=np.mean(np.diff(centers)) if len(centers) > 1 else 0.1)
            axes[1, 1].plot(centers, non_cdf)
            axes[1, 1].plot(centers, blank_cdf)
            report['figHandles'].append(fig)
        except Exception:
            report['figHandles'] = []
    return report, {'maxHammingSphere': maxHammingSphere}


def generate_hamming_sphere_report(*args, **kwargs):
    return GenerateHammingSphereReport(*args, **kwargs)
