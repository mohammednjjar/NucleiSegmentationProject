"""Python translation of `deprecated/reports/GenerateMoleculeStatsReport.m`."""
from __future__ import annotations

import numpy as np

try:
    from ..helpers import histogram
except ImportError:
    from deprecated.helpers import histogram


def _field(obj, name):
    return obj.get(name) if isinstance(obj, dict) else getattr(obj, name)


def _iqr(arr):
    return float(np.nanpercentile(arr, 75) - np.nanpercentile(arr, 25)) if len(arr) else float('nan')


def GenerateMoleculeStatsReport(words, reportsToGenerate=None, printedUpdates: bool = True,
                                saveAndClose: bool = False, useSubFolderForCellReport: bool = True,
                                overwrite: bool = True, figFormats=('png', 'fig'), brightHistBins: int = 100,
                                distHistBins: int = 100):
    if not words:
        raise ValueError('words is empty')
    num_hyb = int(_field(words[0], 'numHyb'))
    stats = {}
    for fld in ['a', 'bg', 'h']:
        data = np.asarray([v for w in words for v in np.ravel(_field(w, fld))], dtype=float)
        hist_arr = np.zeros((num_hyb, 2, int(brightHistBins)), dtype=float)
        for i in range(num_hyb):
            local = data[i::num_hyb]
            local = local[np.isfinite(local)]
            counts, centers = histogram(local, int(brightHistBins))
            hist_arr[i, 0, :] = centers
            hist_arr[i, 1, :] = counts
            stats[f'{fld}Mean'] = stats.get(f'{fld}Mean', []) + [float(np.nanmean(local)) if local.size else float('nan')]
            stats[f'{fld}STD'] = stats.get(f'{fld}STD', []) + [float(np.nanstd(local)) if local.size else float('nan')]
            stats[f'{fld}SEM'] = stats.get(f'{fld}SEM', []) + [float(np.nanstd(local) / np.sqrt(local.size)) if local.size else float('nan')]
            stats[f'{fld}Median'] = stats.get(f'{fld}Median', []) + [float(np.nanmedian(local)) if local.size else float('nan')]
            stats[f'{fld}IQR'] = stats.get(f'{fld}IQR', []) + [_iqr(local)]
            stats[f'{fld}N'] = stats.get(f'{fld}N', []) + [int(local.size)]
        stats[f'{fld}Hist'] = hist_arr
    multi = [w for w in words if _field(w, 'numOnBits') > 1]
    if multi:
        x_pos = np.asarray([v for w in multi for v in np.ravel(_field(w, 'xc'))], dtype=float)
        y_pos = np.asarray([v for w in multi for v in np.ravel(_field(w, 'yc'))], dtype=float)
        x_pos[x_pos == 0] = np.nan
        y_pos[y_pos == 0] = np.nan
        x_cent = np.asarray([_field(w, 'wordCentroidX') for w in multi], dtype=float)
        y_cent = np.asarray([_field(w, 'wordCentroidY') for w in multi], dtype=float)
        x_dist = x_pos.reshape(len(multi), num_hyb).T - x_cent[None, :]
        y_dist = y_pos.reshape(len(multi), num_hyb).T - y_cent[None, :]
    else:
        x_dist = np.empty((num_hyb, 0))
        y_dist = np.empty((num_hyb, 0))
    stats['xDist'] = x_dist
    stats['yDist'] = y_dist
    for name, dist in [('x', x_dist), ('y', y_dist)]:
        hist_arr = np.zeros((num_hyb, 2, int(distHistBins)), dtype=float)
        for i in range(num_hyb):
            data = dist[i, :]
            data = data[np.isfinite(data)]
            counts, centers = histogram(data, int(distHistBins))
            hist_arr[i, 0, :] = centers
            hist_arr[i, 1, :] = counts
            stats[f'{name}Mean'] = stats.get(f'{name}Mean', []) + [float(np.nanmean(data)) if data.size else float('nan')]
            stats[f'{name}STD'] = stats.get(f'{name}STD', []) + [float(np.nanstd(data)) if data.size else float('nan')]
            stats[f'{name}MeanAbs'] = stats.get(f'{name}MeanAbs', []) + [float(np.nanmean(np.abs(data))) if data.size else float('nan')]
        stats[f'{name}Hist'] = hist_arr
    stats['figHandles'] = []
    if reportsToGenerate:
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.bar(np.arange(1, num_hyb + 1), stats['aN'])
            ax.set_xlabel('Hybe Number')
            ax.set_ylabel('Counts')
            stats['figHandles'].append(fig)
        except Exception:
            stats['figHandles'] = []
    return stats, {'brightHistBins': brightHistBins, 'distHistBins': distHistBins}


def generate_molecule_stats_report(*args, **kwargs):
    return GenerateMoleculeStatsReport(*args, **kwargs)
