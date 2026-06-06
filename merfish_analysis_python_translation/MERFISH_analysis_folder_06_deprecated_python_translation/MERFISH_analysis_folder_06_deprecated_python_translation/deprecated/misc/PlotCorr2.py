"""Python translation of `deprecated/misc/PlotCorr2.m`."""
from __future__ import annotations

from typing import Sequence
import numpy as np

try:
    from ..helpers import pearson_corr
except ImportError:
    from deprecated.helpers import pearson_corr


def PlotCorr2(x: Sequence[float], y: Sequence[float], MarkerSize: float = 10, FontSize: float = 6,
              colorMap: str = 'jet', nameBuffer: float = 0.1, figHandle=None, axesHandle=None,
              pointNames: Sequence[str] | None = None, plotFunction: str = 'loglog',
              includeLog10: bool = True, includeLin: bool = True):
    xarr = np.asarray(x, dtype=float).reshape(-1)
    yarr = np.asarray(y, dtype=float).reshape(-1)
    if xarr.shape != yarr.shape:
        raise ValueError('x and y must have the same length')
    finite = np.isfinite(xarr) & np.isfinite(yarr)
    positive = (xarr > 0) & (yarr > 0) & finite
    out = {'log10rho': None, 'log10pvalue': None, 'rho': None, 'pvalue': None, 'figure': None, 'axes': None}
    if positive.sum() == 0:
        return out
    rho, pvalue = pearson_corr(xarr[finite], yarr[finite])
    out['rho'] = rho
    out['pvalue'] = pvalue
    if includeLog10:
        lrho, lpvalue = pearson_corr(np.log10(xarr[positive]), np.log10(yarr[positive]))
        out['log10rho'] = lrho
        out['log10pvalue'] = lpvalue
    try:
        import matplotlib.pyplot as plt
        fig = figHandle if figHandle is not None else plt.figure()
        ax = axesHandle if axesHandle is not None else fig.add_subplot(111)
        if plotFunction == 'plot':
            ax.plot(xarr, yarr, 'k.', markersize=MarkerSize)
        else:
            ax.loglog(xarr[positive], yarr[positive], 'k.', markersize=MarkerSize)
        if pointNames:
            for xi, yi, name in zip(xarr, yarr, pointNames):
                if np.isfinite(xi) and np.isfinite(yi):
                    ax.text(xi + nameBuffer * xi, yi, str(name), fontsize=FontSize)
        title_parts = []
        if includeLog10:
            title_parts.append(f'rho_log10 = {out["log10rho"]:.2g} (p={out["log10pvalue"]:.2g})')
        if includeLin:
            title_parts.append(f'rho = {out["rho"]:.2g} (p={out["pvalue"]:.2g})')
        ax.set_title('  '.join(title_parts))
        out['figure'] = fig
        out['axes'] = ax
    except Exception:
        out['figure'] = None
        out['axes'] = None
    return out


def plot_corr2(*args, **kwargs):
    return PlotCorr2(*args, **kwargs)
