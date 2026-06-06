"""Python translation of `deprecated/reports/GenerateCompositeImage.m`."""
from __future__ import annotations

import numpy as np

try:
    from ..helpers import field_vector
except ImportError:
    from deprecated.helpers import field_vector


def GenerateCompositeImage(words, imageData, reportsToGenerate=None, printedUpdates: bool = True,
                           showNames: bool = False, embedNames: bool = True, saveAndClose: bool = False,
                           useSubFolderForCellReport: bool = True, overwrite: bool = True,
                           figFormats=('png', 'fig'), numImageColumns: int = 4, displayHybLabel: bool = True):
    figs = []
    cells = sorted(set(field_vector(imageData, 'cellNum')))
    try:
        import matplotlib.pyplot as plt
        for cell in cells:
            local_words = [w for w in words if (w.get('cellID') if isinstance(w, dict) else getattr(w, 'cellID')) == cell]
            fig, ax = plt.subplots()
            xs_exact = [w.get('wordCentroidX') for w in local_words if w.get('isExactMatch')]
            ys_exact = [w.get('wordCentroidY') for w in local_words if w.get('isExactMatch')]
            xs_corr = [w.get('wordCentroidX') for w in local_words if w.get('isCorrectedMatch')]
            ys_corr = [w.get('wordCentroidY') for w in local_words if w.get('isCorrectedMatch')]
            xs_other = [w.get('wordCentroidX') for w in local_words if not (w.get('isExactMatch') or w.get('isCorrectedMatch'))]
            ys_other = [w.get('wordCentroidY') for w in local_words if not (w.get('isExactMatch') or w.get('isCorrectedMatch'))]
            if xs_exact:
                ax.plot(xs_exact, ys_exact, 'go')
            if xs_corr:
                ax.plot(xs_corr, ys_corr, 'rx')
            if xs_other:
                ax.plot(xs_other, ys_other, 'b.')
            ax.set_title(f'CellWithWords_Cell_{cell}')
            figs.append(fig)
    except Exception:
        figs = []
    return figs, {'numImageColumns': numImageColumns, 'displayHybLabel': displayHybLabel}


def generate_composite_image(*args, **kwargs):
    return GenerateCompositeImage(*args, **kwargs)
