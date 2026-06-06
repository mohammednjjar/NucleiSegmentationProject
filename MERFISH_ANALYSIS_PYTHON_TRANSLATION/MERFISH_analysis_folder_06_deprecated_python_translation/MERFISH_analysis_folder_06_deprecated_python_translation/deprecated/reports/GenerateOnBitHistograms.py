"""Python translation of `deprecated/reports/GenerateOnBitHistograms.m`."""
from __future__ import annotations


def _field(obj, name):
    return obj.get(name) if isinstance(obj, dict) else getattr(obj, name)


def GenerateOnBitHistograms(words, reportsToGenerate=(('numOnBitsHistByCell', 'on'), ('numOnBitsHistAllCells', 'on')),
                            printedUpdates: bool = True, saveAndClose: bool = False,
                            useSubFolderForCellReport: bool = True, overwrite: bool = True,
                            figFormats=('png', 'fig')):
    figs = []
    params = {'reportsToGenerate': reportsToGenerate, 'saveAndClose': saveAndClose}
    try:
        import matplotlib.pyplot as plt
        report_names = [r[0] if isinstance(r, (list, tuple)) else r for r in reportsToGenerate]
        num_hyb = int(_field(words[0], 'numHyb')) if words else 0
        if 'numOnBitsHistByCell' in report_names:
            for cell in sorted(set(_field(w, 'cellID') for w in words)):
                local = [_field(w, 'numOnBits') for w in words if _field(w, 'cellID') == cell]
                fig, ax = plt.subplots()
                ax.hist(local, bins=range(0, num_hyb + 2))
                ax.set_xlabel('Number of On Bits')
                ax.set_ylabel('Counts')
                ax.set_title(f'Cell {cell}')
                figs.append(fig)
        if 'numOnBitsHistAllCells' in report_names:
            fig, ax = plt.subplots()
            ax.hist([_field(w, 'numOnBits') for w in words], bins=range(0, num_hyb + 2))
            ax.set_xlabel('Number of On Bits')
            ax.set_ylabel('Counts')
            figs.append(fig)
    except Exception:
        figs = []
    return figs, params


def generate_on_bit_histograms(*args, **kwargs):
    return GenerateOnBitHistograms(*args, **kwargs)
