"""Python translation of `deprecated/reports/GenerateFPKMReport.m`."""
from __future__ import annotations

try:
    from ..misc.PlotCorr2 import PlotCorr2
except ImportError:
    from deprecated.misc.PlotCorr2 import PlotCorr2


def _field(obj, name):
    return obj.get(name) if isinstance(obj, dict) else getattr(obj, name)


def _gene_counts(gene_names, target_names):
    remaining = list(gene_names)
    counts = []
    for target in target_names:
        c = sum(1 for g in remaining if g == target)
        counts.append(c)
        remaining = [g for g in remaining if g != target]
    counts.append(len(remaining))
    return counts


def GenerateFPKMReport(words, FPKMData, reportsToGenerate=None, printedUpdates: bool = True,
                       FPKMReportExactMatchOnly: bool = False, FPKMReportEmbedNames: bool = True,
                       showNames: bool = False, embedNames: bool = False, saveAndClose: bool = False,
                       useSubFolderForCellReport: bool = True, overwrite: bool = True, figFormats=('png', 'fig')):
    if FPKMReportExactMatchOnly:
        valid = [bool(_field(w, 'isExactMatch')) for w in words]
        modifier = 'ExactOnly'
    else:
        valid = [bool(_field(w, 'isExactMatch')) or bool(_field(w, 'isCorrectedMatch')) for w in words]
        modifier = 'ExactAndCorrected'
    valid_words = [w for w, ok in zip(words, valid) if ok]
    target_names = [_field(x, 'geneName') for x in FPKMData]
    counts = _gene_counts([_field(w, 'geneName') for w in valid_words], target_names)
    fpkm = [_field(x, 'FPKM') for x in FPKMData]
    corr = PlotCorr2(fpkm, counts[:-1], pointNames=target_names if showNames else None)
    report = {'wordInds': [i for i, ok in enumerate(valid) if ok], 'geneNames': target_names + ['unknown name'],
              'counts': counts, 'countsWOUnknown': counts[:-1], 'FPKM': fpkm, 'pearsonCorr': corr,
              'figNameModifier': modifier}
    if reportsToGenerate and any(r[0] == 'cellByCellFPKMReport' if isinstance(r, (list, tuple)) else r == 'cellByCellFPKMReport' for r in reportsToGenerate):
        cell_ids = sorted(set(_field(w, 'cellID') for w in words))
        report['cellReportCounts'] = []
        report['cellByCellPearsonCorr'] = []
        for cell in cell_ids:
            local = [w for w in valid_words if _field(w, 'cellID') == cell]
            local_counts = _gene_counts([_field(w, 'geneName') for w in local], target_names)
            report['cellReportCounts'].append(local_counts)
            report['cellByCellPearsonCorr'].append(PlotCorr2(fpkm, local_counts[:-1], pointNames=target_names if showNames else None))
    return report, {'FPKMReportExactMatchOnly': FPKMReportExactMatchOnly}


def generate_fpkm_report(*args, **kwargs):
    return GenerateFPKMReport(*args, **kwargs)
