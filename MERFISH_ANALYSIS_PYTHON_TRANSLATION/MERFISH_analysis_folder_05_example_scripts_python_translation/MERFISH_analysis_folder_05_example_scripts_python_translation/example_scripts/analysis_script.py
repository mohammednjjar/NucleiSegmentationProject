"""Python translation of `example_scripts/analysis_script.m`.

Purpose: illustrate an end-to-end MERFISH analysis run on the example data,
then generate summary reports and save decoded outputs.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import shutil

try:
    from codes.CodebookToMap import CodebookToMap
    from codes.SECDEDCorrectableWords import SECDEDCorrectableWords
except Exception:
    CodebookToMap = None
    SECDEDCorrectableWords = None

try:
    from .script_utils import import_symbol, load_bytestream, page_break, save_bytestream, save_mat_file, set_figure_save_path
except ImportError:
    from script_utils import import_symbol, load_bytestream, page_break, save_bytestream, save_mat_file, set_figure_save_path


@dataclass
class AnalysisParameters:
    imageTag: str = "STORM"
    imageMListType: str = "alist"
    fiducialMListType: str = "list"
    numHybs: int = 16
    bitOrder: list[int] = field(default_factory=lambda: list(range(16, 0, -1)))
    wordConstMethod: str = "perLocalization"
    codebookPath: str = ""
    exactMap: Any = None
    errCorrFunc: Any = None
    FPKMData: Any = None
    cellsToAnalyze: list[Any] = field(default_factory=list)
    savePath: str = ""
    reportsToGenerate: list[tuple[str, str]] = field(default_factory=list)
    overwrite: bool = True
    figFormats: list[str] = field(default_factory=lambda: ["fig", "png"])
    useSubFolderForCellReport: bool = True
    saveAndClose: bool = True

    def as_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


def build_parameters(example_data_path: str | Path, analysis_base_path: str | Path) -> AnalysisParameters:
    example_data_path = Path(example_data_path)
    params = AnalysisParameters()
    params.codebookPath = str(example_data_path / "codebook" / "codebook.fasta")
    if CodebookToMap is not None:
        params.exactMap = CodebookToMap(params.codebookPath, keyType="binStr")
    else:
        mapper = import_symbol(["codes.CodebookToMap.CodebookToMap"])
        params.exactMap = mapper(params.codebookPath, keyType="binStr")
    params.errCorrFunc = SECDEDCorrectableWords
    params.FPKMData = load_bytestream(example_data_path / "FPKM_data" / "FPKMData.matb")
    params.savePath = set_figure_save_path(analysis_base_path, make_dir=True)
    params.reportsToGenerate.extend(
        [
            ("fiducialReport2", "off"),
            ("numOnBitsHistByCell", "off"),
            ("focusLockReportByCell", "off"),
            ("totalFPKMReport", "off"),
            ("cellByCellFPKMReport", "off"),
            ("cellWithWordsImage", "off"),
            ("molStats", "off"),
            ("molDistStats", "off"),
            ("compositeHybImage", "off"),
            ("hamming1DReportAllGenes", "off"),
            ("bitFlipProbabilitiesAverage", "off"),
            ("bitFlipProbabilitiesAllGenes", "off"),
            ("confidenceRatioReport", "off"),
        ]
    )
    return params


def run_analysis(example_data_path: str | Path, analysis_base_path: str | Path) -> dict[str, Any]:
    example_data_path = Path(example_data_path)
    parameters = build_parameters(example_data_path, analysis_base_path)

    analyze_merfish = import_symbol([
        "deprecated.AnalyzeMERFISH.AnalyzeMERFISH",
        "analysis.AnalyzeMERFISH.AnalyzeMERFISH",
        "AnalyzeMERFISH.AnalyzeMERFISH",
    ])
    generate_fpkm_report = import_symbol([
        "deprecated.reports.GenerateFPKMReport.GenerateFPKMReport",
        "GenerateFPKMReport.GenerateFPKMReport",
    ])
    generate_on_bit_histograms = import_symbol([
        "deprecated.reports.GenerateOnBitHistograms.GenerateOnBitHistograms",
        "GenerateOnBitHistograms.GenerateOnBitHistograms",
    ])
    generate_molecule_stats_report = import_symbol([
        "deprecated.reports.GenerateMoleculeStatsReport.GenerateMoleculeStatsReport",
        "GenerateMoleculeStatsReport.GenerateMoleculeStatsReport",
    ])
    generate_bit_flip_report = import_symbol([
        "deprecated.reports.GenerateBitFlipReport.GenerateBitFlipReport",
        "GenerateBitFlipReport.GenerateBitFlipReport",
    ])
    generate_hamming_sphere_report = import_symbol([
        "deprecated.reports.GenerateHammingSphereReport.GenerateHammingSphereReport",
        "GenerateHammingSphereReport.GenerateHammingSphereReport",
    ])
    strip_words = import_symbol([
        "deprecated.StripWords.StripWords",
        "analysis.deprecated.StripWords.StripWords",
        "StripWords.StripWords",
    ])

    words, image_data, fiducial_data, parameters_out = analyze_merfish(
        str(example_data_path / "example_data"), parameters=parameters.as_dict()
    )

    generate_fpkm_report(
        words,
        parameters.FPKMData,
        parameters=parameters.as_dict(),
        reportsToGenerate=("totalFPKMReport", "off"),
        FPKMReportExactMatchOnly=False,
        showNames=False,
    )
    generate_on_bit_histograms(words, parameters=parameters.as_dict(), reportsToGenerate=("numOnBitsHistAllCells", "off"))
    molecule_stats = generate_molecule_stats_report(words, parameters=parameters.as_dict())
    bit_flip_report = generate_bit_flip_report(words, parameters.exactMap, parameters=parameters.as_dict())
    generate_hamming_sphere_report(words, parameters.exactMap, parameters=parameters.as_dict())

    save_path = Path(parameters.savePath)
    save_bytestream(save_path / "imageData.matb", image_data)
    save_bytestream(save_path / "fiducialData.matb", fiducial_data)
    try:
        save_mat_file(save_path / "parameters.mat", {"parameters": parameters_out})
    except Exception:
        print("Warning: corrupt parameters file or incompatible MATLAB structure")
    save_bytestream(save_path / "bitFlipReport.matb", bit_flip_report)
    save_bytestream(save_path / "moleculeStats.matb", molecule_stats)

    word_breaks = sorted(set([1] + list(range(300000, len(words) + 1, 300000)) + [len(words)]))
    for j in range(1, len(word_breaks)):
        start = word_breaks[j - 1]
        stop = word_breaks[j]
        save_bytestream(save_path / f"words{j}.matb", words[start:stop])

    stripped_words = strip_words(words)
    save_bytestream(save_path / "strippedWords.matb", stripped_words)
    page_break()
    print(f"Saved words, imageData, parameters, and reports to {parameters.savePath}")
    return {
        "words": words,
        "imageData": image_data,
        "fiducialData": fiducial_data,
        "parameters": parameters_out,
        "moleculeStats": molecule_stats,
        "bitFlipReport": bit_flip_report,
    }


def main(merfish_analysis_path: str | Path = ".") -> dict[str, Any]:
    merfish_analysis_path = Path(merfish_analysis_path)
    example_data_path = merfish_analysis_path / "MERFISH_Examples"
    analysis_base_path = example_data_path / "MERFISH_Demo_Output"
    return run_analysis(example_data_path, analysis_base_path)


if __name__ == "__main__":
    main()
