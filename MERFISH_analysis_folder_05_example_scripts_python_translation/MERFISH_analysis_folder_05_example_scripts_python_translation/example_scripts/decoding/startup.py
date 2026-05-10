"""Python translation of `example_scripts/decoding/startup.m`.

Purpose: configure local paths for STORM and MERFISH analysis dependencies used by
`runMERFISH.py`.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import os
import sys


@dataclass
class DecodingStartupConfig:
    homePath: str = "/n/home06/lsepulvedaduran"
    basePath: str = "/n/home06/lsepulvedaduran/Software/new_merfish_pipeline/"
    scratchPath: str = "/n/home06/lsepulvedaduran/scratch/"
    pythonPath: str = "/n/home06/lsepulvedaduran/.conda/envs/merfish_analysis/lib/python2.7/"
    matlabStormPath: str = "/n/home06/lsepulvedaduran/Software/new_merfish_pipeline/matlab-storm/"
    stormAnalysisPath: str = "/n/home06/lsepulvedaduran/Software/new_merfish_pipeline/storm-analysis/"
    MERFISHAnalysisPath: str = "/n/home06/lsepulvedaduran/Software/new_merfish_pipeline/MERFISH_analysis/"
    exportFigPath: str = "/n/home06/lsepulvedaduran/Software/new_merfish_pipeline/export_fig/"
    daoSTORMexe: str = ""

    def finalize(self) -> "DecodingStartupConfig":
        self.daoSTORMexe = (
            str(Path(self.homePath) / ".conda/envs/merfish_analysis/bin/python")
            + " "
            + str(Path(self.stormAnalysisPath) / "3d_daostorm/mufit_analysis.py")
        )
        return self


def merfish_decoding_startup(config: DecodingStartupConfig | None = None, add_to_sys_path: bool = True) -> DecodingStartupConfig:
    cfg = (config or DecodingStartupConfig()).finalize()
    print("-" * 66)
    print("Adding matlab-storm equivalent path")
    print(f"    {cfg.matlabStormPath}")
    print("Adding path to daoSTORM")
    print(f"    daoSTORMexe = {cfg.daoSTORMexe}")
    print("-" * 66)
    print("Adding MERFISH_analysis")
    print(f"    {cfg.MERFISHAnalysisPath}")
    print("Adding export_fig equivalent path")
    print(f"    {cfg.exportFigPath}")

    os.environ["MERFISH_SCRATCH_PATH"] = cfg.scratchPath
    os.environ["MERFISH_PYTHON_PATH"] = cfg.pythonPath
    os.environ["MERFISH_MATLAB_STORM_PATH"] = cfg.matlabStormPath
    os.environ["MERFISH_STORM_ANALYSIS_PATH"] = cfg.stormAnalysisPath
    os.environ["MERFISH_ANALYSIS_PATH"] = cfg.MERFISHAnalysisPath
    os.environ["MERFISH_DAOSTORM_EXE"] = cfg.daoSTORMexe

    if add_to_sys_path:
        for path in [cfg.MERFISHAnalysisPath, cfg.matlabStormPath, cfg.stormAnalysisPath, cfg.exportFigPath]:
            path_str = str(Path(path).expanduser())
            if path_str not in sys.path:
                sys.path.insert(0, path_str)
    return cfg


def main() -> dict[str, str]:
    return asdict(merfish_decoding_startup())


if __name__ == "__main__":
    main()
