"""Python translation of `example_scripts/decoding/runMERFISH.m`.

Purpose: configure one MERFISH experiment and launch the translated MERFISH
scheduler with the same analysis-control string used by the MATLAB script.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import time

try:
    from ..script_utils import ensure_dir, import_symbol, page_break, tic, toc
except ImportError:
    from script_utils import ensure_dir, import_symbol, page_break, tic, toc


@dataclass
class MERFISHRunPaths:
    repositoryPath: str = "/n/regal/zhuang_lab/lsepulveda"
    experimentPath: str = "181029_BC037_MERFISH"
    samplePath: str = "sample_02"
    scratchPath: str = "/n/regal/zhuang_lab/lsepulveda/scratch/lsepulvedaduran"

    @property
    def dataPath(self) -> str:
        return str(Path(self.repositoryPath) / self.experimentPath)

    @property
    def settingsPath(self) -> str:
        return str(Path(self.dataPath) / "Settings" / self.samplePath)

    @property
    def rDPath(self) -> str:
        return str(Path(self.dataPath) / self.samplePath)

    @property
    def nDPath(self) -> str:
        return str(Path(self.dataPath) / "normalized_data" / f"{self.samplePath}-new-pip-test") + "/"

    @property
    def aPath(self) -> str:
        return "~/Software/code-matlab/merfish/fpkm/FPKMDataPublished.matb"

    @property
    def positionsPath(self) -> str:
        return str(Path(self.settingsPath) / "positions_full-1000_1024_60x_MERFISH.txt")

    @property
    def dataOrganizationPath(self) -> str:
        return str(Path(self.settingsPath) / "data_organization.csv")

    @property
    def codebookPath(self) -> str:
        return str(Path(self.settingsPath) / "L26E1_codebook.csv")


@dataclass
class DecoderParameters:
    pixelSize: float = 107.4
    imageSize: tuple[int, int] = (2048, 2048)
    lowPassKernelSize: int = 1
    crop: int = 40
    minBrightness: int = 1
    minArea: int = 1
    areaThresh: int = 4
    stageOrientation: tuple[int, int] = (-1, 1)
    overwrite: bool = False
    codebookPath: str = ""
    dataOrganizationPath: str = ""
    hal_version: str = "hal2"

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_decoder_parameters(paths: MERFISHRunPaths) -> DecoderParameters:
    return DecoderParameters(codebookPath=paths.codebookPath, dataOrganizationPath=paths.dataOrganizationPath)


def run_merfish(paths: MERFISHRunPaths | None = None, merfish_ENV: str = "merfish_analysis", aControl: str = "mwodscpfnlbu") -> dict[str, Any]:
    paths = paths or MERFISHRunPaths()
    log_dir = ensure_dir(Path(paths.nDPath) / "log")
    page_break()
    print(f"Creating MERFISHDecoder for {paths.rDPath}")
    print(f"Normalized data will be saved in {paths.nDPath}")
    print(f"Started at {time.ctime()}")
    script_timer = tic()

    parameters = build_decoder_parameters(paths)
    scheduler = import_symbol([
        "analysis.SLURM_scripts.MERFISHScheduler.MERFISHScheduler",
        "MERFISHScheduler.MERFISHScheduler",
    ])
    result = scheduler(
        rDPath=paths.rDPath,
        nDPath=paths.nDPath,
        aPath=paths.aPath,
        positionsPath=paths.positionsPath,
        dataOrganizationPath=paths.dataOrganizationPath,
        codebookPath=paths.codebookPath,
        scratchPath=paths.scratchPath,
        merfish_ENV=merfish_ENV,
        parameters=parameters.as_dict(),
        aControl=aControl,
    )
    page_break()
    print(f"...completed in {toc(script_timer)} s")
    print(f"Completed at {time.ctime()}")
    return {"paths": paths, "parameters": parameters, "scheduler_result": result, "log_dir": str(log_dir)}


def main() -> dict[str, Any]:
    return run_merfish()


if __name__ == "__main__":
    main()
