from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

from startup import MerfishStartupConfig, add_path, configure_merfish_environment


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def test_add_path_recursive() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "MERFISH_analysis"
        nested = root / "analysis" / "classes"
        nested.mkdir(parents=True)
        inserted = add_path(str(root), recursive=True, prepend=True, strict=True)
        require(str(root) in inserted, "root directory was not inserted")
        require(str(nested) in inserted, "nested directory was not inserted")
        require(str(root) in sys.path, "root directory missing from sys.path")
        require(str(nested) in sys.path, "nested directory missing from sys.path")


def test_configure_environment() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        matlab_dir = base / "MATLAB"
        startup_dir = base / "startup"
        merfish_dir = base / "MERFISH_analysis"
        export_fig = base / "export_fig"
        storm_startup = base / "matlab-storm" / "Startup"
        for folder in [matlab_dir, startup_dir, merfish_dir / "codes", export_fig, storm_startup]:
            folder.mkdir(parents=True)
        cfg = MerfishStartupConfig(
            matlab_directory=str(matlab_dir),
            startup_directory=str(startup_dir),
            scratch_path=str(base / "Scratch"),
            python_path=str(base / "Python27"),
            matlab_storm_path=str(base / "matlab-storm"),
            storm_analysis_path=str(base / "storm-analysis"),
            insight_exe=str(base / "InsightM.exe"),
            merfish_analysis_path=str(merfish_dir),
            export_fig_path=str(export_fig),
        )
        report = configure_merfish_environment(cfg, strict=False, verbose=False)
        require(os.environ["MERFISH_ANALYSIS_PATH"] == str(merfish_dir), "MERFISH_ANALYSIS_PATH not set")
        require(str(merfish_dir) in sys.path, "MERFISH root not added")
        require(str(merfish_dir / "codes") in sys.path, "recursive MERFISH child path not added")
        require(report["missing_paths"] == {}, "existing test paths should not be missing")


if __name__ == "__main__":
    test_add_path_recursive()
    test_configure_environment()
    print("startup folder tests passed")
