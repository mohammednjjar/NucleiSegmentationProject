"""
Python translation of MERFISH_analysis/startup/merfish_startup.m.

Purpose
-------
Configure paths and environment variables needed by the MERFISH_analysis
software stack. The MATLAB original clears/restores MATLAB paths, defines global
paths for matlab-storm, storm-analysis, legacy BLAST/OligoArray, MERFISH_analysis,
and export_fig, then adds those directories to MATLAB's search path.

Python use
----------
Edit a MerfishStartupConfig object with paths on your machine, then call
configure_merfish_environment(config). It updates os.environ and sys.path.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import os
import sys
from typing import Iterable, Optional


@dataclass
class MerfishStartupConfig:
    """Local path configuration equivalent to variables in merfish_startup.m."""

    matlab_directory: Optional[str] = None
    startup_directory: Optional[str] = None
    scratch_path: Optional[str] = None
    python_path: Optional[str] = None
    matlab_storm_path: Optional[str] = None
    storm_analysis_path: Optional[str] = None
    insight_exe: Optional[str] = None
    legacy_blast_path: Optional[str] = None
    oligo_array_exe: Optional[str] = None
    oligo_array_aux_path: Optional[str] = None
    merfish_analysis_path: Optional[str] = None
    export_fig_path: Optional[str] = None

    @classmethod
    def original_windows_defaults(cls) -> "MerfishStartupConfig":
        """Return the hard-coded Windows paths present in the MATLAB script."""
        return cls(
            matlab_directory="C:\\Users\\Jeff.Morgan0\\Documents\\MATLAB",
            startup_directory="C:\\Users\\Jeff.Morgan0\\Dropbox\\ZhuangLab\\MERFISH_Public\\MERFISH_analysis\\startup",
            scratch_path="C:\\Users\\Jeff.Morgan0\\Documents\\Scratch\\",
            python_path="C:\\Python27\\",
            matlab_storm_path="C:\\Users\\Jeff.Morgan0\\Dropbox\\ZhuangLab\\Coding\\Matlab\\matlab-storm\\",
            storm_analysis_path="C:\\Users\\Jeff.Morgan0\\Dropbox\\ZhuangLab\\Coding\\Python\\storm-analysis\\",
            insight_exe="C:\\Utilities\\STORMAnalysis\\Insight3\\InsightM.exe",
            legacy_blast_path="C:\\Users\\Jeff.Morgan0\\Dropbox\\ZhuangLab\\Coding\\Library\\LegacyBLAST",
            oligo_array_exe="C:\\Users\\Jeff.Morgan0\\Dropbox\\ZhuangLab\\Coding\\Library\\OligoArray\\OligoArray2.jar",
            oligo_array_aux_path="C:\\Program Files\\OligoArrayAux\\bin",
            merfish_analysis_path="C:\\Users\\Jeff.Morgan0\\Dropbox\\ZhuangLab\\MERFISH_Public\\MERFISH_analysis\\",
            export_fig_path="C:\\Users\\Jeff.Morgan0\\Dropbox\\ZhuangLab\\Coding\\Matlab\\matlab-functions\\FromGitHub\\export_fig\\",
        )


def _normalize_path(value: Optional[str]) -> Optional[Path]:
    if value is None:
        return None
    stripped = str(value).strip()
    if stripped == "":
        return None
    return Path(stripped).expanduser()


def _iter_directory_tree(root: Path) -> Iterable[Path]:
    yield root
    for current, dir_names, _file_names in os.walk(root):
        dir_names[:] = [name for name in dir_names if name not in {".git", "__pycache__"}]
        for name in sorted(dir_names):
            yield Path(current) / name


def _add_to_sys_path(path: Path, prepend: bool = True) -> None:
    text_path = str(path)
    if text_path in sys.path:
        sys.path.remove(text_path)
    if prepend:
        sys.path.insert(0, text_path)
    else:
        sys.path.append(text_path)


def add_path(path_value: Optional[str], *, prepend: bool = True, recursive: bool = False, strict: bool = False) -> list[str]:
    """Add one directory or a recursive directory tree to sys.path.

    This is the Python equivalent of MATLAB addpath(path) and addpath(genpath(path)).
    Returns the list of paths actually inserted.
    """
    path = _normalize_path(path_value)
    if path is None:
        return []
    if not path.exists():
        if strict:
            raise FileNotFoundError(f"Required path does not exist: {path}")
        return []
    if not path.is_dir():
        if strict:
            raise NotADirectoryError(f"Expected a directory: {path}")
        return []

    paths = list(_iter_directory_tree(path)) if recursive else [path]
    inserted: list[str] = []
    for item in reversed(paths) if prepend else paths:
        _add_to_sys_path(item, prepend=prepend)
        inserted.append(str(item))
    return inserted


def _set_env_path(name: str, value: Optional[str]) -> Optional[str]:
    path = _normalize_path(value)
    if path is None:
        return None
    os.environ[name] = str(path)
    return str(path)


def configure_merfish_environment(
    config: Optional[MerfishStartupConfig] = None,
    *,
    prepend: bool = True,
    strict: bool = False,
    verbose: bool = True,
) -> dict[str, object]:
    """Configure the Python environment for a local MERFISH_analysis checkout.

    Parameters
    ----------
    config:
        Machine-specific paths. If omitted, the original MATLAB Windows defaults
        are used for documentation-equivalent behavior.
    prepend:
        When True, inserted paths go to the front of sys.path, equivalent to
        MATLAB addpath(..., '-begin').
    strict:
        When True, missing paths raise errors. When False, missing paths are
        reported but skipped.
    verbose:
        Print a short report, matching MATLAB display(...) behavior.

    Returns
    -------
    Dictionary containing environment variables, inserted paths, and missing paths.
    """
    cfg = config if config is not None else MerfishStartupConfig.original_windows_defaults()

    env_names = {
        "MERFISH_SCRATCH_PATH": cfg.scratch_path,
        "MERFISH_PYTHON_PATH": cfg.python_path,
        "MERFISH_MATLAB_STORM_PATH": cfg.matlab_storm_path,
        "MERFISH_STORM_ANALYSIS_PATH": cfg.storm_analysis_path,
        "MERFISH_INSIGHT_EXE": cfg.insight_exe,
        "MERFISH_LEGACY_BLAST_PATH": cfg.legacy_blast_path,
        "MERFISH_OLIGO_ARRAY_EXE": cfg.oligo_array_exe,
        "MERFISH_OLIGO_ARRAY_AUX_PATH": cfg.oligo_array_aux_path,
        "MERFISH_ANALYSIS_PATH": cfg.merfish_analysis_path,
        "MERFISH_EXPORT_FIG_PATH": cfg.export_fig_path,
    }
    environment = {name: _set_env_path(name, value) for name, value in env_names.items()}

    path_fields = [
        ("matlab_directory", cfg.matlab_directory, False),
        ("startup_directory", cfg.startup_directory, False),
        ("matlab_storm_startup", str(Path(cfg.matlab_storm_path) / "Startup") if cfg.matlab_storm_path else None, False),
        ("legacy_blast_path", cfg.legacy_blast_path, False),
        ("oligo_array_aux_path", cfg.oligo_array_aux_path, False),
        ("merfish_analysis_path", cfg.merfish_analysis_path, True),
        ("export_fig_path", cfg.export_fig_path, False),
    ]

    inserted: dict[str, list[str]] = {}
    missing: dict[str, str] = {}
    for field_name, value, recursive in path_fields:
        path = _normalize_path(value)
        if path is None:
            inserted[field_name] = []
            continue
        if not path.exists():
            missing[field_name] = str(path)
        inserted[field_name] = add_path(value, prepend=prepend, recursive=recursive, strict=strict)

    report = {
        "config": asdict(cfg),
        "environment": environment,
        "inserted_paths": inserted,
        "missing_paths": missing,
    }

    if verbose:
        print("Adding MERFISH_analysis")
        print(f"  {cfg.merfish_analysis_path}")
        print("  And all enclosed paths" if inserted.get("merfish_analysis_path") else "  No MERFISH paths inserted")
        if missing:
            print("Missing paths skipped:")
            for name, path_text in missing.items():
                print(f"  {name}: {path_text}")

    return report


def main() -> None:
    """Command-line entry point using original MATLAB default paths."""
    configure_merfish_environment(verbose=True)


if __name__ == "__main__":
    main()
