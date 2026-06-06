"""Shared Python utilities for the MERFISH analysis/SLURM_scripts translation.

The original MATLAB files are cluster driver scripts. These helpers provide the
Python equivalents for SLURM environment parsing, decoder loading, logging,
filesystem creation, and light SLURM script generation.
"""

from __future__ import annotations

import contextlib
import csv
import importlib
import os
import pickle
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence


def page_break() -> None:
    """Print a visual separator equivalent to MATLAB PageBreak()."""
    print("=" * 80)


def now_string() -> str:
    """Return a MATLAB datestr-like timestamp."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def ensure_dir(path: os.PathLike[str] | str) -> Path:
    """Create a directory if needed and return it as a Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def normalized_join(base: os.PathLike[str] | str, *parts: str) -> Path:
    """Join paths while accepting MATLAB-style strings with trailing separators."""
    return Path(str(base)).joinpath(*parts)


@contextlib.contextmanager
def diary(log_file: os.PathLike[str] | str):
    """Mirror MATLAB diary by teeing stdout/stderr to a log file inside a context."""
    log_path = Path(log_file)
    ensure_dir(log_path.parent)
    old_stdout = sys.stdout
    old_stderr = sys.stderr

    class _Tee:
        def __init__(self, *streams: Any) -> None:
            self.streams = streams

        def write(self, text: str) -> int:
            for stream in self.streams:
                stream.write(text)
                stream.flush()
            return len(text)

        def flush(self) -> None:
            for stream in self.streams:
                stream.flush()

    with log_path.open("a", encoding="utf-8") as handle:
        sys.stdout = _Tee(old_stdout, handle)
        sys.stderr = _Tee(old_stderr, handle)
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


def slurm_node_id() -> str:
    """Return the SLURM node list, matching getenv('SLURM_NODELIST')."""
    return os.getenv("SLURM_NODELIST", "")


def slurm_job_id() -> str:
    """Return the SLURM job id, checking both MATLAB-used variants."""
    return os.getenv("SLURM_JOBID") or os.getenv("SLURM_JOB_ID", "")


def slurm_ntasks(default: int = 1) -> int:
    """Return SLURM_NTASKS as an integer."""
    value = os.getenv("SLURM_NTASKS", str(default))
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"SLURM_NTASKS must be an integer, got {value!r}") from exc
    return parsed


def get_fov_id(array_id: Optional[int] = None) -> int:
    """Return the FOV id from an explicit value or SLURM_ARRAY_TASK_ID."""
    if array_id is None:
        raw = os.getenv("SLURM_ARRAY_TASK_ID", "")
        if raw == "":
            raise ValueError("array_id was not provided and SLURM_ARRAY_TASK_ID is empty")
        try:
            return int(raw)
        except ValueError as exc:
            raise ValueError(f"SLURM_ARRAY_TASK_ID must be an integer, got {raw!r}") from exc
    return int(array_id)


def log_run_start(script_name: str, fov_id: Optional[int] = None) -> float:
    """Print the common SLURM start block and return a start time."""
    page_break()
    if fov_id is None:
        print(f"Running {script_name} at {now_string()}")
    else:
        print(f"Running {script_name} on {fov_id} at {now_string()}")
    print(f"Running on {slurm_node_id()}")
    print(f"With job id {slurm_job_id()}")
    return time.perf_counter()


def log_run_end(message: str, start_time: Optional[float] = None) -> None:
    """Print the common completion block."""
    page_break()
    if start_time is not None:
        print(f"...completed in {time.perf_counter() - start_time:.3f} s")
    print(f"{message} at {now_string()}")


def _load_decoder_from_python_module(path: str) -> Any:
    """Load MERFISHDecoder.Load(path) from an available Python translation module."""
    candidate_modules = [
        "MERFISHDecoder",
        "analysis.classes.MERFISHDecoder",
        "python_translation.analysis.classes.MERFISHDecoder",
    ]
    for module_name in candidate_modules:
        try:
            module = importlib.import_module(module_name)
        except ImportError:
            continue
        decoder_class = getattr(module, "MERFISHDecoder", None)
        if decoder_class is not None and hasattr(decoder_class, "Load"):
            return decoder_class.Load(path)
        if hasattr(module, "Load"):
            return module.Load(path)
    raise RuntimeError(
        "Could not import a Python MERFISHDecoder with a Load(path) method. "
        "Pass decoder_loader=callable or provide the translated analysis/classes package."
    )


def load_decoder(
    path: os.PathLike[str] | str,
    decoder_loader: Optional[Callable[[str], Any]] = None,
) -> Any:
    """Load the MERFISH decoder, equivalent to MERFISHDecoder.Load(path)."""
    path_str = str(path)
    if decoder_loader is not None:
        return decoder_loader(path_str)
    pickle_path = Path(path_str) / "mDecoder.pkl"
    if pickle_path.exists():
        with pickle_path.open("rb") as handle:
            return pickle.load(handle)
    return _load_decoder_from_python_module(path_str)


def call_decoder_method(decoder: Any, method_name: str, *args: Any, **kwargs: Any) -> Any:
    """Call a MERFISHDecoder method with a clear error when it is missing."""
    method = getattr(decoder, method_name, None)
    if method is None or not callable(method):
        raise AttributeError(f"Decoder object has no callable method {method_name!r}")
    return method(*args, **kwargs)


def set_decoder_overwrite(decoder: Any, value: bool) -> None:
    """Set decoder.overwrite if the translated decoder supports the attribute."""
    setattr(decoder, "overwrite", bool(value))


def set_parallel_if_available(decoder: Any, workers: Optional[int]) -> None:
    """Call SetParallel(workers) when available; otherwise store workers on the decoder."""
    if workers is None:
        return
    method = getattr(decoder, "SetParallel", None)
    if callable(method):
        method(workers)
    else:
        setattr(decoder, "parallel_workers", workers)


def decoder_fov_ids(decoder: Any) -> List[int]:
    """Return fov ids from decoder.fovIDs or decoder.fov_ids."""
    ids = getattr(decoder, "fovIDs", None)
    if ids is None:
        ids = getattr(decoder, "fov_ids", None)
    if ids is None:
        return []
    return [int(v) for v in ids]


def decoder_normalized_data_path(decoder: Any, fallback: os.PathLike[str] | str) -> str:
    """Return normalized data path from a decoder object, with a fallback path."""
    value = getattr(decoder, "normalizedDataPath", None)
    if value is None:
        value = getattr(decoder, "normalized_data_path", None)
    return str(value if value is not None else fallback)


def decoder_report_path(decoder: Any) -> str:
    """Return the report path suffix from the decoder, matching mDecoder.reportPath."""
    value = getattr(decoder, "reportPath", None)
    if value is None:
        value = getattr(decoder, "report_path", "reports")
    return str(value).strip("/\\")


def nested_get(obj: Any, dotted_path: str, default: Any = None) -> Any:
    """Get nested attributes or mapping keys using dotted names."""
    current = obj
    for part in dotted_path.split("."):
        if isinstance(current, Mapping):
            current = current.get(part, default)
        else:
            current = getattr(current, part, default)
        if current is default:
            return default
    return current


def build_fov_barcode_files(nDPath: os.PathLike[str] | str) -> List[Dict[str, Any]]:
    """Find barcode_fov files matching MATLAB BuildFileStructure usage."""
    barcode_dir = Path(nDPath) / "barcodes" / "barcode_fov"
    pattern = re.compile(r"fov_(?P<fov>[0-9]+)_blist.*\.bin$")
    records: List[Dict[str, Any]] = []
    if not barcode_dir.exists():
        return records
    for file_path in sorted(barcode_dir.glob("*.bin")):
        match = pattern.search(file_path.name)
        if match:
            records.append({"fov": int(match.group("fov")), "filePath": str(file_path)})
    return records


def read_binary_file_header(file_path: os.PathLike[str] | str) -> Dict[str, Any]:
    """Read enough of a binary file to report simple corruption status.

    MATLAB used ReadBinaryFileHeader from the fileIO folder. This Python version
    performs the scheduler-level validation used here: existence and nonzero size.
    """
    path = Path(file_path)
    exists = path.exists()
    size = path.stat().st_size if exists else 0
    return {"isCorrupt": (not exists) or size == 0, "size": size, "filePath": str(path)}


def save_csv_matrix(file_path: os.PathLike[str] | str, rows: Sequence[Sequence[Any]]) -> None:
    """Write a numeric/text matrix to CSV."""
    path = Path(file_path)
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


@dataclass
class SlurmJob:
    """Python representation of the MATLAB SLURMJob used by MERFISHScheduler."""

    command: Sequence[str]
    name: str
    script_path: Path
    script_name: str
    output_log: Path
    error_log: Path
    open_mode: str = "append"
    constraint: str = ""
    time_limit_minutes: int = 60
    memory_limit_mb: int = 4000
    cpus_per_task: int = 1
    ntasks: int = 1
    partition: str = ""
    exclude: str = ""
    number_resubmit: int = 0
    timer_period_seconds: int = 300
    pre_check: bool = False
    verbose: bool = True
    dependencies: List["SlurmJob"] = field(default_factory=list)

    def script_file(self) -> Path:
        return self.script_path / self.script_name

    def render(self) -> str:
        """Render an sbatch script equivalent to the MATLAB SLURMJob settings."""
        lines = ["#!/usr/bin/env bash", f"#SBATCH --job-name={self.name}"]
        lines.append(f"#SBATCH --output={self.output_log}")
        lines.append(f"#SBATCH --error={self.error_log}")
        lines.append(f"#SBATCH --time={self.time_limit_minutes}")
        lines.append(f"#SBATCH --mem={self.memory_limit_mb}")
        if self.ntasks:
            lines.append(f"#SBATCH --ntasks={self.ntasks}")
        if self.cpus_per_task:
            lines.append(f"#SBATCH --cpus-per-task={self.cpus_per_task}")
        if self.partition:
            lines.append(f"#SBATCH --partition={self.partition}")
        if self.constraint:
            lines.append(f"#SBATCH --constraint={self.constraint}")
        if self.exclude:
            lines.append(f"#SBATCH --exclude={self.exclude}")
        lines.extend(["set -euo pipefail", ""])
        lines.extend(self.command)
        lines.append("")
        return "\n".join(lines)

    def write_script(self) -> Path:
        """Write the sbatch script and return its path."""
        ensure_dir(self.script_path)
        ensure_dir(self.output_log.parent)
        ensure_dir(self.error_log.parent)
        script = self.script_file()
        script.write_text(self.render(), encoding="utf-8")
        script.chmod(0o755)
        if self.verbose:
            print(f"Wrote SLURM script: {script}")
        return script

    def submit(self, dry_run: bool = True) -> Optional[str]:
        """Submit the job with sbatch, or only render it when dry_run is True."""
        script = self.write_script()
        if dry_run:
            return None
        result = subprocess.run(["sbatch", str(script)], check=True, text=True, capture_output=True)
        return result.stdout.strip()


@dataclass
class SlurmJobArray:
    """Container equivalent to MATLAB SLURMJobArray."""

    jobs: List[SlurmJob]
    name: str
    verbose: bool = True

    def write_scripts(self) -> List[Path]:
        paths = [job.write_script() for job in self.jobs]
        if self.verbose:
            print(f"Wrote {len(paths)} scripts for job array {self.name!r}")
        return paths

    def submit(self, dry_run: bool = True) -> List[Optional[str]]:
        return [job.submit(dry_run=dry_run) for job in self.jobs]
