"""Translate analysis/SLURM_scripts/MERFISHScheduler.m.

Purpose/use: coordinate the MERFISH analysis workflow on a SLURM cluster by
creating preprocessing, optimization, decoding, segmentation, parsing,
performance, counting, barcode metadata, doublet-score, summation, combine, and
mosaic jobs.
"""

from __future__ import annotations

import argparse
import importlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from .slurm_utils import SlurmJob, SlurmJobArray, decoder_fov_ids, ensure_dir, load_decoder, log_run_end, log_run_start


@dataclass
class SchedulerResult:
    """Return object containing all generated SLURM job arrays."""

    job_arrays: Dict[str, SlurmJobArray] = field(default_factory=dict)
    single_jobs: Dict[str, SlurmJob] = field(default_factory=dict)
    scripts: List[Path] = field(default_factory=list)


def _load_decoder_class() -> Any:
    candidates = [
        "MERFISHDecoder",
        "analysis.classes.MERFISHDecoder",
        "python_translation.analysis.classes.MERFISHDecoder",
    ]
    for module_name in candidates:
        try:
            module = importlib.import_module(module_name)
        except ImportError:
            continue
        cls = getattr(module, "MERFISHDecoder", None)
        if cls is not None:
            return cls
    raise RuntimeError(
        "Could not import translated MERFISHDecoder. Pass decoder_factory or set skip_decoder_construction=True."
    )


def _construct_decoder(
    rDPath: str,
    nDPath: str,
    parameters: Optional[dict[str, Any]],
    decoder_factory: Optional[Callable[[str, str, Optional[dict[str, Any]]], Any]],
) -> Any:
    """Create and save a MERFISHDecoder equivalent to MERFISHDecoder(rDPath,nDPath)."""
    if decoder_factory is not None:
        decoder = decoder_factory(rDPath, nDPath, parameters)
    else:
        cls = _load_decoder_class()
        decoder = cls(rDPath, nDPath, parameters=parameters)
    save_method = getattr(decoder, "Save", None)
    if callable(save_method):
        save_method()
    return decoder


def _python_command(module_name: str, *args: str, options: Optional[dict[str, Any]] = None) -> str:
    quoted = [repr(str(arg)) for arg in args]
    extra: List[str] = []
    if options:
        for key, value in options.items():
            flag = "--" + key.replace("_", "-")
            if isinstance(value, bool):
                if value:
                    extra.append(flag)
            elif value is not None:
                extra.extend([flag, repr(str(value))])
    return "python -m " + module_name + " " + " ".join(quoted + extra)


def _job(
    *,
    command: List[str],
    name: str,
    scheduler_path: Path,
    subdir: str,
    script_name: str,
    output_log: Path,
    error_log: Path,
    partition: str,
    time_limit: int,
    memory: int,
    ntasks: int = 1,
    cpus_per_task: int = 1,
    constraint: str = "holyib",
    exclude: str = "",
    open_mode: str = "append",
    timer_period: int = 300,
    number_resubmit: int = 0,
    pre_check: bool = False,
) -> SlurmJob:
    return SlurmJob(
        command=command,
        name=name,
        script_path=scheduler_path / subdir,
        script_name=script_name,
        output_log=output_log,
        error_log=error_log,
        open_mode=open_mode,
        constraint=constraint,
        time_limit_minutes=time_limit,
        memory_limit_mb=memory,
        cpus_per_task=cpus_per_task,
        ntasks=ntasks,
        partition=partition,
        exclude=exclude,
        timer_period_seconds=timer_period,
        number_resubmit=number_resubmit,
        pre_check=pre_check,
    )


def merfish_scheduler(
    nDPath: str,
    rDPath: str = "",
    aPath: str = "",
    aControl: str = "mwodscpfnrilbu",
    parameters: Optional[dict[str, Any]] = None,
    skip_decoder_construction: bool = False,
    small_job_queue: str = "shared,zhuang",
    big_job_queue: str = "zhuang,general,shared",
    load_matlab_command: str = "module load matlab/R2017a-fasrc02",
    load_python_command: str = "module load python/2.7.14-fasrc01",
    activate_environment_command: str = "",
    decoder_loader: Optional[Callable[[str], Any]] = None,
    decoder_factory: Optional[Callable[[str, str, Optional[dict[str, Any]]], Any]] = None,
    dry_run: bool = True,
) -> SchedulerResult:
    """Create the SLURM workflow scripts equivalent to MERFISHScheduler.m.

    The output is a Python-native job plan. When dry_run is True, scripts are
    written but not submitted with sbatch.
    """
    if not nDPath:
        raise ValueError("A normalized data path must be provided")
    if "m" in aControl and not rDPath:
        raise ValueError("A raw data path must be provided to create a MERFISHDecoder instance")
    if "f" in aControl and not aPath:
        raise ValueError("A path to abundance data must be provided for performance metrics")

    npath = Path(nDPath)
    ensure_dir(npath)
    log_path = ensure_dir(npath / "log")
    scheduler_path = ensure_dir(npath / "scheduler")
    result = SchedulerResult()

    start = log_run_start("MERFISHScheduler")
    print("Running MERFISH Scheduler")
    print(f"Raw data path: {rDPath}")
    print(f"Analyzed data path: {nDPath}")
    print(f"Requested analysis: {aControl}")
    print(f"Small jobs queue: {small_job_queue}")
    print(f"Big jobs queue: {big_job_queue}")
    print(f"Provided parameters: {parameters}")

    if "m" in aControl and not skip_decoder_construction:
        if (npath / "mDecoder").exists():
            raise FileExistsError("A decoder already exists.")
        print("Creating MERFISHDecoder")
        _construct_decoder(rDPath, nDPath, parameters, decoder_factory)

    decoder = load_decoder(nDPath, decoder_loader=decoder_loader)
    fov_ids = decoder_fov_ids(decoder)

    prelude = [cmd for cmd in [load_python_command, activate_environment_command] if cmd]
    legacy_matlab_prelude = [cmd for cmd in [load_matlab_command] if cmd]

    if "w" in aControl:
        jobs: List[SlurmJob] = []
        local_log = ensure_dir(log_path / "process")
        for index, fov_id in enumerate(fov_ids, start=1):
            cmd = prelude + [_python_command("analysis.SLURM_scripts.ProcessFOV", nDPath, options={"array_id": fov_id})]
            jobs.append(_job(command=cmd, name=f"Preprocess task {index} for fov {fov_id}", scheduler_path=scheduler_path, subdir="preprocess", script_name=f"w_f{fov_id}.slurm", output_log=local_log / f"fov_process_{fov_id}.out", error_log=local_log / f"fov_process_{fov_id}.err", partition=small_job_queue, time_limit=12 * 60, memory=12000, exclude="holyzhuang01,holy2c14407", number_resubmit=5))
        result.job_arrays["preprocess"] = SlurmJobArray(jobs, name="preprocess")

    if "o" in aControl:
        cmd = prelude + ["mkdir -p /scratch/$USER/$SLURM_JOB_ID", _python_command("analysis.SLURM_scripts.Optimize", nDPath, options={"overwrite": True}), "rm -rf /scratch/$USER/$SLURM_JOB_ID"]
        result.single_jobs["optimize"] = _job(command=cmd, name="Optimize", scheduler_path=scheduler_path, subdir="optimize", script_name="optimize.slurm", output_log=log_path / "optimize.out", error_log=log_path / "optimize.err", partition=big_job_queue, time_limit=24 * 60, memory=255000, ntasks=10, timer_period=20 * 60)

    if "d" in aControl:
        jobs = []
        local_log = ensure_dir(log_path / "decoding")
        for index, fov_id in enumerate(fov_ids, start=1):
            cmd = prelude + [_python_command("analysis.SLURM_scripts.DecodeFOV", nDPath, options={"array_id": fov_id})]
            jobs.append(_job(command=cmd, name=f"Decoding task {index} for fov {fov_id}", scheduler_path=scheduler_path, subdir="decoding", script_name=f"d_f{fov_id}.slurm", output_log=local_log / f"fov_decode_{fov_id}.out", error_log=local_log / f"fov_decode_{fov_id}.err", partition=small_job_queue, time_limit=6 * 60, memory=32000, exclude="holyzhuang01,holy2c14407", number_resubmit=5))
        result.job_arrays["decoding"] = SlurmJobArray(jobs, name="decoding")

    if "s" in aControl:
        jobs = []
        local_log = ensure_dir(log_path / "segment")
        for index, fov_id in enumerate(fov_ids, start=1):
            cmd = prelude + [_python_command("analysis.SLURM_scripts.SegmentFOV", nDPath, options={"array_id": fov_id})]
            jobs.append(_job(command=cmd, name=f"Segmentation task {index} for fov {fov_id}", scheduler_path=scheduler_path, subdir="segmentation", script_name=f"s_f{fov_id}.slurm", output_log=local_log / f"fov_segment_{fov_id}.out", error_log=local_log / f"fov_segment_{fov_id}.err", partition=small_job_queue, time_limit=30, memory=8000, exclude="holyzhuang01,holy2c14407", number_resubmit=5))
        result.job_arrays["segment"] = SlurmJobArray(jobs, name="segment")

    if "c" in aControl:
        cmd = prelude + [_python_command("analysis.SLURM_scripts.CombineFoundFeatures", nDPath, options={"overwrite": True})]
        result.single_jobs["combine"] = _job(command=cmd, name="Combine", scheduler_path=scheduler_path, subdir="combine", script_name="combine.slurm", output_log=log_path / "combine.out", error_log=log_path / "combine.err", partition=big_job_queue, time_limit=12 * 60, memory=64000)

    if "l" in aControl:
        cmd = prelude + [_python_command("analysis.SLURM_scripts.LowResMosaic", nDPath)]
        result.single_jobs["mosaic"] = _job(command=cmd, name="Mosaic", scheduler_path=scheduler_path, subdir="mosaic", script_name="mosaic.slurm", output_log=log_path / "mosaic.out", error_log=log_path / "mosaic.err", partition=big_job_queue, time_limit=12 * 60, memory=127000, number_resubmit=3, pre_check=True)

    if "p" in aControl:
        jobs = []
        local_log = ensure_dir(log_path / "parse")
        for index, fov_id in enumerate(fov_ids, start=1):
            cmd = prelude + [_python_command("analysis.SLURM_scripts.ParseFOV", nDPath, options={"array_id": fov_id})]
            jobs.append(_job(command=cmd, name=f"Parse task {index} for fov {fov_id}", scheduler_path=scheduler_path, subdir="parse", script_name=f"p_f{fov_id}.slurm", output_log=local_log / f"fov_parse_{fov_id}.out", error_log=local_log / f"fov_parse_{fov_id}.err", partition=small_job_queue, time_limit=3 * 60, memory=24000, exclude="holyzhuang01,holy2c14407", number_resubmit=5))
        result.job_arrays["parse"] = SlurmJobArray(jobs, name="parse")

    if "f" in aControl:
        cmd = prelude + ["mkdir -p /scratch/$USER/$SLURM_JOB_ID", _python_command("analysis.SLURM_scripts.Performance", nDPath, aPath, "/scratch/$USER/"), "rm -rf /scratch/$USER/$SLURM_JOB_ID"]
        result.single_jobs["performance"] = _job(command=cmd, name="Performance", scheduler_path=scheduler_path, subdir="performance", script_name="performance.slurm", output_log=log_path / "performance.out", error_log=log_path / "performance.err", partition=big_job_queue, time_limit=4 * 60, memory=255000, ntasks=12, open_mode="truncate", pre_check=False)

    if "n" in aControl:
        cmd = prelude + ["mkdir -p /scratch/$USER/$SLURM_JOB_ID", _python_command("analysis.SLURM_scripts.CalculateNumbers", nDPath, "/scratch/$USER/"), "rm -rf /scratch/$USER/$SLURM_JOB_ID"]
        result.single_jobs["numbers"] = _job(command=cmd, name="Numbers", scheduler_path=scheduler_path, subdir="numbers", script_name="numbers.slurm", output_log=log_path / "numbers.out", error_log=log_path / "numbers.err", partition=big_job_queue, time_limit=10 * 60, memory=128000, ntasks=12, open_mode="truncate", pre_check=False)

    if "b" in aControl:
        cmd = prelude + ["mkdir -p /scratch/$USER/$SLURM_JOB_ID", _python_command("analysis.SLURM_scripts.ExportBarcodeMetadata", nDPath), "rm -rf /scratch/$USER/$SLURM_JOB_ID"]
        result.single_jobs["barcode_metadata"] = _job(command=cmd, name="Export Barcodes", scheduler_path=scheduler_path, subdir="barcode_metadata", script_name="barcode_metadata.slurm", output_log=log_path / "barcode_metadata.out", error_log=log_path / "barcode_metadata.err", partition=big_job_queue, time_limit=10 * 60, memory=250000, ntasks=12, open_mode="truncate", pre_check=False)

    if "u" in aControl:
        cmd = prelude + [_python_command("analysis.SLURM_scripts.CalculateDoubletScore", nDPath)]
        result.single_jobs["doublet_score"] = _job(command=cmd, name="Calculate Doublet Score", scheduler_path=scheduler_path, subdir="doublet_score", script_name="doublet_score.slurm", output_log=log_path / "doublet_score.out", error_log=log_path / "doublet_score.err", partition=big_job_queue, time_limit=10 * 60, memory=32000, open_mode="truncate", pre_check=False)

    if "r" in aControl:
        jobs = []
        local_log = ensure_dir(log_path / "sum")
        for index, fov_id in enumerate(fov_ids, start=1):
            cmd = prelude + [_python_command("analysis.SLURM_scripts.SumFOV", nDPath, options={"array_id": fov_id})]
            jobs.append(_job(command=cmd, name=f"Sum task {index} for fov {fov_id}", scheduler_path=scheduler_path, subdir="sum", script_name=f"r_f{fov_id}.slurm", output_log=local_log / f"fov_sum_{fov_id}.out", error_log=local_log / f"fov_sum_{fov_id}.err", partition=small_job_queue, time_limit=3 * 60, memory=24000, exclude="holyzhuang01,holy2c14407", number_resubmit=5))
        result.job_arrays["sum"] = SlurmJobArray(jobs, name="sum")

    if "i" in aControl:
        cmd = prelude + [_python_command("analysis.SLURM_scripts.CombineSum", nDPath, options={"overwrite": True})]
        result.single_jobs["combine_sum"] = _job(command=cmd, name="Combine sum", scheduler_path=scheduler_path, subdir="combine_sum", script_name="combine_sum.slurm", output_log=log_path / "combine_sum.out", error_log=log_path / "combine_sum.err", partition=big_job_queue, time_limit=3 * 60, memory=64000)

    for array in result.job_arrays.values():
        result.scripts.extend(array.write_scripts())
        if not dry_run:
            array.submit(dry_run=False)
    for job in result.single_jobs.values():
        result.scripts.append(job.write_script())
        if not dry_run:
            job.submit(dry_run=False)

    log_run_end("Completed MERFISH Scheduler", start)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Create MERFISH SLURM workflow scripts.")
    parser.add_argument("nDPath")
    parser.add_argument("--raw-data-path", dest="rDPath", default="")
    parser.add_argument("--abundance-path", dest="aPath", default="")
    parser.add_argument("--analysis-control", dest="aControl", default="mwodscpfnrilbu")
    parser.add_argument("--skip-decoder-construction", action="store_true")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args()
    merfish_scheduler(
        nDPath=args.nDPath,
        rDPath=args.rDPath,
        aPath=args.aPath,
        aControl=args.aControl,
        skip_decoder_construction=args.skip_decoder_construction,
        dry_run=not args.submit,
    )


if __name__ == "__main__":
    main()
