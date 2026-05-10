"""Python translation of analysis/classes/SLURMJob.m.

SLURMJob writes shell scripts, submits them with sbatch, polls status via squeue,
can cancel/requeue jobs, records command history, and exposes MATLAB-style method
names for compatibility with the original pipeline design.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Tuple, Any
import datetime as _dt
import os
import subprocess
import threading
import time
import uuid


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


@dataclass
class SLURMJob:
    scriptText: Iterable[str]
    name: str = ""
    jobName: str = ""
    queue: str = ""
    memory: str = "1G"
    nodes: int = 1
    ntasks: int = 1
    cpusPerTask: int = 1
    walltime: str = "01:00:00"
    logPath: str = ""
    scriptPath: str = "."
    uniqueID: str = field(default_factory=lambda: uuid.uuid4().hex[:12])
    email: str = ""
    emailType: str = "FAIL"
    preCheck: bool = False
    completeFcn: Optional[Callable[["SLURMJob"], bool]] = None
    verbose: bool = True
    dryRun: bool = False
    timerPeriod: float = 0.0

    scriptName: str = ""
    jobID: str = ""
    submitted: bool = False
    numFailures: int = 0
    failed: bool = False
    completed: bool = False
    user_canceled: bool = False
    state: str = ""
    duration: float = 0.0
    startTime: str = ""
    endTime: str = ""
    history: List[Tuple[str, str]] = field(default_factory=list)
    command_line_history: List[Tuple[str, int, str, str]] = field(default_factory=list)
    numStatusChecks: int = 0
    _timer: Optional[threading.Timer] = field(default=None, repr=False)
    _submit_epoch: Optional[float] = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self.scriptText = list(self.scriptText) if not isinstance(self.scriptText, str) else [self.scriptText]
        if not self.jobName:
            self.jobName = self.name or f"slurm_job_{self.uniqueID}"
        if not self.scriptName:
            self.scriptName = f"{self.jobName}_{self.uniqueID}.sh".replace(os.sep, "_")
        Path(self.scriptPath).mkdir(parents=True, exist_ok=True)
        if self.logPath:
            Path(self.logPath).mkdir(parents=True, exist_ok=True)

    @property
    def script_file(self) -> Path:
        return Path(self.scriptPath) / self.scriptName

    def _record(self, state: str) -> None:
        self.state = state
        self.history.append((state, _now()))

    def _run(self, cmd: List[str]) -> subprocess.CompletedProcess:
        proc = subprocess.run(cmd, text=True, capture_output=True)
        self.command_line_history.append((" ".join(cmd), proc.returncode, proc.stdout, proc.stderr))
        return proc

    def WriteScript(self, overwrite: bool = True) -> Path:
        path = self.script_file
        if path.exists() and not overwrite:
            return path
        lines = ["#!/usr/bin/env bash", f"#SBATCH --job-name={self.jobName}"]
        if self.queue:
            lines.append(f"#SBATCH --partition={self.queue}")
        if self.walltime:
            lines.append(f"#SBATCH --time={self.walltime}")
        if self.nodes:
            lines.append(f"#SBATCH --nodes={int(self.nodes)}")
        if self.ntasks:
            lines.append(f"#SBATCH --ntasks={int(self.ntasks)}")
        if self.cpusPerTask:
            lines.append(f"#SBATCH --cpus-per-task={int(self.cpusPerTask)}")
        if self.memory:
            lines.append(f"#SBATCH --mem={self.memory}")
        if self.logPath:
            lines.append(f"#SBATCH --output={Path(self.logPath) / (self.jobName + '_%j.out')}")
            lines.append(f"#SBATCH --error={Path(self.logPath) / (self.jobName + '_%j.err')}")
        if self.email:
            lines.append(f"#SBATCH --mail-user={self.email}")
            lines.append(f"#SBATCH --mail-type={self.emailType}")
        lines.append("set -euo pipefail")
        lines.extend(str(x) for x in self.scriptText)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        path.chmod(0o755)
        return path

    write_script = WriteScript

    def Submit(self, dryRun: Optional[bool] = None) -> bool:
        if self.preCheck and self.completeFcn is not None and bool(self.completeFcn(self)):
            self.completed = True
            self._record("COMPLETED")
            return True
        path = self.WriteScript()
        dry = self.dryRun if dryRun is None else dryRun
        self.startTime = _now()
        self._submit_epoch = time.time()
        if dry:
            self.submitted = True
            self.jobID = f"DRYRUN-{self.uniqueID}"
            self._record("SUBMITTED")
            return True
        proc = self._run(["sbatch", str(path)])
        success = proc.returncode == 0
        if success:
            self.submitted = True
            text = (proc.stdout or "").strip().split()
            self.jobID = text[-1] if text else ""
            self._record("SUBMITTED")
            self.StartTimer()
        else:
            self.failed = True
            self.numFailures += 1
            self._record("FAILED_TO_SUBMIT")
        return success

    submit = Submit

    def Resubmit(self, **kwargs) -> bool:
        self.submitted = False
        self.completed = False
        self.failed = False
        self.user_canceled = False
        self.numFailures += 1
        for k, v in kwargs.items():
            setattr(self, k, v)
        return self.Submit()

    resubmit = Resubmit

    def Requeue(self) -> bool:
        if not self.jobID:
            return self.Resubmit()
        proc = self._run(["scontrol", "requeue", str(self.jobID)])
        ok = proc.returncode == 0
        self._record("REQUEUED" if ok else "REQUEUE_FAILED")
        return ok

    requeue = Requeue

    def StartTimer(self) -> None:
        if self.timerPeriod and self.timerPeriod > 0:
            self.StopTimer()
            self._timer = threading.Timer(float(self.timerPeriod), self._timer_tick)
            self._timer.daemon = True
            self._timer.start()

    start_timer = StartTimer

    def _timer_tick(self) -> None:
        self.UpdateStatus()
        if self.submitted and not (self.completed or self.failed or self.user_canceled):
            self.StartTimer()

    def StopTimer(self) -> None:
        if self._timer is not None:
            self._timer.cancel()
            self._timer = None

    stop_timer = StopTimer

    def UpdateStatus(self) -> str:
        self.numStatusChecks += 1
        if self.completeFcn is not None and bool(self.completeFcn(self)):
            self.completed = True
            self.endTime = _now()
            self.duration = time.time() - (self._submit_epoch or time.time())
            self._record("COMPLETED")
            self.StopTimer()
            return self.state
        if self.jobID.startswith("DRYRUN"):
            self._record("DRYRUN")
            return self.state
        if not self.jobID:
            self._record("UNSUBMITTED")
            return self.state
        proc = self._run(["squeue", "-j", str(self.jobID), "-h", "-o", "%T"])
        state = (proc.stdout or "").strip().splitlines()[0] if proc.returncode == 0 and proc.stdout.strip() else "COMPLETED_OR_NOT_FOUND"
        self.HandleState(state)
        return self.state

    update_status = UpdateStatus

    def Cancel(self, userCanceled: bool = True) -> bool:
        self.user_canceled = bool(userCanceled)
        if not self.jobID or self.jobID.startswith("DRYRUN"):
            self._record("CANCELED")
            self.StopTimer()
            return True
        proc = self._run(["scancel", str(self.jobID)])
        ok = proc.returncode == 0
        self._record("CANCELED" if ok else "CANCEL_FAILED")
        self.StopTimer()
        return ok

    cancel = Cancel

    def Display(self, fileToDisplay: str) -> str:
        text = Path(fileToDisplay).read_text(encoding="utf-8", errors="replace")
        if self.verbose:
            print(text)
        return text

    display = Display

    def HandleState(self, state: str) -> None:
        state_upper = (state or "").upper()
        terminal_ok = {"COMPLETED", "COMPLETING", "COMPLETED_OR_NOT_FOUND"}
        terminal_bad = {"FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "BOOT_FAIL"}
        if state_upper in terminal_ok:
            self.completed = True
            self.endTime = _now()
            self.duration = time.time() - (self._submit_epoch or time.time())
            self.StopTimer()
        elif state_upper in terminal_bad:
            self.failed = True
            self.numFailures += 1
            self.endTime = _now()
            self.duration = time.time() - (self._submit_epoch or time.time())
            self.StopTimer()
        self._record(state_upper)

    handle_state = HandleState

    def Status(self):
        short = f"{self.jobName}: {self.state or 'UNKNOWN'}"
        long = [
            f"jobID: {self.jobID}",
            f"submitted: {self.submitted}",
            f"completed: {self.completed}",
            f"failed: {self.failed}",
            f"numFailures: {self.numFailures}",
            f"duration_seconds: {self.duration:.1f}",
        ]
        return short, "\n".join(long)

    status = Status
