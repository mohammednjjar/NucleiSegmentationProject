"""Python translation of analysis/classes/SLURMJobArray.m.

Coordinates related SLURMJob objects, tracks submitted/completed/failed flags,
submits jobs, resubmits failed jobs, cancels jobs, and exposes a process-local
registry like the MATLAB global registry.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Tuple
import datetime as _dt
import time

from .SLURMJob import SLURMJob

_JOB_ARRAY_REGISTRY: List["SLURMJobArray"] = []


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


@dataclass
class SLURMJobArray:
    jobs: Iterable[SLURMJob]
    name: str = ""
    verbose: bool = True

    completed: bool = False
    submitted: bool = False
    failed: bool = False
    jobsSubmitted: List[bool] = field(default_factory=list)
    jobsCompleted: List[bool] = field(default_factory=list)
    jobsFailed: List[bool] = field(default_factory=list)
    startTime: str = ""
    endTime: str = ""
    duration: float = 0.0
    _start_epoch: float = 0.0

    def __post_init__(self) -> None:
        self.jobs = list(self.jobs)
        self.numJobs = len(self.jobs)
        self.jobsSubmitted = [False] * self.numJobs
        self.jobsCompleted = [False] * self.numJobs
        self.jobsFailed = [False] * self.numJobs
        if not self.name:
            self.name = f"job_array_{id(self):x}"
        self.RegisterJobArray(self)

    def Submit(self) -> None:
        self.submitted = True
        self.startTime = _now()
        self._start_epoch = time.time()
        for idx, job in enumerate(self.jobs):
            ok = job.Submit()
            self.jobsSubmitted[idx] = bool(ok or job.submitted)
            self.jobsFailed[idx] = bool(job.failed)
            self.jobsCompleted[idx] = bool(job.completed)
        self._refresh_state()

    submit = Submit

    def Resubmit(self) -> None:
        for idx, job in enumerate(self.jobs):
            if job.failed or not job.completed:
                ok = job.Resubmit()
                self.jobsSubmitted[idx] = bool(ok or job.submitted)
                self.jobsFailed[idx] = bool(job.failed)
                self.jobsCompleted[idx] = bool(job.completed)
        self._refresh_state()

    resubmit = Resubmit

    def Cancel(self) -> None:
        for idx, job in enumerate(self.jobs):
            job.Cancel()
            self.jobsFailed[idx] = bool(job.failed)
        self.failed = True
        self.endTime = _now()
        self.duration = time.time() - (self._start_epoch or time.time())

    cancel = Cancel

    def AddSubmitTrigger(self, sourceObj, eventName: str) -> None:
        listeners = getattr(self, "jobSubmitListeners", [])
        listeners.append((sourceObj, eventName))
        self.jobSubmitListeners = listeners

    add_submit_trigger = AddSubmitTrigger

    def HandleSubmitTrigger(self, source=None, event=None, listenerID=None) -> None:
        self.Submit()

    handle_submit_trigger = HandleSubmitTrigger

    def AddJob(self, job: SLURMJob) -> None:
        if not isinstance(job, SLURMJob):
            raise TypeError("job must be an SLURMJob")
        self.jobs.append(job)
        self.numJobs = len(self.jobs)
        self.jobsSubmitted.append(False)
        self.jobsCompleted.append(False)
        self.jobsFailed.append(False)

    add_job = AddJob

    def MarkJobSubmitted(self, jobID: int) -> None:
        self.jobsSubmitted[int(jobID)] = True
        self._refresh_state()

    mark_job_submitted = MarkJobSubmitted

    def MarkJobComplete(self, jobID: int) -> None:
        self.jobsCompleted[int(jobID)] = True
        self._refresh_state()

    mark_job_complete = MarkJobComplete

    def MarkJobFailed(self, jobID: int) -> None:
        self.jobsFailed[int(jobID)] = True
        self._refresh_state()

    mark_job_failed = MarkJobFailed

    def _refresh_state(self) -> None:
        for idx, job in enumerate(self.jobs):
            self.jobsSubmitted[idx] = bool(job.submitted)
            self.jobsCompleted[idx] = bool(job.completed)
            self.jobsFailed[idx] = bool(job.failed)
        self.completed = all(self.jobsCompleted) if self.jobsCompleted else False
        self.failed = any(self.jobsFailed) if self.jobsFailed else False
        if self.completed or self.failed:
            self.endTime = _now()
            self.duration = time.time() - (self._start_epoch or time.time())

    def Status(self) -> Tuple[str, str]:
        self._refresh_state()
        n_sub = sum(self.jobsSubmitted)
        n_done = sum(self.jobsCompleted)
        n_failed = sum(self.jobsFailed)
        short = f"{self.name}: {n_done}/{self.numJobs} complete, {n_failed} failed"
        long = [
            f"submitted: {self.submitted}",
            f"completed: {self.completed}",
            f"failed: {self.failed}",
            f"jobs_submitted: {n_sub}",
            f"jobs_completed: {n_done}",
            f"jobs_failed: {n_failed}",
            f"duration_seconds: {self.duration:.1f}",
        ]
        return short, "\n".join(long)

    status = Status

    @staticmethod
    def RegisterJobArray(jobArrayObject: "SLURMJobArray") -> None:
        if jobArrayObject not in _JOB_ARRAY_REGISTRY:
            _JOB_ARRAY_REGISTRY.append(jobArrayObject)

    @staticmethod
    def GetJobArrayRegistry() -> List["SLURMJobArray"]:
        return list(_JOB_ARRAY_REGISTRY)

    @staticmethod
    def ClearRegistry() -> None:
        _JOB_ARRAY_REGISTRY.clear()

    @staticmethod
    def CancelAllJobArrays() -> None:
        for array in list(_JOB_ARRAY_REGISTRY):
            array.Cancel()
