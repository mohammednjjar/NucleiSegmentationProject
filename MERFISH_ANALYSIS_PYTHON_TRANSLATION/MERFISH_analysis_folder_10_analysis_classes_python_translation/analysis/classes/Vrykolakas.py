"""Python translation of analysis/classes/Vrykolakas.m.

Keeps an interactive SSH/shell session active by printing a periodic message.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import threading


@dataclass
class Vrykolakas:
    timerPeriod: float = 10 * 60
    aliveString: str = "I'm alive!"
    verbose: bool = True
    statusTimer: threading.Timer | None = field(default=None, init=False, repr=False)
    numStatusChecks: int = 0

    def __post_init__(self) -> None:
        self.StartTimer()

    def _tick(self) -> None:
        self.numStatusChecks += 1
        print(self.aliveString)
        self.StartTimer()

    def StartTimer(self) -> None:
        self.StopTimer()
        self.statusTimer = threading.Timer(float(self.timerPeriod), self._tick)
        self.statusTimer.daemon = True
        self.statusTimer.start()

    start_timer = StartTimer

    def StopTimer(self) -> None:
        if self.statusTimer is not None:
            self.statusTimer.cancel()
            self.statusTimer = None

    stop_timer = StopTimer

    def delete(self) -> None:
        self.StopTimer()

    def __del__(self):
        try:
            self.StopTimer()
        except Exception:
            return None
