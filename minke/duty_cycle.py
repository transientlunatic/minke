"""Detector duty-cycle simulation.

Real gravitational-wave detectors are not continuously operating: each one
alternates between "locked" (observing) and "unlocked" (down) stretches.
This module models that as a two-state renewal process -- a detector
alternates locked/unlocked segments with exponentially-distributed
durations -- so that which detectors are "on" varies realistically (and
with time-correlation, not independently event-by-event) across an
observation window.
"""

from __future__ import annotations

import bisect
import functools
import math
from dataclasses import dataclass

import numpy as np

__all__ = [
    "DutyCycleSegment",
    "DutyCycleSchedule",
    "generate_duty_cycle_schedule",
    "active_detectors",
]


@dataclass(frozen=True)
class DutyCycleSegment:
    """One contiguous locked/unlocked interval.

    Parameters
    ----------
    start, end : float
        GPS start/end time of the segment, in seconds.
    active : bool
        ``True`` if the detector is locked (observing) during this segment.
    """

    start: float
    end: float
    active: bool


@dataclass(frozen=True)
class DutyCycleSchedule:
    """A detector's locked/unlocked segments over an observation window.

    ``segments`` is sorted by start time, contiguous (each segment's end
    equals the next segment's start), and covers the full window passed to
    :func:`generate_duty_cycle_schedule`.
    """

    segments: tuple[DutyCycleSegment, ...]

    @functools.cached_property
    def _starts(self) -> list[float]:
        # Built once per schedule so each ``is_active`` query is O(log n).
        return [segment.start for segment in self.segments]

    def is_active(self, t: float) -> bool:
        """Return whether the detector is locked at GPS time *t*.

        Parameters
        ----------
        t : float
            GPS time in seconds. Must fall within the schedule's window
            (``segments[0].start <= t <= segments[-1].end``).

        Raises
        ------
        ValueError
            If *t* falls outside the schedule's window.
        """
        starts = self._starts
        if t < starts[0] or t > self.segments[-1].end:
            raise ValueError(
                f"t={t!r} is outside the schedule's window "
                f"[{starts[0]!r}, {self.segments[-1].end!r}]"
            )
        # bisect_right(starts, t) - 1 is the index of the last segment whose
        # start is <= t; segments are contiguous, so that's the segment
        # containing t (or, at t == segments[-1].end exactly, the last one).
        index = min(bisect.bisect_right(starts, t) - 1, len(self.segments) - 1)
        return self.segments[index].active


def generate_duty_cycle_schedule(
    duty_cycle: float,
    mean_lock_duration: float,
    t_start: float,
    t_end: float,
    rng: np.random.Generator,
) -> DutyCycleSchedule:
    r"""Generate a two-state renewal-process duty-cycle schedule.

    The detector alternates locked/unlocked segments with
    exponentially-distributed durations. The mean locked duration is
    *mean_lock_duration*; the mean unlocked duration is derived so the
    long-run locked fraction equals *duty_cycle*:

    .. math::

        \overline{T}_{\text{unlocked}} = \overline{T}_{\text{locked}}
        \frac{1 - d}{d}

    where :math:`d` is *duty_cycle*. The initial state at *t_start* is
    drawn as locked with probability *duty_cycle*.

    Parameters
    ----------
    duty_cycle : float
        Target long-run fraction of time locked, in ``(0, 1)``. ``1.0``
        gives a single locked segment covering the whole window; ``0.0``
        gives a single unlocked segment.
    mean_lock_duration : float
        Mean duration in seconds of a locked segment.
    t_start, t_end : float
        GPS start/end of the observation window in seconds.
    rng : numpy.random.Generator
        Random number generator, for reproducibility.

    Returns
    -------
    DutyCycleSchedule
        Segments covering ``[t_start, t_end]``.

    Raises
    ------
    ValueError
        If ``duty_cycle`` is outside ``[0, 1]``, ``mean_lock_duration`` is
        not positive, ``t_end <= t_start``, or any argument is non-finite.
    """
    for label, value in (
        ("duty_cycle", duty_cycle),
        ("mean_lock_duration", mean_lock_duration),
        ("t_start", t_start),
        ("t_end", t_end),
    ):
        if not math.isfinite(value):
            raise ValueError(f"{label} must be finite, got {value!r}")
    if not (0.0 <= duty_cycle <= 1.0):
        raise ValueError(f"duty_cycle must be in [0, 1], got {duty_cycle!r}")
    if mean_lock_duration <= 0.0:
        raise ValueError(
            f"mean_lock_duration must be positive, got {mean_lock_duration!r}"
        )
    if t_end <= t_start:
        raise ValueError(f"t_end ({t_end!r}) must be > t_start ({t_start!r})")

    if duty_cycle in (0.0, 1.0):
        return DutyCycleSchedule(
            segments=(DutyCycleSegment(t_start, t_end, bool(duty_cycle)),)
        )

    mean_unlock_duration = mean_lock_duration * (1.0 - duty_cycle) / duty_cycle

    active = bool(rng.random() < duty_cycle)
    segments: list[DutyCycleSegment] = []
    t = t_start
    while t < t_end:
        mean_duration = mean_lock_duration if active else mean_unlock_duration
        duration = rng.exponential(mean_duration)
        end = min(t + duration, t_end)
        segments.append(DutyCycleSegment(start=t, end=end, active=active))
        t = end
        active = not active

    return DutyCycleSchedule(segments=tuple(segments))


def active_detectors(
    schedules: dict[str, DutyCycleSchedule],
    detectors: dict[str, str],
    t: float,
) -> dict[str, str]:
    """Restrict a detector network to those active (locked) at GPS time *t*.

    Parameters
    ----------
    schedules : dict[str, DutyCycleSchedule]
        One schedule per detector name, as generated by
        :func:`generate_duty_cycle_schedule`. Must have an entry for every
        key in *detectors*.
    detectors : dict[str, str]
        Detector network, mapping detector name to PSD name (as passed to
        :func:`minke.injection.make_injection`).
    t : float
        GPS time in seconds to query.

    Returns
    -------
    dict[str, str]
        The subset of *detectors* whose schedule is active at *t*, in the
        same ``{name: psd_name}`` form.
    """
    return {
        name: psd_name
        for name, psd_name in detectors.items()
        if schedules[name].is_active(t)
    }
