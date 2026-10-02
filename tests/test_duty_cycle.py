import unittest

import numpy as np

from minke.duty_cycle import (
    active_detectors,
    generate_duty_cycle_schedule,
)

_YEAR = 365.25 * 86400.0


class TestGenerateDutyCycleSchedule(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(1234)
        self.t_start = 1_238_166_018.0
        self.t_end = self.t_start + _YEAR

    def test_segments_are_sorted_and_contiguous(self):
        schedule = generate_duty_cycle_schedule(
            duty_cycle=0.75,
            mean_lock_duration=8 * 3600.0,
            t_start=self.t_start,
            t_end=self.t_end,
            rng=self.rng,
        )
        segments = schedule.segments
        self.assertEqual(segments[0].start, self.t_start)
        self.assertEqual(segments[-1].end, self.t_end)
        for a, b in zip(segments[:-1], segments[1:]):
            self.assertEqual(a.end, b.start)
            # Alternating state
            self.assertNotEqual(a.active, b.active)

    def test_realised_duty_cycle_converges_to_target(self):
        # Long window + many segments -> the realised locked fraction
        # should be close to the target duty cycle.
        schedule = generate_duty_cycle_schedule(
            duty_cycle=0.75,
            mean_lock_duration=8 * 3600.0,
            t_start=self.t_start,
            t_end=self.t_start + 50 * _YEAR,
            rng=self.rng,
        )
        locked_time = sum(
            segment.end - segment.start
            for segment in schedule.segments
            if segment.active
        )
        total_time = schedule.segments[-1].end - schedule.segments[0].start
        self.assertAlmostEqual(locked_time / total_time, 0.75, delta=0.02)

    def test_reproducible_with_same_seed(self):
        schedule_a = generate_duty_cycle_schedule(
            duty_cycle=0.7,
            mean_lock_duration=3600.0,
            t_start=0.0,
            t_end=1e6,
            rng=np.random.default_rng(42),
        )
        schedule_b = generate_duty_cycle_schedule(
            duty_cycle=0.7,
            mean_lock_duration=3600.0,
            t_start=0.0,
            t_end=1e6,
            rng=np.random.default_rng(42),
        )
        self.assertEqual(schedule_a.segments, schedule_b.segments)

    def test_is_active_matches_segments(self):
        schedule = generate_duty_cycle_schedule(
            duty_cycle=0.6,
            mean_lock_duration=1000.0,
            t_start=0.0,
            t_end=100_000.0,
            rng=self.rng,
        )
        for segment in schedule.segments:
            midpoint = (segment.start + segment.end) / 2.0
            self.assertEqual(schedule.is_active(midpoint), segment.active)

    def test_is_active_out_of_range_raises(self):
        schedule = generate_duty_cycle_schedule(
            duty_cycle=0.5,
            mean_lock_duration=1000.0,
            t_start=0.0,
            t_end=10_000.0,
            rng=self.rng,
        )
        with self.assertRaises(ValueError):
            schedule.is_active(-1.0)
        with self.assertRaises(ValueError):
            schedule.is_active(10_001.0)

    def test_non_finite_inputs_are_rejected(self):
        base = dict(
            duty_cycle=0.5, mean_lock_duration=1000.0, t_start=0.0, t_end=10_000.0
        )
        for key in base:
            for bad in (float("nan"), float("inf"), float("-inf")):
                with self.subTest(key=key, value=bad):
                    with self.assertRaises(ValueError):
                        generate_duty_cycle_schedule(
                            **{**base, key: bad}, rng=self.rng
                        )

    def test_starts_are_cached_between_queries(self):
        schedule = generate_duty_cycle_schedule(
            duty_cycle=0.5,
            mean_lock_duration=1000.0,
            t_start=0.0,
            t_end=10_000.0,
            rng=self.rng,
        )
        schedule.is_active(5.0)
        self.assertIs(schedule._starts, schedule._starts)

    def test_duty_cycle_of_one_is_always_active(self):
        schedule = generate_duty_cycle_schedule(
            duty_cycle=1.0,
            mean_lock_duration=1000.0,
            t_start=0.0,
            t_end=10_000.0,
            rng=self.rng,
        )
        self.assertTrue(schedule.is_active(5000.0))

    def test_duty_cycle_of_zero_is_never_active(self):
        schedule = generate_duty_cycle_schedule(
            duty_cycle=0.0,
            mean_lock_duration=1000.0,
            t_start=0.0,
            t_end=10_000.0,
            rng=self.rng,
        )
        self.assertFalse(schedule.is_active(5000.0))

    def test_invalid_duty_cycle_raises(self):
        with self.assertRaises(ValueError):
            generate_duty_cycle_schedule(
                duty_cycle=1.5,
                mean_lock_duration=1000.0,
                t_start=0.0,
                t_end=10_000.0,
                rng=self.rng,
            )

    def test_invalid_window_raises(self):
        with self.assertRaises(ValueError):
            generate_duty_cycle_schedule(
                duty_cycle=0.5,
                mean_lock_duration=1000.0,
                t_start=10_000.0,
                t_end=0.0,
                rng=self.rng,
            )


class TestActiveDetectors(unittest.TestCase):
    def test_restricts_to_active_schedules(self):
        rng = np.random.default_rng(7)
        t_start, t_end = 0.0, 100_000.0
        schedules = {
            "AlwaysOn": generate_duty_cycle_schedule(1.0, 1000.0, t_start, t_end, rng),
            "AlwaysOff": generate_duty_cycle_schedule(0.0, 1000.0, t_start, t_end, rng),
        }
        detectors = {"AlwaysOn": "SomePSD", "AlwaysOff": "SomePSD"}
        result = active_detectors(schedules, detectors, t=50_000.0)
        self.assertEqual(result, {"AlwaysOn": "SomePSD"})


if __name__ == "__main__":
    unittest.main()
