"""Unit tests for the exact midpoint-assigned Stage 2 preflight regression."""

import unittest

import numpy as np

from experiments.stage2_pooled_preflight_core import midpoint_w80_profile


class Stage2PooledPreflightTests(unittest.TestCase):
    def test_midpoint_profile_matches_naive_window_enumeration(self):
        x = np.array([0, 0, 10, 20, 30, 40, 50], dtype=np.int64)
        y = np.array([10, 20, 20, 40, 40, 60, 50], dtype=np.int64)
        values = np.array([2, 3, 4, 5, 6, 7, 9], dtype=np.int64)

        width, totals, n_records = midpoint_w80_profile(
            x, y, values,
            start=0,
            end=80,
            resolution=10,
            window_bins=3,
            max_distance_bins=3,
        )
        self.assertEqual(n_records, 6)
        self.assertEqual(len(width), 6)

        expected_width = []
        expected_totals = []
        for start_bin in range(6):
            bands = np.zeros(3, dtype=np.int64)
            for x_bp, y_bp, count in zip(x, y, values):
                distance = (y_bp - x_bp) // 10
                midpoint = (x_bp + y_bp) // 20
                if (
                    1 <= distance <= 3
                    and start_bin <= midpoint < start_bin + 3
                ):
                    bands[distance - 1] += count
            total = int(bands.sum())
            expected_totals.append(total)
            if total == 0:
                expected_width.append(-1)
            else:
                expected_width.append(
                    int(np.searchsorted(
                        np.cumsum(bands),
                        0.8 * total,
                    ) + 1)
                )

        np.testing.assert_array_equal(width, expected_width)
        np.testing.assert_array_equal(totals, expected_totals)

    def test_no_off_diagonal_signal_is_explicitly_unscored(self):
        width, totals, count = midpoint_w80_profile(
            np.array([0], dtype=np.int64),
            np.array([0], dtype=np.int64),
            np.array([5], dtype=np.int64),
            start=0, end=80,
            resolution=10, window_bins=3, max_distance_bins=3,
        )
        self.assertEqual(count, 0)
        self.assertTrue(np.all(width == -1))
        self.assertTrue(np.all(totals == 0))

    def test_invalid_records_rejected(self):
        for y, counts in (([15], [1]), ([10], [-1]), ([10], [1.5])):
            with self.subTest(y=y, counts=counts):
                with self.assertRaises(ValueError):
                    midpoint_w80_profile(
                        np.array([0]),
                        np.array(y),
                        np.array(counts),
                        start=0, end=80,
                        resolution=10, window_bins=3,
                        max_distance_bins=3,
                    )


if __name__ == "__main__":
    unittest.main()
