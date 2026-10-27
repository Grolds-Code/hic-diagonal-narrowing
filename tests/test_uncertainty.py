"""Tests for rank-based W80 uncertainty estimation."""

import unittest

import numpy as np

from hic_diagonal_narrowing.uncertainty import (
    rank_interval_w80_equal_support,
)


class TestRankIntervalW80(unittest.TestCase):

    def setUp(self):
        self.matrix = np.zeros((12, 12), dtype=np.int64)

        for distance, count in ((1, 20), (2, 5), (3, 3)):
            for row in range(12 - distance):
                midpoint = row + distance / 2

                if 3 <= midpoint < 9:
                    col = row + distance
                    self.matrix[row, col] = count
                    self.matrix[col, row] = count

    def test_returns_valid_interval(self):
        result = rank_interval_w80_equal_support(
            self.matrix,
            (3, 9),
            max_separation_bins=3,
        )

        self.assertEqual(result["observed_w80"], 2)
        self.assertEqual(
            result["eligible_pairs_per_distance"], 6
        )
        self.assertLessEqual(
            result["rank_lower"], result["observed_w80"]
        )
        self.assertGreaterEqual(
            result["rank_upper"], result["observed_w80"]
        )

    def test_single_contact_produces_full_range(self):
        matrix = np.zeros((12, 12), dtype=np.int64)
        matrix[4, 5] = 1
        matrix[5, 4] = 1

        result = rank_interval_w80_equal_support(
            matrix,
            (3, 9),
            max_separation_bins=3,
        )

        self.assertEqual(result["rank_lower"], 1)
        self.assertEqual(result["rank_upper"], 3)
        self.assertEqual(result["total_contacts"], 1)

    def test_rejects_unequal_distance_support(self):
        matrix = self.matrix.copy()

        # Provide genuine signal within the selected window.
        # Eligible-pair counts still differ across distances.
        matrix[0, 1] = 1
        matrix[1, 0] = 1

        with self.assertRaisesRegex(
            ValueError, "equal eligible-pair counts"
        ):
            rank_interval_w80_equal_support(
                matrix,
                (0, 3),
                max_separation_bins=3,
            )

    def test_rejects_zero_contact_signal(self):
        matrix = np.zeros((12, 12), dtype=np.int64)

        with self.assertRaises(ValueError):
            rank_interval_w80_equal_support(
                matrix,
                (3, 9),
                max_separation_bins=3,
            )

    def test_rejects_noninteger_counts(self):
        with self.assertRaises(TypeError):
            rank_interval_w80_equal_support(
                self.matrix.astype(float),
                (3, 9),
                max_separation_bins=3,
            )

    def test_rejects_invalid_confidence_level(self):
        for confidence_level in (0, 1, -0.1, 1.1):
            with self.subTest(confidence_level=confidence_level):
                with self.assertRaises(ValueError):
                    rank_interval_w80_equal_support(
                        self.matrix,
                        (3, 9),
                        max_separation_bins=3,
                        confidence_level=confidence_level,
                    )


if __name__ == "__main__":
    unittest.main()
