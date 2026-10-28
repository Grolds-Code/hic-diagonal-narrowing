"""Unit tests for unequal-support Poisson W80 intervals."""

import unittest

import numpy as np

from hic_diagonal_narrowing.poisson_uncertainty import (
    poisson_w80_interval_from_bands,
)


class TestPoissonW80Interval(unittest.TestCase):

    def test_equal_support(self):
        result = poisson_w80_interval_from_bands(
            np.array([20, 5, 3]),
            np.array([6, 6, 6]),
        )

        self.assertEqual(result["observed_w80"], 2)
        self.assertLessEqual(result["poisson_lower"], 2)
        self.assertGreaterEqual(result["poisson_upper"], 2)

    def test_unequal_support(self):
        result = poisson_w80_interval_from_bands(
            np.array([100, 40, 20, 10]),
            np.array([10, 5, 8, 4]),
        )

        self.assertEqual(result["observed_w80"], 3)
        self.assertLessEqual(result["poisson_lower"], 3)
        self.assertGreaterEqual(result["poisson_upper"], 3)

    def test_single_contact(self):
        result = poisson_w80_interval_from_bands(
            np.array([1, 0, 0, 0]),
            np.array([5, 3, 2, 1]),
        )

        self.assertEqual(result["poisson_lower"], 1)
        self.assertEqual(result["poisson_upper"], 4)

    def test_higher_confidence_not_narrower(self):
        counts = np.array([100, 40, 20, 10])
        eligible = np.array([10, 5, 8, 4])

        lower_confidence = poisson_w80_interval_from_bands(
            counts, eligible, confidence_level=0.90
        )
        higher_confidence = poisson_w80_interval_from_bands(
            counts, eligible, confidence_level=0.99
        )

        self.assertLessEqual(
            higher_confidence["poisson_lower"],
            lower_confidence["poisson_lower"],
        )
        self.assertGreaterEqual(
            higher_confidence["poisson_upper"],
            lower_confidence["poisson_upper"],
        )

    def test_zero_contacts_rejected(self):
        with self.assertRaises(ValueError):
            poisson_w80_interval_from_bands(
                np.array([0, 0, 0]),
                np.array([3, 3, 3]),
            )

    def test_missing_distance_rejected(self):
        with self.assertRaisesRegex(
            ValueError, "Every separation"
        ):
            poisson_w80_interval_from_bands(
                np.array([5, 0, 2]),
                np.array([3, 0, 3]),
            )

    def test_negative_counts_rejected(self):
        with self.assertRaises(ValueError):
            poisson_w80_interval_from_bands(
                np.array([5, -1]),
                np.array([3, 3]),
            )

    def test_noninteger_counts_rejected(self):
        with self.assertRaises(TypeError):
            poisson_w80_interval_from_bands(
                np.array([5.0, 2.0]),
                np.array([3, 3]),
            )

    def test_invalid_confidence_rejected(self):
        for level in (0, 1, -0.1, 1.1):
            with self.subTest(level=level):
                with self.assertRaises(ValueError):
                    poisson_w80_interval_from_bands(
                        np.array([5, 2]),
                        np.array([3, 3]),
                        confidence_level=level,
                    )


if __name__ == "__main__":
    unittest.main()
