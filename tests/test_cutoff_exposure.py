"""Tests for deterministic cutoff-specific W80 intervals."""

import unittest

import numpy as np
from scipy.stats import binom, multinomial

from hic_diagonal_narrowing.poisson_uncertainty import (
    cutoff_exposure_w80_interval_from_bands,
)


class TestCutoffExposureW80(unittest.TestCase):

    def test_equal_support_matches_rank_quantiles(self):
        counts = np.array([20, 5, 3, 2], dtype=np.int64)
        exposure = np.full(4, 6, dtype=np.int64)

        result = cutoff_exposure_w80_interval_from_bands(
            counts, exposure
        )

        n = int(counts.sum())
        cumulative = np.cumsum(counts)

        lower_rank = int(binom.ppf(0.025, n, 0.80))
        upper_rank = int(binom.isf(0.025, n, 0.80)) + 1

        expected_lower = (
            1 if lower_rank < 1
            else int(
                np.searchsorted(cumulative, lower_rank) + 1
            )
        )

        expected_upper = (
            4 if upper_rank > n
            else int(
                np.searchsorted(cumulative, upper_rank) + 1
            )
        )

        self.assertEqual(
            result["cutoff_lower"], expected_lower
        )
        self.assertEqual(
            result["cutoff_upper"], expected_upper
        )

    def test_unequal_exposure(self):
        counts = np.array([100, 40, 20, 10])
        exposure = np.array([10, 5, 8, 4])

        result = cutoff_exposure_w80_interval_from_bands(
            counts, exposure
        )

        self.assertEqual(result["observed_w80"], 3)
        self.assertLessEqual(result["cutoff_lower"], 3)
        self.assertGreaterEqual(result["cutoff_upper"], 3)

    def test_single_contact(self):
        result = cutoff_exposure_w80_interval_from_bands(
            np.array([1, 0, 0, 0]),
            np.array([10, 3, 5, 2]),
        )

        self.assertEqual(result["cutoff_lower"], 1)
        self.assertEqual(result["cutoff_upper"], 4)

    def test_deterministic(self):
        counts = np.array([50, 20, 9, 3])
        exposure = np.array([10, 7, 3, 1])

        first = cutoff_exposure_w80_interval_from_bands(
            counts, exposure
        )

        second = cutoff_exposure_w80_interval_from_bands(
            counts, exposure
        )

        self.assertEqual(first, second)

    def test_confidence_nesting(self):
        counts = np.array([100, 40, 20, 10])
        exposure = np.array([10, 5, 8, 4])

        low = cutoff_exposure_w80_interval_from_bands(
            counts, exposure, confidence_level=0.90
        )

        high = cutoff_exposure_w80_interval_from_bands(
            counts, exposure, confidence_level=0.99
        )

        self.assertLessEqual(
            high["cutoff_lower"], low["cutoff_lower"]
        )
        self.assertGreaterEqual(
            high["cutoff_upper"], low["cutoff_upper"]
        )

    def test_zero_contacts_rejected(self):
        with self.assertRaises(ValueError):
            cutoff_exposure_w80_interval_from_bands(
                np.zeros(4, dtype=int),
                np.ones(4, dtype=int),
            )

    def test_zero_exposure_rejected(self):
        with self.assertRaises(ValueError):
            cutoff_exposure_w80_interval_from_bands(
                np.array([5, 3, 2]),
                np.array([4, 0, 3]),
            )

    def test_noninteger_counts_rejected(self):
        with self.assertRaises(TypeError):
            cutoff_exposure_w80_interval_from_bands(
                np.array([5.0, 3.0]),
                np.array([4, 4]),
            )

    def test_invalid_confidence_rejected(self):
        for level in (0, 1, -0.1, 1.1):
            with self.subTest(level=level):
                with self.assertRaises(ValueError):
                    cutoff_exposure_w80_interval_from_bands(
                        np.array([5, 3]),
                        np.array([4, 4]),
                        confidence_level=level,
                    )

    def test_exact_small_sample_coverage(self):
        # Independently enumerate the full multinomial
        # sampling distribution for two four-band models.
        models = [
            (
                np.array([4, 4, 4, 4]),
                np.array([4., 1., 3., 2.]),
            ),
            (
                np.array([1, 10, 2, 8]),
                np.array([4., 1., 3., 2.]),
            ),
        ]

        n = 8

        for exposure, intensity in models:
            truth = int(
                np.searchsorted(
                    np.cumsum(intensity),
                    0.8 * intensity.sum(),
                ) + 1
            )

            probabilities = exposure * intensity
            probabilities = (
                probabilities / probabilities.sum()
            )

            coverage = 0.0

            for a in range(n + 1):
                for b in range(n - a + 1):
                    for c in range(n - a - b + 1):
                        counts = np.array([
                            a, b, c, n - a - b - c
                        ])

                        result = (
                            cutoff_exposure_w80_interval_from_bands(
                                counts, exposure
                            )
                        )

                        mass = multinomial.pmf(
                            counts, n=n, p=probabilities
                        )

                        if (
                            result["cutoff_lower"]
                            <= truth
                            <= result["cutoff_upper"]
                        ):
                            coverage += mass

            self.assertGreaterEqual(
                coverage, 0.95 - 1e-10
            )


if __name__ == "__main__":
    unittest.main()
