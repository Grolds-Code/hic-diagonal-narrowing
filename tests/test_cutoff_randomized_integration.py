"""Independent randomized validation of genomic W80 uncertainty.

Tests midpoint assignments, masked pairs, distance support,
chromosome boundaries, and matrix/band consistency.
"""

import unittest

import numpy as np

from hic_diagonal_narrowing.width import estimate_w80
from hic_diagonal_narrowing.qc import summarize_local_support
from hic_diagonal_narrowing.poisson_uncertainty import (
    cutoff_exposure_w80_interval,
    cutoff_exposure_w80_interval_from_bands,
)


def independent_band_enumeration(matrix, region, mask, d_max):
    """Enumerate upper-triangle pairs using explicit Python loops."""
    n = matrix.shape[0]
    start, end = region

    counts = [0] * d_max
    exposures = [0] * d_max
    nonzero = 0

    for i in range(n):
        for j in range(i + 1, n):
            distance = j - i

            if distance > d_max:
                continue

            midpoint = (i + j) / 2

            if not start <= midpoint < end:
                continue

            if not (mask[i] and mask[j]):
                continue

            value = int(matrix[i, j])

            counts[distance - 1] += value
            exposures[distance - 1] += 1
            nonzero += int(value > 0)

    return (
        np.asarray(counts, dtype=np.int64),
        np.asarray(exposures, dtype=np.int64),
        nonzero,
    )


class TestRandomizedCutoffIntegration(unittest.TestCase):

    def test_randomized_independent_enumeration(self):
        rng = np.random.default_rng(20261028)

        checked = 0
        missing_support = 0

        n = 32
        windows = (
            (0, 8),
            (8, 24),
            (24, 32),
            (0, 32),
        )

        for trial in range(200):
            d_max = int(rng.integers(2, 9))
            region = windows[trial % len(windows)]

            mask = rng.random(n) > 0.18

            matrix = np.zeros((n, n), dtype=np.int64)

            for i in range(n):
                for j in range(i + 1, n):
                    distance = j - i
                    rate = 1.0 + 7.0 * np.exp(-distance / 3.0)

                    value = int(rng.poisson(rate))

                    matrix[i, j] = value
                    matrix[j, i] = value

            counts, exposures, nonzero = (
                independent_band_enumeration(
                    matrix, region, mask, d_max
                )
            )

            qc = summarize_local_support(
                matrix,
                region,
                max_separation_bins=d_max,
                valid_bins=mask,
            )

            self.assertEqual(
                qc["eligible_pairs"],
                int(exposures.sum()),
            )

            self.assertEqual(
                qc["total_contact_signal"],
                float(counts.sum()),
            )

            self.assertEqual(
                qc["nonzero_pairs"],
                nonzero,
            )

            if np.any(exposures == 0) or counts.sum() == 0:
                missing_support += 1

                with self.assertRaises(ValueError):
                    cutoff_exposure_w80_interval(
                        matrix,
                        region,
                        max_separation_bins=d_max,
                        valid_bins=mask,
                    )
                continue

            # Independently calculate the normalized W80.
            profile = counts / exposures

            expected_w80 = int(
                np.searchsorted(
                    np.cumsum(profile),
                    0.80 * profile.sum(),
                ) + 1
            )

            matrix_result = cutoff_exposure_w80_interval(
                matrix,
                region,
                max_separation_bins=d_max,
                valid_bins=mask,
            )

            band_result = (
                cutoff_exposure_w80_interval_from_bands(
                    counts, exposures
                )
            )

            authoritative_w80 = estimate_w80(
                matrix,
                region,
                max_separation_bins=d_max,
                valid_bins=mask,
            )

            self.assertEqual(
                matrix_result["eligible_pairs_by_distance"],
                exposures.tolist(),
            )

            self.assertEqual(
                matrix_result["contact_counts_by_distance"],
                counts.tolist(),
            )

            self.assertEqual(
                matrix_result["observed_w80"],
                expected_w80,
            )

            self.assertEqual(
                authoritative_w80,
                expected_w80,
            )

            for key in (
                "observed_w80",
                "cutoff_lower",
                "cutoff_upper",
                "total_contacts",
            ):
                self.assertEqual(
                    matrix_result[key],
                    band_result[key],
                )

            self.assertLessEqual(
                matrix_result["cutoff_lower"],
                expected_w80,
            )

            self.assertGreaterEqual(
                matrix_result["cutoff_upper"],
                expected_w80,
            )

            checked += 1

        self.assertGreaterEqual(checked, 100)

        print(
            "\nRandomized validation: "
            f"{checked} valid cases, "
            f"{missing_support} correctly rejected cases"
        )

    def test_explicit_chromosome_boundaries(self):
        n = 20
        matrix = np.zeros((n, n), dtype=np.int64)

        for i in range(n):
            for j in range(i + 1, n):
                value = max(1, 10 - (j - i))
                matrix[i, j] = value
                matrix[j, i] = value

        mask = np.ones(n, dtype=bool)

        for region in ((0, 5), (15, 20)):
            with self.subTest(region=region):
                counts, exposures, _ = (
                    independent_band_enumeration(
                        matrix, region, mask, 4
                    )
                )

                result = cutoff_exposure_w80_interval(
                    matrix,
                    region,
                    max_separation_bins=4,
                    valid_bins=mask,
                )

                self.assertEqual(
                    result["contact_counts_by_distance"],
                    counts.tolist(),
                )

                self.assertEqual(
                    result["eligible_pairs_by_distance"],
                    exposures.tolist(),
                )

    def test_missing_distance_support_is_rejected(self):
        matrix = np.ones((12, 12), dtype=np.int64)
        np.fill_diagonal(matrix, 0)

        mask = np.zeros(12, dtype=bool)
        mask[[0, 1, 2]] = True

        with self.assertRaisesRegex(
            ValueError, "No eligible pairs"
        ):
            cutoff_exposure_w80_interval(
                matrix,
                (0, 5),
                max_separation_bins=4,
                valid_bins=mask,
            )


if __name__ == "__main__":
    unittest.main()
