"""Verify mask-aware genomic integration of W80 uncertainty."""

import unittest

import numpy as np

from hic_diagonal_narrowing.width import estimate_w80
from hic_diagonal_narrowing.qc import summarize_local_support
from hic_diagonal_narrowing.poisson_uncertainty import (
    cutoff_exposure_w80_interval,
    cutoff_exposure_w80_interval_from_bands,
)


class TestCutoffMatrixIntegration(unittest.TestCase):

    def setUp(self):
        self.matrix = np.zeros((20, 20), dtype=np.int64)
        self.region = (5, 15)
        self.d_max = 4

        for distance in range(1, self.d_max + 1):
            for row in range(20 - distance):
                col = row + distance

                value = 6 - distance
                self.matrix[row, col] = value
                self.matrix[col, row] = value

        self.mask = np.ones(20, dtype=bool)
        self.mask[[6, 9, 12]] = False

    def test_masked_matrix_matches_width_and_qc(self):
        result = cutoff_exposure_w80_interval(
            self.matrix,
            self.region,
            max_separation_bins=self.d_max,
            valid_bins=self.mask,
        )

        width = estimate_w80(
            self.matrix,
            self.region,
            max_separation_bins=self.d_max,
            valid_bins=self.mask,
        )

        qc = summarize_local_support(
            self.matrix,
            self.region,
            max_separation_bins=self.d_max,
            valid_bins=self.mask,
        )

        self.assertEqual(result["observed_w80"], width)

        self.assertEqual(
            result["eligible_pairs_by_distance"],
            [4, 4, 7, 5],
        )

        self.assertEqual(
            sum(result["eligible_pairs_by_distance"]),
            qc["eligible_pairs"],
        )

        self.assertEqual(
            sum(result["contact_counts_by_distance"]),
            int(qc["total_contact_signal"]),
        )

        self.assertLessEqual(
            result["cutoff_lower"], width
        )
        self.assertGreaterEqual(
            result["cutoff_upper"], width
        )

    def test_matrix_and_band_interfaces_agree(self):
        matrix_result = cutoff_exposure_w80_interval(
            self.matrix,
            self.region,
            max_separation_bins=self.d_max,
            valid_bins=self.mask,
        )

        band_result = (
            cutoff_exposure_w80_interval_from_bands(
                np.asarray(
                    matrix_result["contact_counts_by_distance"],
                    dtype=np.int64,
                ),
                np.asarray(
                    matrix_result["eligible_pairs_by_distance"],
                    dtype=np.int64,
                ),
            )
        )

        for key in (
            "observed_w80",
            "cutoff_lower",
            "cutoff_upper",
            "total_contacts",
        ):
            self.assertEqual(
                matrix_result[key], band_result[key]
            )

    def test_masked_bins_do_not_affect_results(self):
        altered = self.matrix.copy()

        # Values involving invalid bins must be excluded.
        altered[6, :] = 100000
        altered[:, 6] = 100000

        original = cutoff_exposure_w80_interval(
            self.matrix,
            self.region,
            max_separation_bins=self.d_max,
            valid_bins=self.mask,
        )

        modified = cutoff_exposure_w80_interval(
            altered,
            self.region,
            max_separation_bins=self.d_max,
            valid_bins=self.mask,
        )

        self.assertEqual(original, modified)

    def test_rejects_noninteger_matrix(self):
        with self.assertRaises(TypeError):
            cutoff_exposure_w80_interval(
                self.matrix.astype(float),
                self.region,
                max_separation_bins=self.d_max,
                valid_bins=self.mask,
            )


if __name__ == "__main__":
    unittest.main()
