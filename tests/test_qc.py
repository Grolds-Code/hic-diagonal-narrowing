import unittest
import numpy as np

from hic_diagonal_narrowing.qc import summarize_local_support


class TestQualityControl(unittest.TestCase):
    def test_contact_support(self):
        matrix = np.zeros((6, 6), dtype=int)

        for i, j, count in [(0, 1, 2), (2, 3, 3), (1, 3, 4)]:
            matrix[i, j] = count
            matrix[j, i] = count

        qc = summarize_local_support(
            matrix, (0, 6), max_separation_bins=2
        )

        self.assertEqual(qc["eligible_pairs"], 9)
        self.assertEqual(qc["nonzero_pairs"], 3)
        self.assertEqual(qc["total_contact_signal"], 9)
        self.assertAlmostEqual(qc["pair_occupancy"], 1 / 3)
        self.assertEqual(qc["supported_distances"], 2)
        self.assertEqual(qc["eligible_distances"], 2)



    def test_explicit_valid_bin_mask(self):
        """Pairs touching invalid genomic bins must be excluded."""
        matrix = np.zeros((6, 6), dtype=int)

        for i, j, count in [(0, 1, 2), (2, 3, 3), (1, 3, 4)]:
            matrix[i, j] = count
            matrix[j, i] = count

        valid_bins = np.array(
            [True, True, False, True, True, True]
        )

        qc = summarize_local_support(
            matrix,
            (0, 6),
            max_separation_bins=2,
            valid_bins=valid_bins,
        )

        self.assertEqual(qc["eligible_pairs"], 5)
        self.assertEqual(qc["nonzero_pairs"], 2)
        self.assertEqual(qc["total_contact_signal"], 6)
        self.assertAlmostEqual(qc["pair_occupancy"], 2 / 5)



    def test_nan_values_in_masked_bins(self):
        """Missing values in excluded bins must not affect QC."""
        matrix = np.zeros((6, 6), dtype=float)

        for i, j, count in [(0, 1, 2), (1, 3, 4)]:
            matrix[i, j] = count
            matrix[j, i] = count

        matrix[2, :] = np.nan
        matrix[:, 2] = np.nan

        valid_bins = np.array(
            [True, True, False, True, True, True]
        )

        qc = summarize_local_support(
            matrix,
            (0, 6),
            max_separation_bins=2,
            valid_bins=valid_bins,
        )

        self.assertEqual(qc["eligible_pairs"], 5)
        self.assertEqual(qc["nonzero_pairs"], 2)
        self.assertEqual(qc["total_contact_signal"], 6)
        self.assertAlmostEqual(qc["pair_occupancy"], 2 / 5)



    def test_reject_invalid_contacts_between_valid_bins(self):
        """Masking must not conceal invalid data in valid bins."""
        valid_bins = np.array(
            [True, True, False, True, True, True]
        )

        cases = [
            ("nonfinite", np.nan, np.nan, "non-finite"),
            ("negative", -2.0, -2.0, "negative"),
            ("asymmetric", 2.0, 3.0, "symmetric"),
        ]

        for name, forward, reverse, message in cases:
            with self.subTest(case=name):
                matrix = np.zeros((6, 6), dtype=float)
                matrix[0, 1] = forward
                matrix[1, 0] = reverse

                with self.assertRaisesRegex(ValueError, message):
                    summarize_local_support(
                        matrix,
                        (0, 6),
                        max_separation_bins=2,
                        valid_bins=valid_bins,
                    )



    def test_reject_invalid_bin_masks(self):
        """Reject incorrectly shaped or non-Boolean genomic masks."""
        matrix = np.zeros((6, 6), dtype=float)

        cases = [
            (
                "wrong_length",
                np.ones(5, dtype=bool),
                "one value per genomic bin",
            ),
            (
                "wrong_dimensions",
                np.ones((6, 1), dtype=bool),
                "one value per genomic bin",
            ),
            (
                "integer_mask",
                np.ones(6, dtype=int),
                "boolean array",
            ),
        ]

        for name, mask, expected_error in cases:
            with self.subTest(case=name):
                with self.assertRaisesRegex(
                    ValueError, expected_error
                ):
                    summarize_local_support(
                        matrix,
                        (0, 6),
                        max_separation_bins=2,
                        valid_bins=mask,
                    )



    def test_midpoint_and_chromosome_boundaries(self):
        """Use half-open midpoint windows at chromosome edges."""
        matrix = np.zeros((6, 6), dtype=float)

        for i, j, count in [
            (0, 1, 2),
            (0, 2, 7),
            (1, 2, 3),
            (4, 5, 5),
        ]:
            matrix[i, j] = count
            matrix[j, i] = count

        first = summarize_local_support(
            matrix, (0, 1), max_separation_bins=2
        )

        middle = summarize_local_support(
            matrix, (1, 2), max_separation_bins=2
        )

        last = summarize_local_support(
            matrix, (5, 6), max_separation_bins=2
        )

        # Midpoint 0.5 belongs to [0, 1).
        self.assertEqual(first["eligible_pairs"], 1)
        self.assertEqual(first["total_contact_signal"], 2)

        # Midpoints 1.0 and 1.5 belong to [1, 2).
        self.assertEqual(middle["eligible_pairs"], 2)
        self.assertEqual(middle["total_contact_signal"], 10)
        self.assertEqual(middle["supported_distances"], 2)

        # No valid pair midpoint lies in [5, 6).
        self.assertEqual(last["eligible_pairs"], 0)
        self.assertEqual(last["total_contact_signal"], 0)
        self.assertEqual(last["pair_occupancy"], 0)
        self.assertEqual(last["distance_support_fraction"], 0)


if __name__ == "__main__":
    unittest.main()
