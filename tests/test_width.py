import unittest
import numpy as np

from hic_diagonal_narrowing.width import estimate_w80


class TestW80(unittest.TestCase):
    def test_adjacent_contacts(self):
        matrix = np.zeros((8, 8), dtype=int)

        for i in range(7):
            matrix[i, i + 1] = 2
            matrix[i + 1, i] = 2

        result = estimate_w80(
            matrix,
            (0, 8),
            max_separation_bins=3,
        )

        self.assertEqual(result, 1)



    def test_contacts_at_four_bins(self):
        matrix = np.zeros((12, 12), dtype=int)

        for i in range(8):
            matrix[i, i + 4] = 3
            matrix[i + 4, i] = 3

        result = estimate_w80(
            matrix,
            (0, 12),
            max_separation_bins=6,
        )

        self.assertEqual(result, 4)



    def test_midpoint_window_assignment(self):
        matrix = np.zeros((12, 12), dtype=int)

        matrix[5, 7] = matrix[7, 5] = 10
        matrix[0, 1] = matrix[1, 0] = 1000

        result = estimate_w80(
            matrix,
            (4, 8),
            max_separation_bins=4,
        )

        self.assertEqual(result, 2)



    def test_distance_normalization(self):
        matrix = np.zeros((12, 12), dtype=int)

        for i in range(11):
            matrix[i, i + 1] = matrix[i + 1, i] = 3

        for i in range(8):
            matrix[i, i + 4] = matrix[i + 4, i] = 1

        result = estimate_w80(
            matrix, (0, 12), max_separation_bins=4
        )

        self.assertEqual(result, 4)



    def test_w80_excludes_invalid_bin_pairs(self):
        """Masked contacts must not contribute to W80."""
        matrix = np.zeros((6, 6), dtype=float)

        matrix[0, 1] = matrix[1, 0] = 60
        matrix[1, 3] = matrix[3, 1] = 10

        # The strong adjacent contact makes unmasked W80 = 1.
        unmasked = estimate_w80(
            matrix,
            (0, 6),
            max_separation_bins=2,
        )
        self.assertEqual(unmasked, 1)

        # Bin 0 is invalid. Its contacts must be excluded.
        valid_bins = np.array(
            [False, True, True, True, True, True]
        )

        masked = estimate_w80(
            matrix,
            (0, 6),
            max_separation_bins=2,
            valid_bins=valid_bins,
        )

        self.assertEqual(masked, 2)



    def test_missing_distance_band_is_not_zero_signal(self):
        """Distinguish observed zero contacts from missing support."""
        matrix = np.zeros((6, 6), dtype=float)

        matrix[0, 1] = matrix[1, 0] = 20

        # Separation 2 has eligible pairs but zero contacts.
        # This is a valid observation.
        result = estimate_w80(
            matrix,
            (0, 6),
            max_separation_bins=2,
        )
        self.assertEqual(result, 1)

        # Only bins 0 and 1 remain valid.
        # Separation 2 now has no eligible pairs at all.
        valid_bins = np.array(
            [True, True, False, False, False, False]
        )

        with self.assertRaisesRegex(
            ValueError,
            "No eligible pairs at separation 2",
        ):
            estimate_w80(
                matrix,
                (0, 6),
                max_separation_bins=2,
                valid_bins=valid_bins,
            )


if __name__ == "__main__":
    unittest.main()
