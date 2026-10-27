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


if __name__ == "__main__":
    unittest.main()
