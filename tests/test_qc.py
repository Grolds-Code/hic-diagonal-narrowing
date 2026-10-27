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


if __name__ == "__main__":
    unittest.main()
