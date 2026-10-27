"""Independent reference checks for masked W80 and QC."""

import unittest

import numpy as np

from hic_diagonal_narrowing.width import estimate_w80
from hic_diagonal_narrowing.qc import summarize_local_support


class TestMaskedConsistency(unittest.TestCase):

    def test_against_independent_pair_enumeration(self):
        rng = np.random.default_rng(20261027)

        n = 18
        region = (3, 15)
        max_distance = 5

        for replicate in range(25):
            with self.subTest(replicate=replicate):

                upper = np.triu(
                    rng.poisson(4, size=(n, n)),
                    k=1,
                ).astype(float)

                matrix = upper + upper.T

                valid = rng.random(n) > 0.20

                # Missing genomic bins may contain NaN.
                matrix[~valid, :] = np.nan
                matrix[:, ~valid] = np.nan

                band_values = [[] for _ in range(max_distance)]

                # Independent pair-by-pair reference calculation.
                for i in range(n):
                    for j in range(i + 1, n):
                        d = j - i
                        midpoint = (i + j) / 2

                        if not 1 <= d <= max_distance:
                            continue

                        if not region[0] <= midpoint < region[1]:
                            continue

                        if not (valid[i] and valid[j]):
                            continue

                        band_values[d - 1].append(matrix[i, j])

                expected_pairs = sum(map(len, band_values))

                expected_nonzero = sum(
                    sum(v > 0 for v in band)
                    for band in band_values
                )

                expected_signal = sum(
                    sum(band) for band in band_values
                )

                supported_bands = sum(
                    any(v > 0 for v in band)
                    for band in band_values
                )

                qc = summarize_local_support(
                    matrix,
                    region,
                    max_separation_bins=max_distance,
                    valid_bins=valid,
                )

                self.assertEqual(
                    qc["eligible_pairs"], expected_pairs
                )
                self.assertEqual(
                    qc["nonzero_pairs"], expected_nonzero
                )
                self.assertAlmostEqual(
                    qc["total_contact_signal"], expected_signal
                )
                self.assertEqual(
                    qc["supported_distances"], supported_bands
                )

                if any(len(band) == 0 for band in band_values):
                    with self.assertRaises(ValueError):
                        estimate_w80(
                            matrix,
                            region,
                            max_separation_bins=max_distance,
                            valid_bins=valid,
                        )
                    continue

                means = [
                    sum(band) / len(band)
                    for band in band_values
                ]

                total = sum(means)
                cumulative = 0
                expected_width = None

                for distance, value in enumerate(means, start=1):
                    cumulative += value
                    if cumulative >= 0.8 * total:
                        expected_width = distance
                        break

                observed_width = estimate_w80(
                    matrix,
                    region,
                    max_separation_bins=max_distance,
                    valid_bins=valid,
                )

                self.assertEqual(observed_width, expected_width)


if __name__ == "__main__":
    unittest.main()
