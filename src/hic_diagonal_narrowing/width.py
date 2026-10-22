"""Initial local Hi-C diagonal-width estimation.

Research prototype v0.1.
Not yet validated for biological interpretation.
"""

import numpy as np


def estimate_w80(
    matrix: np.ndarray,
    region: tuple[int, int],
    *,
    max_separation_bins: int = 64,
) -> int:
    """Estimate local diagonal width in genomic bins.

    Contact pairs are assigned to a region by their midpoint.
    The local contact-distance profile is normalized by the
    number of available bin pairs at each separation.

    W80 is the smallest separation containing at least 80%
    of the accumulated mean contact-distance intensity.

    Assumes a complete, symmetric, nonnegative matrix.
    Missing-bin handling will be added separately.
    """
    matrix = np.asarray(matrix)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    if not np.all(np.isfinite(matrix)):
        raise ValueError("Matrix contains non-finite values")

    if np.any(matrix < 0):
        raise ValueError("Matrix contains negative counts")

    if not np.allclose(matrix, matrix.T):
        raise ValueError("Matrix must be symmetric")

    n_bins = matrix.shape[0]
    start, end = region

    if not (0 <= start < end <= n_bins):
        raise ValueError("Invalid genomic window")

    if not (1 <= max_separation_bins < n_bins):
        raise ValueError("Invalid maximum separation")

    profile = np.zeros(max_separation_bins, dtype=float)

    for distance in range(1, max_separation_bins + 1):
        rows = np.arange(n_bins - distance)
        midpoint = rows + distance / 2

        selected = (
            (midpoint >= start)
            & (midpoint < end)
        )

        if not np.any(selected):
            raise ValueError(
                f"No eligible pairs at separation {distance}"
            )

        contacts = matrix[
            rows[selected],
            rows[selected] + distance,
        ]

        profile[distance - 1] = contacts.mean()

    total_intensity = profile.sum()

    if total_intensity <= 0:
        raise ValueError("No informative contact signal")

    cumulative = np.cumsum(profile)

    w80 = np.searchsorted(
        cumulative,
        0.80 * total_intensity,
    ) + 1

    return int(w80)
