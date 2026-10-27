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
    valid_bins: np.ndarray | None = None,
) -> int:
    """Estimate local diagonal width in genomic bins.

    Contact pairs are assigned to a region by their midpoint.
    The local contact-distance profile is normalized by the
    number of available bin pairs at each separation.

    W80 is the smallest separation containing at least 80%
    of the accumulated mean contact-distance intensity.

    Supports an optional Boolean genomic-bin validity mask.
    Pairs touching invalid bins are excluded from both sums
    and eligible-pair denominators.

    Every separation from 1 to max_separation_bins must retain
    at least one eligible pair; otherwise W80 is undefined
    under this complete-distance-support policy.
    """
    matrix = np.asarray(matrix)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    n_bins = matrix.shape[0]

    if valid_bins is None:
        valid_bins = np.ones(n_bins, dtype=bool)
    else:
        valid_bins = np.asarray(valid_bins)

        if valid_bins.shape != (n_bins,):
            raise ValueError(
                "valid_bins must contain one value per genomic bin"
            )

        if valid_bins.dtype != np.dtype(bool):
            raise ValueError("valid_bins must be a boolean array")

    valid_idx = np.flatnonzero(valid_bins)
    valid_matrix = matrix[np.ix_(valid_idx, valid_idx)]

    if not np.all(np.isfinite(valid_matrix)):
        raise ValueError(
            "Valid genomic bins contain non-finite values"
        )

    if np.any(valid_matrix < 0):
        raise ValueError(
            "Valid genomic bins contain negative counts"
        )

    if not np.allclose(valid_matrix, valid_matrix.T):
        raise ValueError(
            "Contacts between valid bins must be symmetric"
        )

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
            & valid_bins[rows]
            & valid_bins[rows + distance]
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
