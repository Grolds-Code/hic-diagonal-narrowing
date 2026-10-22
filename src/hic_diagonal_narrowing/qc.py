"""Coverage diagnostics for local Hi-C diagonal analysis.

Research prototype v0.1.

These metrics describe observed contact support.
They are not validated biological confidence measures.
"""

import numpy as np


def summarize_local_support(
    matrix: np.ndarray,
    region: tuple[int, int],
    *,
    max_separation_bins: int = 64,
) -> dict:
    """Summarize contact support for a genomic window.

    Contacts are assigned to the region by their midpoint.
    Only unique upper-triangle, off-diagonal pairs are used.

    This version assumes a complete matrix without explicit
    missing-bin masks.
    """
    matrix = np.asarray(matrix)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    if not np.all(np.isfinite(matrix)):
        raise ValueError("Matrix contains non-finite values")

    if np.any(matrix < 0):
        raise ValueError("Matrix contains negative values")

    if not np.allclose(matrix, matrix.T):
        raise ValueError("Matrix must be symmetric")

    n_bins = matrix.shape[0]
    start, end = region

    if not (0 <= start < end <= n_bins):
        raise ValueError("Invalid genomic region")

    if not (1 <= max_separation_bins < n_bins):
        raise ValueError("Invalid maximum separation")

    eligible_pairs = 0
    nonzero_pairs = 0
    total_signal = 0.0
    eligible_distances = 0
    supported_distances = 0

    for distance in range(1, max_separation_bins + 1):
        rows = np.arange(n_bins - distance)
        midpoints = rows + distance / 2

        selected = (
            (midpoints >= start)
            & (midpoints < end)
        )

        if not np.any(selected):
            continue

        contacts = matrix[
            rows[selected],
            rows[selected] + distance,
        ]

        eligible_distances += 1
        eligible_pairs += contacts.size
        nonzero_pairs += np.count_nonzero(contacts)
        total_signal += float(np.sum(contacts))

        if np.any(contacts > 0):
            supported_distances += 1

    return {
        "total_contact_signal": total_signal,
        "eligible_pairs": int(eligible_pairs),
        "nonzero_pairs": int(nonzero_pairs),
        "pair_occupancy": (
            float(nonzero_pairs / eligible_pairs)
            if eligible_pairs else 0.0
        ),
        "eligible_distances": eligible_distances,
        "supported_distances": supported_distances,
        "distance_support_fraction": (
            float(supported_distances / eligible_distances)
            if eligible_distances else 0.0
        ),
    }
