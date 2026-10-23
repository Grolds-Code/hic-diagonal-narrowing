"""Bootstrap uncertainty for local Hi-C diagonal width.

Research prototype v0.1.

Estimates conditional contact-sampling variability.
Does not account for missing bins, mapping bias, or
systematic biological and technical uncertainty.
"""

import numpy as np

from hic_diagonal_narrowing.width import estimate_w80


def bootstrap_w80(
    matrix: np.ndarray,
    region: tuple[int, int],
    *,
    max_separation_bins: int = 64,
    n_bootstraps: int = 300,
    confidence_level: float = 0.95,
    seed: int = 42,
) -> dict:
    """Estimate W80 variability by multinomial resampling."""

    matrix = np.asarray(matrix)

    if not np.issubdtype(matrix.dtype, np.integer):
        raise TypeError("Raw integer contact counts required")

    if n_bootstraps < 20:
        raise ValueError("At least 20 bootstrap replicates required")

    if not 0 < confidence_level < 1:
        raise ValueError("Invalid confidence level")

    observed = estimate_w80(
        matrix,
        region,
        max_separation_bins=max_separation_bins,
    )

    n_bins = matrix.shape[0]
    start, end = region

    band_counts = np.zeros(
        max_separation_bins, dtype=np.int64
    )
    eligible_pairs = np.zeros(
        max_separation_bins, dtype=np.int64
    )

    for distance in range(1, max_separation_bins + 1):
        rows = np.arange(n_bins - distance)
        midpoints = rows + distance / 2

        selected = (
            (midpoints >= start)
            & (midpoints < end)
        )

        contacts = matrix[
            rows[selected],
            rows[selected] + distance,
        ]

        band_counts[distance - 1] = contacts.sum()
        eligible_pairs[distance - 1] = contacts.size

    if np.any(eligible_pairs == 0):
        raise ValueError("Insufficient eligible pairs")

    total_contacts = int(band_counts.sum())

    if total_contacts <= 0:
        raise ValueError("No informative contact signal")

    probabilities = band_counts / total_contacts

    rng = np.random.default_rng(seed)

    samples = rng.multinomial(
        total_contacts,
        probabilities,
        size=n_bootstraps,
    )

    profiles = samples / eligible_pairs
    cumulative = np.cumsum(profiles, axis=1)
    targets = 0.80 * profiles.sum(axis=1)

    widths = (
        np.sum(
            cumulative < targets[:, None],
            axis=1,
        ) + 1
    )

    alpha = 1 - confidence_level

    lower, upper = np.quantile(
        widths,
        [alpha / 2, 1 - alpha / 2],
    )

    return {
        "observed_w80": int(observed),
        "bootstrap_lower": float(lower),
        "bootstrap_upper": float(upper),
        "bootstrap_std": float(
            np.std(widths, ddof=1)
        ),
        "total_contacts": total_contacts,
        "n_bootstraps": n_bootstraps,
    }
