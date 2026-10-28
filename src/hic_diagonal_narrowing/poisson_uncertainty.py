"""Simultaneous Poisson confidence bounds for W80.

Research methodology for raw contact counts.

Uses Bonferroni-adjusted, exact Poisson mean intervals
followed by inversion of the cumulative distance profile.

Supports unequal eligible-pair counts across separations.

Coverage guarantees require the specified Poisson sampling
model with fixed eligible-pair exposures. They do not
automatically apply to biological Hi-C measurements.
"""

import numpy as np
from scipy.stats import chi2


def poisson_w80_interval_from_bands(
    band_counts: np.ndarray,
    eligible_pairs: np.ndarray,
    *,
    confidence_level: float = 0.95,
) -> dict:
    """Bound W80 using simultaneous Poisson intensity intervals."""

    counts = np.asarray(band_counts)
    eligible = np.asarray(eligible_pairs)

    if (
        counts.ndim != 1
        or eligible.ndim != 1
        or counts.shape != eligible.shape
        or counts.size == 0
    ):
        raise ValueError(
            "Band counts and eligible pairs must be matching "
            "nonempty one-dimensional arrays"
        )

    if not np.issubdtype(counts.dtype, np.integer):
        raise TypeError("Raw integer contact counts required")

    if not np.issubdtype(eligible.dtype, np.integer):
        raise TypeError("Integer eligible-pair counts required")

    if np.any(counts < 0):
        raise ValueError("Contact counts must be nonnegative")

    if np.any(eligible <= 0):
        raise ValueError(
            "Every separation must contain eligible pairs"
        )

    if not np.isfinite(confidence_level) or not (
        0 < confidence_level < 1
    ):
        raise ValueError("Invalid confidence level")

    if int(counts.sum()) == 0:
        raise ValueError("No informative contact signal")

    counts = counts.astype(np.float64)
    eligible = eligible.astype(np.float64)

    n_distances = counts.size
    alpha = 1.0 - confidence_level

    # Bonferroni correction over both tails of all bands.
    tail_alpha = alpha / (2 * n_distances)

    # Garwood exact Poisson confidence limits.
    lambda_lower = np.zeros(n_distances)

    positive = counts > 0

    lambda_lower[positive] = (
        0.5 * chi2.ppf(
            tail_alpha,
            2 * counts[positive],
        )
    )

    lambda_upper = (
        0.5 * chi2.isf(
            tail_alpha,
            2 * (counts + 1),
        )
    )

    mean_lower = lambda_lower / eligible
    mean_upper = lambda_upper / eligible

    # Lower bound on cumulative intensity fraction:
    # use lower intensity in the prefix and upper in the tail.
    cumulative_lower = np.cumsum(mean_lower)
    cumulative_upper = np.cumsum(mean_upper)

    remaining_lower = np.maximum(
        0.0,
        mean_lower.sum() - cumulative_lower,
    )
    remaining_upper = np.maximum(
        0.0,
        mean_upper.sum() - cumulative_upper,
    )

    fraction_lower = np.divide(
        cumulative_lower,
        cumulative_lower + remaining_upper,
        out=np.zeros(n_distances),
        where=(
            cumulative_lower + remaining_upper
        ) > 0,
    )

    fraction_upper = np.divide(
        cumulative_upper,
        cumulative_upper + remaining_lower,
        out=np.ones(n_distances),
        where=(
            cumulative_upper + remaining_lower
        ) > 0,
    )

    # At the final band, the cumulative fraction is exactly 1
    # whenever the population has positive total intensity.
    fraction_lower[-1] = 1.0
    fraction_upper[-1] = 1.0

    profile = counts / eligible
    observed_w80 = int(
        np.searchsorted(
            np.cumsum(profile),
            0.80 * profile.sum(),
        ) + 1
    )

    # Earliest separation at which crossing 80% is possible.
    lower_w80 = int(
        np.searchsorted(
            fraction_upper,
            0.80,
        ) + 1
    )

    # Earliest separation at which crossing 80% is guaranteed.
    upper_w80 = int(
        np.searchsorted(
            fraction_lower,
            0.80,
        ) + 1
    )

    assert (
        1 <= lower_w80
        <= observed_w80
        <= upper_w80
        <= n_distances
    )

    return {
        "observed_w80": observed_w80,
        "poisson_lower": lower_w80,
        "poisson_upper": upper_w80,
        "total_contacts": int(counts.sum()),
        "n_distances": n_distances,
        "confidence_level": confidence_level,
        "method": "bonferroni_poisson",
    }


def cutoff_exposure_w80_interval_from_bands(
    band_counts: np.ndarray,
    eligible_pairs: np.ndarray,
    *,
    confidence_level: float = 0.95,
) -> dict:
    """Deterministic W80 confidence interval with unequal exposures.

    Population target:
        W80 of the distance-specific mean contact intensities.

    Sampling assumptions:
        Independent Poisson band totals with means K_d * mu_d,
        where positive eligible-pair exposures K_d are fixed.

    Conditional on a positive total count, the interval has
    at least the requested coverage under this model.

    This guarantee does not automatically extend to biological
    Hi-C counts, data-dependent masks, or dependent contacts.
    """
    from scipy.stats import binom

    counts = np.asarray(band_counts)
    exposures = np.asarray(eligible_pairs)

    if (
        counts.ndim != 1
        or exposures.ndim != 1
        or counts.shape != exposures.shape
        or counts.size == 0
    ):
        raise ValueError(
            "Counts and exposures must be matching "
            "nonempty one-dimensional arrays"
        )

    if not np.issubdtype(counts.dtype, np.integer):
        raise TypeError("Raw integer contact counts required")

    if not np.issubdtype(exposures.dtype, np.integer):
        raise TypeError("Integer eligible-pair counts required")

    if np.any(counts < 0):
        raise ValueError("Negative contact counts")

    if np.any(exposures <= 0):
        raise ValueError(
            "Every separation must have positive exposure"
        )

    if not np.isfinite(confidence_level) or not (
        0 < confidence_level < 1
    ):
        raise ValueError("Invalid confidence level")

    total = sum(int(value) for value in counts)

    if total == 0:
        raise ValueError("No informative contact signal")

    # Avoid silently losing integer precision when calculating
    # statistical quantiles or cumulative contact counts.
    if total > 2**53 - 1:
        raise ValueError(
            "Total contacts exceed exact floating-point "
            "integer precision"
        )

    if any(int(value) > 2**53 - 1 for value in exposures):
        raise ValueError(
            "Exposures exceed exact floating-point "
            "integer precision"
        )

    counts = counts.astype(np.int64)
    exposures = exposures.astype(np.float64)

    d_max = counts.size
    alpha = 1.0 - confidence_level
    target = 0.80

    profile = counts / exposures
    observed = int(
        np.searchsorted(
            np.cumsum(profile),
            target * profile.sum(),
        ) + 1
    )

    # Candidate cutoffs are 1,...,D_max-1.
    # The last separation always has cumulative fraction 1.
    prefix_min = np.minimum.accumulate(exposures)[:-1]
    prefix_max = np.maximum.accumulate(exposures)[:-1]

    suffix_min = np.minimum.accumulate(
        exposures[::-1]
    )[::-1][1:]

    suffix_max = np.maximum.accumulate(
        exposures[::-1]
    )[::-1][1:]

    p_low = (
        target * prefix_min
        / (
            target * prefix_min
            + (1 - target) * suffix_max
        )
    )

    p_high = (
        target * prefix_max
        / (
            target * prefix_max
            + (1 - target) * suffix_min
        )
    )

    if not (
        np.all(np.isfinite(p_low))
        and np.all(np.isfinite(p_high))
        and np.all(p_low > 0)
        and np.all(p_high < 1)
        and np.all(p_low <= p_high)
    ):
        raise ValueError("Invalid exposure probability bounds")

    lower_thresholds = binom.ppf(
        alpha / 2, total, p_low
    )

    upper_thresholds = binom.isf(
        alpha / 2, total, p_high
    )

    cumulative = np.cumsum(counts)[:-1]

    lower_candidates = np.flatnonzero(
        cumulative >= lower_thresholds
    )

    upper_candidates = np.flatnonzero(
        cumulative <= upper_thresholds
    )

    lower = (
        int(lower_candidates[0]) + 1
        if lower_candidates.size
        else d_max
    )

    upper = (
        int(upper_candidates[-1]) + 2
        if upper_candidates.size
        else 1
    )

    if not 1 <= lower <= upper <= d_max:
        raise AssertionError(
            "Confidence-bound ordering invariant violated"
        )

    return {
        "observed_w80": observed,
        "cutoff_lower": lower,
        "cutoff_upper": upper,
        "total_contacts": total,
        "n_distances": d_max,
        "exposure_min": int(exposures.min()),
        "exposure_max": int(exposures.max()),
        "confidence_level": confidence_level,
        "method": "cutoff_exposure_rank",
    }


def cutoff_exposure_w80_interval(
    matrix: np.ndarray,
    region: tuple[int, int],
    *,
    max_separation_bins: int = 64,
    valid_bins: np.ndarray | None = None,
    confidence_level: float = 0.95,
) -> dict:
    """Compute cutoff-specific W80 bounds from a genomic contact matrix.

    Uses midpoint-based window assignment and excludes every pair
    touching an invalid genomic bin.

    Requires raw integer contact counts. Missing genomic bins may
    be masked, but the mask must be specified independently of
    the observed contact-count values for the stated fixed-exposure
    statistical guarantee to apply.

    The Poisson coverage guarantee is conditional on the sampling
    assumptions, not established for biological Hi-C.
    """
    from hic_diagonal_narrowing.width import estimate_w80

    matrix = np.asarray(matrix)

    if not np.issubdtype(matrix.dtype, np.integer):
        raise TypeError("Raw integer contact counts required")

    # Reuse the authoritative estimator's matrix, region,
    # symmetry, mask, and complete-distance-support validation.
    observed = estimate_w80(
        matrix,
        region,
        max_separation_bins=max_separation_bins,
        valid_bins=valid_bins,
    )

    n_bins = matrix.shape[0]

    if valid_bins is None:
        mask = np.ones(n_bins, dtype=bool)
    else:
        mask = np.asarray(valid_bins)

    start, end = region

    counts = np.zeros(max_separation_bins, dtype=np.int64)
    exposures = np.zeros(max_separation_bins, dtype=np.int64)

    for distance in range(1, max_separation_bins + 1):
        rows = np.arange(n_bins - distance)
        midpoints = rows + distance / 2

        selected = (
            (midpoints >= start)
            & (midpoints < end)
            & mask[rows]
            & mask[rows + distance]
        )

        selected_rows = rows[selected]

        contacts = matrix[
            selected_rows,
            selected_rows + distance,
        ]

        exposures[distance - 1] = len(contacts)

        # Python integer accumulation avoids intermediate
        # NumPy integer overflow in the contact sum.
        band_total = sum(int(value) for value in contacts)

        if band_total > 2**53 - 1:
            raise ValueError(
                "Band total exceeds exact floating-point "
                "integer precision"
            )

        counts[distance - 1] = band_total

    result = cutoff_exposure_w80_interval_from_bands(
        counts,
        exposures,
        confidence_level=confidence_level,
    )

    if result["observed_w80"] != observed:
        raise AssertionError(
            "Matrix and band-level W80 estimators disagree"
        )

    return {
        **result,
        "eligible_pairs_by_distance": exposures.tolist(),
        "contact_counts_by_distance": counts.tolist(),
    }
