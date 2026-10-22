"""Controlled synthetic Hi-C matrices for method development.

These simulations are algorithmic test fixtures, not biological data.
"""

from typing import Literal

import numpy as np


Scenario = Literal["baseline", "narrowing", "coverage_drop"]


def simulate_contact_map(
    scenario: Scenario,
    *,
    n_bins: int = 256,
    region: tuple[int, int] = (96, 160),
    seed: int = 42,
) -> np.ndarray:
    """Generate a symmetric synthetic cis-contact matrix.

    Contact pairs are assigned to a region by their genomic midpoint.

    baseline:
        Exponential contact-distance decay.

    narrowing:
        Shorter decay scale within the target region, with
        expected regional total contacts approximately preserved.

    coverage_drop:
        Uniform 75% reduction in expected regional contacts,
        without changing the underlying decay shape.
    """
    if scenario not in ("baseline", "narrowing", "coverage_drop"):
        raise ValueError(f"Unknown scenario: {scenario}")

    if n_bins < 2:
        raise ValueError("n_bins must be at least 2")

    start, end = region
    if not (0 <= start < end <= n_bins):
        raise ValueError("Invalid target region")

    rng = np.random.default_rng(seed)

    rows, cols = np.triu_indices(n_bins, k=1)

    distance = (cols - rows).astype(float)
    midpoint = (rows + cols) / 2

    in_region = (midpoint >= start) & (midpoint < end)

    # Baseline expected contact intensity.
    expected = 90.0 * np.exp(-distance / 18.0)

    if scenario == "narrowing":
        narrowed = 90.0 * np.exp(-distance / 5.0)

        # Preserve the expected total within the target region.
        scale = expected[in_region].sum() / narrowed[in_region].sum()

        expected[in_region] = narrowed[in_region] * scale

    elif scenario == "coverage_drop":
        expected[in_region] *= 0.25

    counts = rng.poisson(expected).astype(np.int32)

    matrix = np.zeros((n_bins, n_bins), dtype=np.int32)

    matrix[rows, cols] = counts
    matrix[cols, rows] = counts

    return matrix
