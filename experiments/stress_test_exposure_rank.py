"""Stress-test deterministic W80 intervals across exposure profiles.

Independent Poisson band counts with fixed, unequal exposures.
Tests coverage, interval width and synthetic discrimination.

No biological Hi-C validation is implied.
"""

import numpy as np
from scipy.stats import binom

D_MAX = 64
DISTANCES = np.arange(1, D_MAX + 1)
N_REPLICATES = 300

SCENARIOS = ("baseline", "narrowing", "coverage_drop")
DEPTHS = (1.0, 0.01, 0.001)

mask_rng = np.random.default_rng(20261028)

valid_bins = np.ones(256, dtype=bool)
valid_bins[
    mask_rng.choice(256, size=58, replace=False)
] = False

masked_exposure = []

for d in DISTANCES:
    rows = np.arange(256 - d)

    selected = (
        (rows + d / 2 >= 96)
        & (rows + d / 2 < 160)
        & valid_bins[rows]
        & valid_bins[rows + d]
    )

    masked_exposure.append(int(selected.sum()))

exposures = {
    "equal": np.full(D_MAX, 40, dtype=np.int64),
    "masked": np.array(masked_exposure, dtype=np.int64),
    "decreasing": np.rint(
        np.linspace(70, 10, D_MAX)
    ).astype(np.int64),
    "increasing": np.rint(
        np.linspace(10, 70, D_MAX)
    ).astype(np.int64),
}

assert exposures["masked"].min() == 36
assert exposures["masked"].max() == 44

baseline = 90 * np.exp(-DISTANCES / 18)
narrowing = 90 * np.exp(-DISTANCES / 5)

intensities = {
    "baseline": baseline,
    "narrowing": narrowing,
    "coverage_drop": 0.25 * baseline,
}

truths = {
    "baseline": 28,
    "narrowing": 9,
    "coverage_drop": 28,
}


def interval(counts, eligible):
    n = int(counts.sum())

    if n == 0:
        return None

    k_min = eligible.min()
    k_max = eligible.max()

    p_low = (
        0.8 * k_min
        / (0.8 * k_min + 0.2 * k_max)
    )

    p_high = (
        0.8 * k_max
        / (0.8 * k_max + 0.2 * k_min)
    )

    lower_rank = int(binom.ppf(0.025, n, p_low))
    upper_rank = int(binom.isf(0.025, n, p_high)) + 1

    cumulative = np.cumsum(counts)

    lower = (
        1 if lower_rank < 1
        else int(np.searchsorted(cumulative, lower_rank) + 1)
    )

    upper = (
        D_MAX if upper_rank > n
        else int(np.searchsorted(cumulative, upper_rank) + 1)
    )

    assert 1 <= lower <= upper <= D_MAX

    return lower, upper


rng = np.random.default_rng(20261101)

print("DETERMINISTIC W80 INTERVAL — EXPOSURE STRESS TEST")
print("=" * 105)
print(f"Replicates per condition: {N_REPLICATES}")
print("Nominal confidence level: 95%")
print()

print(
    f"{'Exposure':<13} {'Scenario':<15} {'Depth':>8} "
    f"{'K ratio':>8} {'Valid':>9} "
    f"{'Coverage':>10} {'Discrim.':>10} "
    f"{'Med span':>10}"
)
print("-" * 105)

for exposure_name, eligible in exposures.items():
    assert np.all(eligible > 0)

    ratio = eligible.max() / eligible.min()

    for scenario in SCENARIOS:
        truth = truths[scenario]
        competitor = 28 if scenario == "narrowing" else 9

        for depth in DEPTHS:
            expected = (
                eligible * intensities[scenario] * depth
            )

            valid = 0
            covered = 0
            discriminating = 0
            spans = []

            for _ in range(N_REPLICATES):
                counts = rng.poisson(expected)

                bounds = interval(counts, eligible)

                if bounds is None:
                    continue

                valid += 1
                lower, upper = bounds

                contains_truth = (
                    lower <= truth <= upper
                )

                excludes_competitor = not (
                    lower <= competitor <= upper
                )

                covered += int(contains_truth)

                discriminating += int(
                    contains_truth and excludes_competitor
                )

                spans.append(upper - lower)

            coverage = (
                f"{100 * covered / valid:.1f}%"
                if valid else "NA"
            )

            discrimination = (
                100 * discriminating / N_REPLICATES
            )

            median_span = (
                f"{np.median(spans):.2f}"
                if spans else "NA"
            )

            print(
                f"{exposure_name:<13} "
                f"{scenario:<15} "
                f"{depth:8.4f} "
                f"{ratio:8.2f} "
                f"{valid:3d}/{N_REPLICATES:<5d} "
                f"{coverage:>10} "
                f"{discrimination:9.1f}% "
                f"{median_span:>10}"
            )

    print()

print("Coverage is conditional on a valid interval.")
print("Discrimination is measured across all attempts.")
print("K ratio = maximum exposure / minimum exposure.")
print("Synthetic model assumes independent Poisson band counts.")
print("PASS: Exposure-imbalance stress test completed")
