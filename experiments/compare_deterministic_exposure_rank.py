"""Evaluate deterministic exposure-bounded rank intervals.

Finite-sample justification assumes independent Poisson counts,
fixed positive eligible-pair exposures, and a positive total count.

This is a synthetic calibration experiment, not a biological
Hi-C confidence-interval validation.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

original = pd.read_csv(
    "research_artifacts/stage1/"
    "poisson_unequal_support_calibration.csv"
).set_index(["scenario", "depth", "replicate"])

randomized = pd.read_csv(
    "research_artifacts/stage1/"
    "equalized_poisson_rank_comparison.csv"
).set_index(["scenario", "depth", "replicate"])

assert len(original) == len(randomized) == 2400

D_MAX = 64
distances = np.arange(1, D_MAX + 1)

mask_rng = np.random.default_rng(20261028)
valid_bins = np.ones(256, dtype=bool)

valid_bins[
    mask_rng.choice(256, size=58, replace=False)
] = False

eligible = np.zeros(D_MAX, dtype=np.int64)

for d in distances:
    rows = np.arange(256 - d)
    midpoint = rows + d / 2

    selected = (
        (midpoint >= 96)
        & (midpoint < 160)
        & valid_bins[rows]
        & valid_bins[rows + d]
    )

    eligible[d - 1] = selected.sum()

assert eligible.min() == 36
assert eligible.max() == 44

k_min = int(eligible.min())
k_max = int(eligible.max())

p_low = 0.8 * k_min / (
    0.8 * k_min + 0.2 * k_max
)

p_high = 0.8 * k_max / (
    0.8 * k_max + 0.2 * k_min
)

assert p_low < 0.8 < p_high

baseline = 90 * np.exp(-distances / 18)
narrowing = 90 * np.exp(-distances / 5)

narrowing *= (
    np.sum(eligible * baseline)
    / np.sum(eligible * narrowing)
)

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

def deterministic_bounds(counts):
    n = int(counts.sum())

    if n == 0:
        return None

    lower_rank = int(
        binom.ppf(0.025, n, p_low)
    )

    upper_rank = int(
        binom.isf(0.025, n, p_high)
    ) + 1

    cumulative = np.cumsum(counts)

    lower = (
        1 if lower_rank < 1
        else int(
            np.searchsorted(cumulative, lower_rank) + 1
        )
    )

    upper = (
        D_MAX if upper_rank > n
        else int(
            np.searchsorted(cumulative, upper_rank) + 1
        )
    )

    assert 1 <= lower <= upper <= D_MAX

    return lower, upper


rng = np.random.default_rng(810241)
records = []

for scenario in ("baseline", "narrowing", "coverage_drop"):
    truth = truths[scenario]
    competitor = 28 if scenario == "narrowing" else 9

    for depth in (1.0, 0.01, 0.001, 0.0001):
        expected = eligible * intensities[scenario] * depth

        for replicate in range(200):
            key = (scenario, depth, replicate)

            counts = rng.poisson(expected).astype(np.int64)

            prior = original.loc[key]
            prior_randomized = randomized.loc[key]

            assert int(counts.sum()) == int(prior.contacts)
            assert int(counts.sum()) == int(
                prior_randomized.original_contacts
            )

            bounds = deterministic_bounds(counts)

            valid = bounds is not None

            if valid:
                lower, upper = bounds
                covered = lower <= truth <= upper
                excludes_competitor = not (
                    lower <= competitor <= upper
                )
                discriminates = covered and excludes_competitor
                span = upper - lower
            else:
                lower = np.nan
                upper = np.nan
                covered = False
                discriminates = False
                span = np.nan

            records.append({
                "scenario": scenario,
                "depth": depth,
                "replicate": replicate,
                "contacts": int(counts.sum()),
                "det_valid": valid,
                "det_covered": covered,
                "det_discriminates": discriminates,
                "det_lower": lower,
                "det_upper": upper,
                "det_span": span,
                "random_valid": bool(
                    prior_randomized.rank_valid
                ),
                "random_covered": bool(
                    prior_randomized.rank_covered
                ),
                "random_discriminates": bool(
                    prior_randomized.rank_discriminates
                ),
                "random_span": prior_randomized.rank_span,
                "poisson_span": prior.span,
            })

df = pd.DataFrame(records)
assert len(df) == 2400

output = Path(
    "research_artifacts/stage1/"
    "deterministic_exposure_rank_comparison.csv"
)

output.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output, index=False)

print("DETERMINISTIC EXPOSURE-BOUNDED RANK INTERVAL")
print("=" * 110)
print(f"Eligible-pair range: {k_min}–{k_max}")
print(f"Conservative quantile probabilities: {p_low:.5f}–{p_high:.5f}")
print("Nominal confidence: 95%")
print()

print(
    f"{'Scenario':<15} {'Depth':>8} {'Valid':>9} "
    f"{'Det cover':>10} {'Det discr':>10} "
    f"{'Rand discr':>11} {'Det span':>10} "
    f"{'Rand span':>11} {'Pois span':>11}"
)
print("-" * 110)

for (scenario, depth), group in df.groupby(
    ["scenario", "depth"], sort=False
):
    valid = group[group.det_valid]

    coverage = (
        f"{100*valid.det_covered.mean():.1f}%"
        if len(valid) else "NA"
    )

    det_span = (
        f"{valid.det_span.median():.2f}"
        if len(valid) else "NA"
    )

    print(
        f"{scenario:<15} {depth:8.4f} "
        f"{len(valid):3d}/{len(group):<5d} "
        f"{coverage:>10} "
        f"{100*group.det_discriminates.mean():9.1f}% "
        f"{100*group.random_discriminates.mean():10.1f}% "
        f"{det_span:>10} "
        f"{group.random_span.median():11.2f} "
        f"{group.poisson_span.median():11.2f}"
    )

print()
print(f"Output: {output}")
print("Discrimination uses known truth only for synthetic evaluation.")
print("Deterministic intervals use all observed contacts.")
print("No biological Hi-C coverage claim is implied.")
print("PASS: Deterministic interval comparison completed")
