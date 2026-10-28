"""Compare Poisson W80 intervals with equal-exposure rank intervals.

Both methods use the same synthetic Poisson realizations.

Equal-exposure thinning produces independent Poisson counts with
identical effective exposure at each separation.

Rank-interval coverage is justified conditional on a positive
retained contact total under the independent Poisson model.

No biological Hi-C calibration is implied.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

reference = pd.read_csv(
    "research_artifacts/stage1/"
    "poisson_unequal_support_calibration.csv"
)

assert len(reference) == 2400

keys = ["scenario", "depth", "replicate"]
assert not reference.duplicated(keys).any()
reference = reference.set_index(keys)

N_BINS = 256
REGION = (96, 160)
D_MAX = 64
N_REPLICATES = 200

SCENARIOS = ("baseline", "narrowing", "coverage_drop")
DEPTHS = (1.0, 0.01, 0.001, 0.0001)

distances = np.arange(1, D_MAX + 1)

mask_rng = np.random.default_rng(20261028)
valid_bins = np.ones(N_BINS, dtype=bool)

invalid = mask_rng.choice(
    N_BINS, size=58, replace=False
)
valid_bins[invalid] = False

eligible = np.zeros(D_MAX, dtype=np.int64)

for d in distances:
    rows = np.arange(N_BINS - d)
    midpoint = rows + d / 2

    selected = (
        (midpoint >= REGION[0])
        & (midpoint < REGION[1])
        & valid_bins[rows]
        & valid_bins[rows + d]
    )

    eligible[d - 1] = selected.sum()

assert eligible.min() == 36
assert eligible.max() == 44

effective_exposure = int(eligible.min())
retention_probabilities = effective_exposure / eligible

baseline = 90.0 * np.exp(-distances / 18.0)
narrowing = 90.0 * np.exp(-distances / 5.0)

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

def rank_bounds(counts):
    total = int(counts.sum())

    if total <= 0:
        return None

    alpha = 0.05
    target = 0.80

    lower_rank = int(
        binom.ppf(alpha / 2, total, target)
    )

    upper_rank = int(
        binom.isf(alpha / 2, total, target)
    ) + 1

    cumulative = np.cumsum(counts)

    lower = (
        1 if lower_rank < 1
        else int(np.searchsorted(cumulative, lower_rank) + 1)
    )

    upper = (
        D_MAX if upper_rank > total
        else int(np.searchsorted(cumulative, upper_rank) + 1)
    )

    assert 1 <= lower <= upper <= D_MAX

    return lower, upper


count_rng = np.random.default_rng(810241)
thin_rng = np.random.default_rng(20261029)

records = []

for scenario in SCENARIOS:
    truth = truths[scenario]
    alternative = 28 if scenario == "narrowing" else 9

    intensity = intensities[scenario]

    for depth in DEPTHS:
        expected = eligible * intensity * depth

        for rep in range(N_REPLICATES):
            counts = count_rng.poisson(expected).astype(
                np.int64
            )

            key = (scenario, depth, rep)
            original = reference.loc[key]

            total = int(counts.sum())

            assert total == int(original["contacts"])
            assert bool(original["valid"]) == (total > 0)

            thinned = thin_rng.binomial(
                counts,
                retention_probabilities,
            )

            retained = int(thinned.sum())
            bounds = rank_bounds(thinned)

            rank_valid = bounds is not None

            if rank_valid:
                lower, upper = bounds
                rank_covered = lower <= truth <= upper
                rank_discriminates = (
                    rank_covered
                    and not (lower <= alternative <= upper)
                )
                rank_span = upper - lower
            else:
                lower = np.nan
                upper = np.nan
                rank_covered = False
                rank_discriminates = False
                rank_span = np.nan

            poisson_valid = bool(original["valid"])

            if poisson_valid:
                p_lower = int(original["lower"])
                p_upper = int(original["upper"])

                poisson_discriminates = (
                    p_lower <= truth <= p_upper
                    and not (
                        p_lower <= alternative <= p_upper
                    )
                )
            else:
                poisson_discriminates = False

            records.append({
                "scenario": scenario,
                "depth": depth,
                "replicate": rep,
                "original_contacts": total,
                "retained_contacts": retained,
                "poisson_valid": poisson_valid,
                "rank_valid": rank_valid,
                "poisson_covered": bool(
                    original["covered"]
                ),
                "rank_covered": rank_covered,
                "poisson_discriminates":
                    poisson_discriminates,
                "rank_discriminates":
                    rank_discriminates,
                "poisson_span": original["span"],
                "rank_span": rank_span,
                "rank_lower": lower,
                "rank_upper": upper,
            })

df = pd.DataFrame(records)

assert len(df) == 2400

output = Path(
    "research_artifacts/stage1/"
    "equalized_poisson_rank_comparison.csv"
)

output.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output, index=False)

print("UNEQUAL-SUPPORT W80 UNCERTAINTY COMPARISON")
print("=" * 112)
print(f"Effective exposure per distance: {effective_exposure}")
print(
    "Thinning probabilities: "
    f"{retention_probabilities.min():.3f}–"
    f"{retention_probabilities.max():.3f}"
)
print("Nominal confidence: 95%")
print()

print(
    f"{'Scenario':<15} {'Depth':>8} "
    f"{'Rank valid':>10} "
    f"{'Rank cover':>11} "
    f"{'Pois discr':>11} "
    f"{'Rank discr':>11} "
    f"{'Pois span':>11} "
    f"{'Rank span':>11}"
)
print("-" * 112)

for (scenario, depth), group in df.groupby(
    ["scenario", "depth"], sort=False
):
    valid = group[group.rank_valid]

    coverage = (
        f"{100*valid.rank_covered.mean():.1f}%"
        if len(valid) else "NA"
    )

    rank_span = (
        f"{valid.rank_span.median():.2f}"
        if len(valid) else "NA"
    )

    print(
        f"{scenario:<15} "
        f"{depth:8.4f} "
        f"{len(valid):3d}/{len(group):<5d} "
        f"{coverage:>11} "
        f"{100*group.poisson_discriminates.mean():10.1f}% "
        f"{100*group.rank_discriminates.mean():10.1f}% "
        f"{group.poisson_span.median():11.2f} "
        f"{rank_span:>11}"
    )

print()
print(f"Output: {output}")
print("Coverage is conditional on a valid rank interval.")
print("Discrimination is measured across all trials.")
print("Spans are in genomic bins.")
print("No biological inference or two-window significance test.")
print("PASS: Equal-exposure rank comparison completed")
