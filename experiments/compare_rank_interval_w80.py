"""Compare three interval methods for synthetic W80.

The rank interval uses exact binomial tail probabilities.
Its >=95% finite-sample coverage guarantee applies to the
synthetic independent-contact model with equal eligible-pair
counts per separation, conditional on a positive contact total.

No biological Hi-C coverage guarantee is implied.
"""

from math import log, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

from hic_diagonal_narrowing.synthetic import simulate_contact_map

ALPHA = 0.05
P = 0.80
MAX_DISTANCE = 64

trials = pd.read_csv(
    "research_artifacts/stage1/bootstrap_calibration_trials.csv"
)

assert len(trials) == 1200

rows, cols = np.triu_indices(256, k=1)
distance = cols - rows
midpoint = (rows + cols) / 2

local = (
    (midpoint >= 96)
    & (midpoint < 160)
    & (distance <= MAX_DISTANCE)
)

eligible = np.bincount(
    distance[local], minlength=65
)[1:]

assert np.all(eligible == 64)

scenarios = [
    ("baseline", 28),
    ("narrowing", 9),
    ("coverage_drop", 28),
]

depths = (1.0, 0.01, 0.001, 0.0001)

records = []

for si, (scenario, truth) in enumerate(scenarios):
    for rep in range(100):
        matrix = simulate_contact_map(
            scenario,
            seed=42000 + si * 10000 + rep,
        )

        source = matrix[rows, cols].astype(np.int64)

        for di, depth in enumerate(depths):
            selected = trials[
                (trials.scenario == scenario)
                & (trials.replicate == rep)
                & np.isclose(trials.depth, depth)
            ]

            assert len(selected) == 1
            trial = selected.iloc[0]

            rng = np.random.default_rng(
                500000 + si * 100000 + di * 1000 + rep
            )

            sampled = rng.binomial(source, depth)

            bands = np.bincount(
                distance[local],
                weights=sampled[local],
                minlength=65,
            )[1:].astype(np.int64)

            n = int(bands.sum())

            assert n == int(trial.contact_count)

            if n == 0:
                assert not bool(trial.valid)
                continue

            assert bool(trial.valid)

            counts_cumulative = np.cumsum(bands)
            empirical_cdf = counts_cumulative / n

            # DKW inversion: simultaneous CDF confidence band.
            epsilon = sqrt(
                log(2 / ALPHA) / (2 * n)
            )

            if P - epsilon <= 0:
                dkw_lower = 1
            else:
                dkw_lower = int(
                    np.searchsorted(
                        empirical_cdf, P - epsilon
                    ) + 1
                )

            if P + epsilon > 1:
                dkw_upper = MAX_DISTANCE
            else:
                dkw_upper = int(
                    np.searchsorted(
                        empirical_cdf, P + epsilon
                    ) + 1
                )

            # Equal-tailed binomial rank construction.
            # Each one-sided noncoverage probability is
            # bounded by alpha/2, including with ties.
            lower_rank = int(
                binom.ppf(ALPHA / 2, n, P)
            )

            upper_rank = int(
                binom.isf(ALPHA / 2, n, P)
            ) + 1

            if lower_rank < 1:
                rank_lower = 1
            else:
                rank_lower = int(
                    np.searchsorted(
                        counts_cumulative, lower_rank
                    ) + 1
                )

            if upper_rank > n:
                rank_upper = MAX_DISTANCE
            else:
                rank_upper = int(
                    np.searchsorted(
                        counts_cumulative, upper_rank
                    ) + 1
                )

            assert (
                1 <= rank_lower <= rank_upper <= MAX_DISTANCE
            )

            assert (
                1 <= dkw_lower <= dkw_upper <= MAX_DISTANCE
            )

            boot_lower = float(trial.lower_bins)
            boot_upper = float(trial.upper_bins)

            records.append({
                "scenario": scenario,
                "depth": depth,
                "replicate": rep,
                "contacts": n,
                "population_w80": truth,
                "bootstrap_covered": (
                    boot_lower <= truth <= boot_upper
                ),
                "dkw_covered": (
                    dkw_lower <= truth <= dkw_upper
                ),
                "rank_covered": (
                    rank_lower <= truth <= rank_upper
                ),
                "bootstrap_span": boot_upper - boot_lower,
                "dkw_span": dkw_upper - dkw_lower,
                "rank_span": rank_upper - rank_lower,
                "rank_lower": rank_lower,
                "rank_upper": rank_upper,
            })

result = pd.DataFrame(records)
assert len(result) == 1190

output = Path(
    "research_artifacts/stage1/"
    "w80_interval_method_comparison.csv"
)

output.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(output, index=False)

print("W80 INTERVAL METHOD COMPARISON")
print("=" * 110)
print(
    f"{'Scenario':<15} {'Depth':>8} {'Valid':>6} "
    f"{'Boot %':>8} {'DKW %':>8} {'Rank %':>8} "
    f"{'Boot span':>11} {'DKW span':>11} {'Rank span':>11}"
)

for (scenario, depth), group in result.groupby(
    ["scenario", "depth"], sort=False
):
    print(
        f"{scenario:<15} "
        f"{depth:8.4f} "
        f"{len(group):6d} "
        f"{100*group.bootstrap_covered.mean():7.1f}% "
        f"{100*group.dkw_covered.mean():7.1f}% "
        f"{100*group.rank_covered.mean():7.1f}% "
        f"{group.bootstrap_span.median():11.2f} "
        f"{group.dkw_span.median():11.2f} "
        f"{group.rank_span.median():11.2f}"
    )

print("\nOVERALL")
for method in ("bootstrap", "dkw", "rank"):
    print(
        f"{method:<10} "
        f"coverage={100*result[f'{method}_covered'].mean():5.1f}% "
        f"median_span={result[f'{method}_span'].median():6.2f} bins"
    )

print(f"\nOutput: {output}")
print("Failures from zero observed contacts are excluded.")
print("Finite-sample rank validity depends on the simulation model.")
print("PASS: Three-method interval comparison completed")
