"""Compare Monte Carlo bootstrap intervals with exact conditional quantiles.

Applies to the synthetic calibration design with equal eligible-pair
counts at every separation. These are exact quantiles of the fitted
empirical multinomial bootstrap, not exact population confidence limits.
"""

import numpy as np
import pandas as pd
from scipy.stats import binom

from hic_diagonal_narrowing.synthetic import simulate_contact_map

data = pd.read_csv(
    "research_artifacts/stage1/bootstrap_calibration_trials.csv"
)
assert len(data) == 1200

rows, cols = np.triu_indices(256, k=1)
distance = cols - rows
midpoint = (rows + cols) / 2

local = (
    (midpoint >= 96)
    & (midpoint < 160)
    & (distance <= 64)
)

local_distance = distance[local]

eligible = np.bincount(
    local_distance, minlength=65
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

        source_counts = matrix[rows, cols].astype(np.int64)

        for di, depth in enumerate(depths):
            trial = data[
                (data.scenario == scenario)
                & (data.replicate == rep)
                & np.isclose(data.depth, depth)
            ]

            assert len(trial) == 1
            trial = trial.iloc[0]

            rng = np.random.default_rng(
                500000 + si * 100000 + di * 1000 + rep
            )

            sampled = rng.binomial(source_counts, depth)

            bands = np.bincount(
                local_distance,
                weights=sampled[local],
                minlength=65,
            )[1:].astype(np.int64)

            n = int(bands.sum())
            assert n == int(trial.contact_count)

            if n == 0:
                assert not bool(trial.valid)
                continue

            assert bool(trial.valid)

            probabilities = bands / n
            cumulative_probabilities = np.cumsum(probabilities)

            # Exact CDF of W80 under the fitted multinomial bootstrap.
            # W80 <= d iff the first d bands contain at least ceil(0.8*N).
            threshold = (4 * n + 4) // 5

            cdf = binom.sf(
                threshold - 1,
                n,
                np.clip(cumulative_probabilities, 0.0, 1.0),
            )

            assert np.all(np.diff(cdf) >= -1e-10)
            assert np.isclose(cdf[-1], 1.0)

            exact_lower = int(np.searchsorted(cdf, 0.025) + 1)
            exact_upper = int(np.searchsorted(cdf, 0.975) + 1)

            original_lower = float(trial.lower_bins)
            original_upper = float(trial.upper_bins)

            original_too_low = original_upper < truth
            exact_too_low = exact_upper < truth

            distal_probability = float(
                probabilities[truth - 1:].sum()
            )
            recovery_probability = float(
                binom.sf(n // 5, n, distal_probability)
            )

            records.append({
                "scenario": scenario,
                "depth": depth,
                "replicate": rep,
                "contacts": n,
                "original_covered": (
                    original_lower <= truth <= original_upper
                ),
                "exact_covered": (
                    exact_lower <= truth <= exact_upper
                ),
                "original_too_low": original_too_low,
                "exact_too_low": exact_too_low,
                "recovery_probability": recovery_probability,
                "exact_lower": exact_lower,
                "exact_upper": exact_upper,
            })

result = pd.DataFrame(records)
assert len(result) == 1190

print("EXACT CONDITIONAL BOOTSTRAP QUANTILE AUDIT")
print("=" * 84)
print(
    f"{'Scenario':<15} {'Depth':>8} {'Valid':>6} "
    f"{'MC cover':>10} {'Exact cover':>12} "
    f"{'MC low':>8} {'Exact low':>10}"
)

for (scenario, depth), g in result.groupby(
    ["scenario", "depth"], sort=False
):
    print(
        f"{scenario:<15} {depth:8.4f} {len(g):6d} "
        f"{100*g.original_covered.mean():9.1f}% "
        f"{100*g.exact_covered.mean():11.1f}% "
        f"{int(g.original_too_low.sum()):8d} "
        f"{int(g.exact_too_low.sum()):10d}"
    )

print("\nOVERALL")
print(
    "Original downward failures:",
    int(result.original_too_low.sum())
)
print(
    "Analytic-quantile downward failures:",
    int(result.exact_too_low.sum())
)

residual = result[
    result.original_too_low
    & (result.recovery_probability >= 0.025)
]

assert len(residual) == 6

print("\nSIX PREVIOUSLY UNEXPLAINED CASES")
print(
    residual[
        [
            "scenario", "depth", "replicate",
            "recovery_probability", "exact_upper", "exact_covered"
        ]
    ].to_string(index=False))

print()
print("MC = original 300-resample percentile interval.")
print("Exact = analytic quantiles of the fitted empirical bootstrap.")
print("Neither method is biologically validated.")
print("PASS: Conditional quantile comparison completed")
