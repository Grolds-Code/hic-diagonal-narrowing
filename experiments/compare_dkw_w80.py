"""Compare empirical-bootstrap and DKW intervals for synthetic W80.

The DKW construction applies here because the synthetic W80 is
the 80th percentile of a contact-separation distribution with
equal eligible-pair counts per distance.

This is not a validated confidence interval for biological Hi-C.
"""

from math import log, sqrt

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map

df = pd.read_csv(
    "research_artifacts/stage1/bootstrap_calibration_trials.csv"
)

assert len(df) == 1200

rows, cols = np.triu_indices(256, k=1)
distances = cols - rows
midpoints = (rows + cols) / 2

local = (
    (midpoints >= 96)
    & (midpoints < 160)
    & (distances <= 64)
)

eligible = np.bincount(
    distances[local], minlength=65
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
            trial = df[
                (df.scenario == scenario)
                & (df.replicate == rep)
                & np.isclose(df.depth, depth)
            ]

            assert len(trial) == 1
            trial = trial.iloc[0]

            rng = np.random.default_rng(
                500000 + si * 100000 + di * 1000 + rep
            )

            sampled = rng.binomial(source, depth)

            bands = np.bincount(
                distances[local],
                weights=sampled[local],
                minlength=65,
            )[1:].astype(np.int64)

            n = int(bands.sum())

            assert n == int(trial.contact_count)

            if n == 0:
                assert not bool(trial.valid)
                continue

            assert bool(trial.valid)

            empirical_cdf = np.cumsum(bands) / n

            alpha = 0.05
            target = 0.80

            epsilon = sqrt(
                log(2.0 / alpha) / (2.0 * n)
            )

            if target - epsilon <= 0:
                lower = 1
            else:
                lower = int(
                    np.searchsorted(
                        empirical_cdf,
                        target - epsilon,
                        side="left",
                    ) + 1
                )

            if target + epsilon > 1:
                upper = 64
            else:
                upper = int(
                    np.searchsorted(
                        empirical_cdf,
                        target + epsilon,
                        side="left",
                    ) + 1
                )

            assert 1 <= lower <= upper <= 64

            records.append({
                "scenario": scenario,
                "depth": depth,
                "replicate": rep,
                "contacts": n,
                "bootstrap_covered": (
                    float(trial.lower_bins)
                    <= truth
                    <= float(trial.upper_bins)
                ),
                "bootstrap_span": (
                    float(trial.upper_bins)
                    - float(trial.lower_bins)
                ),
                "dkw_covered": lower <= truth <= upper,
                "dkw_span": upper - lower,
                "dkw_lower": lower,
                "dkw_upper": upper,
            })

result = pd.DataFrame(records)

assert len(result) == 1190

print("W80 UNCERTAINTY METHOD COMPARISON")
print("=" * 99)
print(
    f"{'Scenario':<15} {'Depth':>8} {'Valid':>7} "
    f"{'Boot cover':>11} {'DKW cover':>11} "
    f"{'Boot span':>11} {'DKW span':>11}"
)

for (scenario, depth), g in result.groupby(
    ["scenario", "depth"], sort=False
):
    print(
        f"{scenario:<15} {depth:8.4f} "
        f"{len(g):7d} "
        f"{100*g.bootstrap_covered.mean():10.1f}% "
        f"{100*g.dkw_covered.mean():10.1f}% "
        f"{g.bootstrap_span.median():11.2f} "
        f"{g.dkw_span.median():11.2f}"
    )

print("\nOVERALL")
print(f"Valid trials: {len(result)}")
print(
    "Bootstrap coverage: "
    f"{100*result.bootstrap_covered.mean():.1f}%"
)
print(
    "DKW coverage: "
    f"{100*result.dkw_covered.mean():.1f}%"
)
print(
    "Median bootstrap span: "
    f"{result.bootstrap_span.median():.2f} bins"
)
print(
    "Median DKW span: "
    f"{result.dkw_span.median():.2f} bins"
)

print("\nCoverage is conditional on a positive contact count.")
print("Interval spans are reported in 500-bp genomic bins.")
print("DKW validity here relies on the synthetic sampling model.")
print("PASS: DKW comparison completed")
