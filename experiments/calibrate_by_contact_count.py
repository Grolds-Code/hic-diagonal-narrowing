"""Explore W80 bootstrap reliability by informative contact count.

Synthetic research experiment, not biological validation.
"""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.uncertainty import bootstrap_w80


REGION = (96, 160)
MAX_DISTANCE = 64
N_REPLICATES = 100
N_BOOTSTRAPS = 200

SCENARIOS = {
    "baseline": 28,
    "narrowing": 9,
    "coverage_drop": 28,
}

DEPTHS = (1.0, 0.01, 0.001, 0.0001)

ROWS, COLS = np.triu_indices(256, k=1)

MIDPOINTS = (ROWS + COLS) / 2
DISTANCES = COLS - ROWS

LOCAL = (
    (MIDPOINTS >= REGION[0])
    & (MIDPOINTS < REGION[1])
    & (DISTANCES >= 1)
    & (DISTANCES <= MAX_DISTANCE)
)


def thin_matrix(matrix, fraction, rng):
    counts = matrix[ROWS, COLS].astype(np.int64)

    sampled = rng.binomial(counts, fraction)

    result = np.zeros_like(matrix)
    result[ROWS, COLS] = sampled
    result[COLS, ROWS] = sampled

    return result


def main():
    records = []

    for scenario_index, (scenario, truth) in enumerate(
        SCENARIOS.items()
    ):
        for replicate in range(N_REPLICATES):

            matrix = simulate_contact_map(
                scenario,
                seed=70000 + scenario_index * 10000 + replicate,
            )

            for depth_index, depth in enumerate(DEPTHS):

                seed = (
                    800000
                    + scenario_index * 100000
                    + depth_index * 1000
                    + replicate
                )

                rng = np.random.default_rng(seed)

                thinned = thin_matrix(matrix, depth, rng)

                contacts = int(
                    thinned[
                        ROWS[LOCAL],
                        COLS[LOCAL],
                    ].sum()
                )

                record = {
                    "scenario": scenario,
                    "depth": depth,
                    "contacts": contacts,
                    "valid": False,
                    "covered": np.nan,
                    "singleton": np.nan,
                    "absolute_error": np.nan,
                }

                try:
                    result = bootstrap_w80(
                        thinned,
                        REGION,
                        max_separation_bins=MAX_DISTANCE,
                        n_bootstraps=N_BOOTSTRAPS,
                        seed=seed + 10000000,
                    )
                except ValueError:
                    records.append(record)
                    continue

                lower = result["bootstrap_lower"]
                upper = result["bootstrap_upper"]

                record.update({
                    "valid": True,
                    "covered": float(lower <= truth <= upper),
                    "singleton": float(lower == upper),
                    "absolute_error": abs(
                        result["observed_w80"] - truth
                    ),
                })

                records.append(record)

    df = pd.DataFrame(records)

    edges = [
        -0.5, 0.5, 5, 10, 25, 50,
        100, 250, 1000, 10000, np.inf
    ]

    labels = [
        "0", "1-4", "5-9", "10-24", "25-49",
        "50-99", "100-249", "250-999",
        "1000-9999", "10000+"
    ]

    df["contact_group"] = pd.cut(
        df["contacts"],
        bins=edges,
        labels=labels,
        right=False,
    )

    summary = (
        df.groupby(
            ["scenario", "contact_group"],
            observed=True,
            sort=False,
        )
        .agg(
            cases=("valid", "size"),
            valid=("valid", "sum"),
            coverage_valid=("covered", "mean"),
            singleton_valid=("singleton", "mean"),
            mean_absolute_error=("absolute_error", "mean"),
        )
        .reset_index()
    )

    print("W80 reliability by observed contact count")
    print(f"Replicates per scenario/depth: {N_REPLICATES}")
    print(f"Bootstrap resamples: {N_BOOTSTRAPS}")
    print("Nominal confidence level: 95%")
    print()
    print(summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    ))

    print()
    print("Contact-count calibration: COMPLETE")
    print("Coverage excludes invalid estimates.")
    print("Small groups have high sampling uncertainty.")
    print("No QC threshold is established.")


if __name__ == "__main__":
    main()
