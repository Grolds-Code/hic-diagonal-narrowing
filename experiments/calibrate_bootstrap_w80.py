"""Evaluate bootstrap W80 interval calibration.

Population targets are derived from the synthetic generative model.
Results are not claims about biological Hi-C confidence intervals.
"""

import numpy as np

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.uncertainty import bootstrap_w80


REGION = (96, 160)
MAX_DISTANCE = 64
N_REPLICATES = 100
N_BOOTSTRAPS = 300

DEPTHS = (1.0, 0.01, 0.001, 0.0001)
SCENARIOS = ("baseline", "narrowing", "coverage_drop")

ROWS, COLS = np.triu_indices(256, k=1)

MIDPOINTS = (ROWS + COLS) / 2
DISTANCES = COLS - ROWS

LOCAL_PAIRS = (
    (MIDPOINTS >= REGION[0])
    & (MIDPOINTS < REGION[1])
    & (DISTANCES <= MAX_DISTANCE)
)


def population_w80(scenario):
    """Exact W80 of the model's expected distance profile."""
    scale = 5.0 if scenario == "narrowing" else 18.0

    distances = np.arange(1, MAX_DISTANCE + 1)
    intensity = np.exp(-distances / scale)

    cumulative = np.cumsum(intensity)
    target = 0.8 * intensity.sum()

    return int(np.searchsorted(cumulative, target) + 1)


def thin_matrix(matrix, depth, rng):
    counts = matrix[ROWS, COLS].astype(np.int64)
    sampled = rng.binomial(counts, depth)

    result = np.zeros_like(matrix)
    result[ROWS, COLS] = sampled
    result[COLS, ROWS] = sampled

    return result


def main():
    assert population_w80("baseline") == 28
    assert population_w80("narrowing") == 9
    assert population_w80("coverage_drop") == 28

    print("Bootstrap W80 calibration")
    print(f"Independent synthetic replicates: {N_REPLICATES}")
    print(f"Bootstrap resamples per estimate: {N_BOOTSTRAPS}")
    print("Nominal confidence level: 95%")
    print()

    header = (
        f"{'Scenario':15s} {'Depth':>8s} "
        f"{'True':>5s} {'Contacts':>9s} "
        f"{'Valid':>8s} {'Coverage':>10s} "
        f"{'Singleton':>10s} {'MAE':>8s}"
    )

    print(header)
    print("-" * len(header))

    for scenario_index, scenario in enumerate(SCENARIOS):
        truth = population_w80(scenario)

        results = {
            depth: {
                "contacts": [],
                "covered": [],
                "singleton": [],
                "errors": [],
            }
            for depth in DEPTHS
        }

        for replicate in range(N_REPLICATES):
            matrix = simulate_contact_map(
                scenario,
                seed=42000 + scenario_index * 10000 + replicate,
            )

            for depth_index, depth in enumerate(DEPTHS):
                rng = np.random.default_rng(
                    500000
                    + scenario_index * 100000
                    + depth_index * 1000
                    + replicate
                )

                thinned = thin_matrix(matrix, depth, rng)

                record = results[depth]

                contacts = int(
                    thinned[
                        ROWS[LOCAL_PAIRS],
                        COLS[LOCAL_PAIRS],
                    ].sum()
                )

                record["contacts"].append(contacts)

                try:
                    estimate = bootstrap_w80(
                        thinned,
                        REGION,
                        max_separation_bins=MAX_DISTANCE,
                        n_bootstraps=N_BOOTSTRAPS,
                        seed=900000
                        + scenario_index * 100000
                        + depth_index * 1000
                        + replicate,
                    )
                except ValueError:
                    continue

                lower = estimate["bootstrap_lower"]
                upper = estimate["bootstrap_upper"]
                observed = estimate["observed_w80"]

                record["covered"].append(
                    lower <= truth <= upper
                )
                record["singleton"].append(
                    lower == upper
                )
                record["errors"].append(
                    abs(observed - truth)
                )

        for depth in DEPTHS:
            record = results[depth]
            valid = len(record["covered"])

            coverage = (
                f"{100 * np.mean(record['covered']):.1f}%"
                if valid else "NA"
            )

            singleton = (
                f"{100 * np.mean(record['singleton']):.1f}%"
                if valid else "NA"
            )

            mae = (
                f"{np.mean(record['errors']):.2f}"
                if valid else "NA"
            )

            print(
                f"{scenario:15s} {depth:8.4f} "
                f"{truth:5d} "
                f"{np.median(record['contacts']):9.0f} "
                f"{valid:3d}/{N_REPLICATES:<4d} "
                f"{coverage:>10s} "
                f"{singleton:>10s} "
                f"{mae:>8s}"
            )

        print()

    print("Bootstrap calibration experiment: COMPLETE")
    print("Coverage and singleton rates exclude invalid estimates.")
    print("Failures must be assessed separately.")
    print("No confidence interval is biologically validated.")


if __name__ == "__main__":
    main()
