"""Coverage stress test for the provisional Hi-C W80 estimator.

Uses independent binomial thinning of synthetic contact counts.
Results measure robustness within the synthetic model only.
"""

import numpy as np

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.width import estimate_w80


N_REPLICATES = 50
REGION = (96, 160)
MAX_SEPARATION = 64
DEPTHS = (1.0, 0.1, 0.01, 0.001, 0.0001)
SCENARIOS = ("baseline", "narrowing", "coverage_drop")
SCORE_THRESHOLD = 0.5

N_BINS = 256
ROWS, COLS = np.triu_indices(N_BINS, k=1)

MIDPOINT = (ROWS + COLS) / 2
DISTANCE = COLS - ROWS

REGIONAL_PAIRS = (
    (MIDPOINT >= REGION[0])
    & (MIDPOINT < REGION[1])
    & (DISTANCE <= MAX_SEPARATION)
)


def thin_matrix(matrix, fraction, rng):
    """Binomially thin unique contacts, preserving symmetry."""
    counts = matrix[ROWS, COLS]

    sampled = rng.binomial(
        counts.astype(np.int64),
        fraction,
    )

    result = np.zeros_like(matrix)
    result[ROWS, COLS] = sampled
    result[COLS, ROWS] = sampled

    return result


def local_contact_count(matrix):
    """Count upper-triangle contacts used by the local analysis."""
    return int(matrix[ROWS[REGIONAL_PAIRS],
                      COLS[REGIONAL_PAIRS]].sum())


def main():
    records = []

    for replicate in range(N_REPLICATES):
        seed = 1000 + replicate

        baseline = simulate_contact_map(
            "baseline", seed=seed
        )

        reference_width = estimate_w80(
            baseline, REGION,
            max_separation_bins=MAX_SEPARATION,
        )

        for scenario_index, scenario in enumerate(SCENARIOS):
            matrix = simulate_contact_map(
                scenario, seed=seed
            )

            for depth_index, fraction in enumerate(DEPTHS):
                rng = np.random.default_rng(
                    1000000
                    + replicate * 1000
                    + scenario_index * 100
                    + depth_index
                )

                thinned = thin_matrix(matrix, fraction, rng)
                contacts = local_contact_count(thinned)

                try:
                    width = estimate_w80(
                        thinned,
                        REGION,
                        max_separation_bins=MAX_SEPARATION,
                    )
                    score = float(
                        np.log2(reference_width / width)
                    )
                except ValueError:
                    width = np.nan
                    score = np.nan

                records.append(
                    (scenario, fraction, contacts, width, score)
                )

    print(f"Replicates per scenario/depth: {N_REPLICATES}")
    print(f"Exploratory narrowing threshold: {SCORE_THRESHOLD}")
    print("Reference: full-depth baseline from matching seed")
    print()

    header = (
        f"{'Scenario':15s} {'Depth':>8s} "
        f"{'Median contacts':>16s} {'Valid':>7s} "
        f"{'Median W80':>12s} {'Calls':>9s}"
    )

    print(header)
    print("-" * len(header))

    for scenario in SCENARIOS:
        for fraction in DEPTHS:
            subset = [
                r for r in records
                if r[0] == scenario and r[1] == fraction
            ]

            contacts = np.array([r[2] for r in subset])
            widths = np.array([r[3] for r in subset])
            scores = np.array([r[4] for r in subset])

            valid = np.isfinite(widths)
            n_valid = int(valid.sum())

            median_width = (
                f"{np.median(widths[valid]):.1f}"
                if n_valid else "NA"
            )

            n_calls = int(
                np.sum(scores[valid] > SCORE_THRESHOLD)
            )

            print(
                f"{scenario:15s} {fraction:8.4f} "
                f"{np.median(contacts):16.0f} "
                f"{n_valid:3d}/{N_REPLICATES:<3d} "
                f"{median_width:>12s} "
                f"{n_calls:3d}/{N_REPLICATES:<3d}"
            )

        print()

    print("Coverage stress test: COMPLETE")
    print("Note: Calls use an ideal full-depth baseline reference.")
    print("No coverage-quality threshold has been validated.")


if __name__ == "__main__":
    main()
