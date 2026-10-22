"""Multi-seed evaluation of the provisional W80 estimator.

Synthetic experiments only; not biological validation.
"""

import numpy as np

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.width import estimate_w80


SCENARIOS = ("baseline", "narrowing", "coverage_drop")
N_REPLICATES = 100
REGION = (96, 160)
MAX_SEPARATION = 64


def main():
    results = {scenario: [] for scenario in SCENARIOS}

    for i in range(N_REPLICATES):
        for j, scenario in enumerate(SCENARIOS):
            # Separate random streams for each experimental condition.
            seed = 1000 + i + j * 10000

            matrix = simulate_contact_map(
                scenario,
                seed=seed,
            )

            width = estimate_w80(
                matrix,
                region=REGION,
                max_separation_bins=MAX_SEPARATION,
            )

            results[scenario].append(width)

    for scenario in SCENARIOS:
        results[scenario] = np.asarray(results[scenario])

    print(f"Replicates per scenario: {N_REPLICATES}")
    print(f"Region: {REGION}")
    print(f"Maximum separation: {MAX_SEPARATION} bins")
    print()

    for scenario, widths in results.items():
        print(
            f"{scenario:15s} "
            f"mean={widths.mean():.2f}, "
            f"std={widths.std(ddof=1):.2f}, "
            f"median={np.median(widths):.1f}, "
            f"range=[{widths.min()}, {widths.max()}], "
            f"2.5-97.5%=[{np.quantile(widths, 0.025):.1f}, "
            f"{np.quantile(widths, 0.975):.1f}]"
        )

    baseline = results["baseline"]

    # Exploratory threshold, not a validated detection rule.
    threshold = 0.5

    for scenario in ("narrowing", "coverage_drop"):
        scores = np.log2(baseline / results[scenario])

        print()
        print(f"Scenario: {scenario}")
        print(f"Mean narrowing score: {scores.mean():.4f}")
        print(
            "Fraction above exploratory threshold: "
            f"{np.mean(scores > threshold):.3f}"
        )

    print()
    print("Multi-seed W80 evaluation: COMPLETE")


if __name__ == "__main__":
    main()
