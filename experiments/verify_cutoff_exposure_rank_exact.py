"""Exact finite-sample coverage audit for cutoff-specific W80 intervals.

Enumerates every possible contact-count allocation for small
four-band models, conditional on the total contact count.

Uses the cutoff-specific implementation from the existing
experiment without executing its Monte Carlo simulation.

Coverage calculations are exact up to numerical precision.
They do not establish validity for biological Hi-C data.
"""

import ast
from pathlib import Path

import numpy as np
from scipy.special import gammaln
from scipy.stats import binom

D_MAX = 4
ALPHA = 0.05
TARGET = 0.80

# Extract the actual interval function without executing
# the full simulation script at import time.
source_path = Path(
    "experiments/compare_cutoff_exposure_rank.py"
)

tree = ast.parse(source_path.read_text())

functions = [
    node for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name == "cutoff_specific_bounds"
]

assert len(functions) == 1

namespace = {
    "np": np,
    "binom": binom,
    "D_MAX": D_MAX,
    "ALPHA": ALPHA,
    "TARGET": TARGET,
}

module = ast.Module(
    body=functions,
    type_ignores=[],
)

exec(
    compile(module, str(source_path), "exec"),
    namespace,
)

interval = namespace["cutoff_specific_bounds"]

exposures = {
    "equal": [4, 4, 4, 4],
    "decreasing": [10, 7, 3, 1],
    "increasing": [1, 3, 7, 10],
    "irregular": [1, 10, 2, 8],
}

intensities = {
    "early_peak": [20, 1, 1, 1],
    "second_peak": [9, 2, 1, 1],
    "third_peak": [1, 2, 8, 1],
    "late_peak": [1, 1, 1, 9],
    "boundary_second": [4, 4, 1, 1],
    "boundary_third": [4, 1, 3, 2],
}

sample_sizes = (1, 3, 5, 10, 20, 40)


def enumerate_counts(n):
    """All four-band nonnegative count vectors summing to n."""
    return np.array([
        (a, b, c, n - a - b - c)
        for a in range(n + 1)
        for b in range(n - a + 1)
        for c in range(n - a - b + 1)
    ], dtype=np.int64)


print("EXACT CONDITIONAL COVERAGE AUDIT")
print("=" * 75)
print("Nominal confidence: 95%")
print("Distance bands: 4")
print("Every possible count allocation is enumerated.")
print()

results = []

for exposure_name, exposure_values in exposures.items():
    exposure = np.asarray(
        exposure_values, dtype=np.int64
    )

    for profile_name, intensity_values in intensities.items():
        intensity = np.asarray(
            intensity_values, dtype=float
        )

        population_w80 = int(
            np.searchsorted(
                np.cumsum(intensity),
                TARGET * intensity.sum(),
            ) + 1
        )

        probabilities = (
            exposure * intensity
            / np.sum(exposure * intensity)
        )

        for n in sample_sizes:
            count_vectors = enumerate_counts(n)

            # Exact multinomial probabilities.
            log_probabilities = (
                gammaln(n + 1)
                - gammaln(count_vectors + 1).sum(axis=1)
                + (
                    count_vectors
                    * np.log(probabilities)
                ).sum(axis=1)
            )

            probability_mass = np.exp(
                log_probabilities
            )

            assert abs(
                probability_mass.sum() - 1
            ) < 1e-10

            covered_probability = 0.0

            for counts, mass in zip(
                count_vectors, probability_mass
            ):
                lower, upper = interval(
                    counts, exposure
                )

                assert (
                    1 <= lower <= upper <= D_MAX
                )

                if (
                    lower <= population_w80 <= upper
                ):
                    covered_probability += mass

            results.append({
                "exposure": exposure_name,
                "profile": profile_name,
                "n": n,
                "true_w80": population_w80,
                "coverage": covered_probability,
            })

print(
    f"{'Exposure':<15} {'Cases':>7} "
    f"{'Minimum coverage':>19}"
)
print("-" * 48)

for exposure_name in exposures:
    group = [
        r["coverage"] for r in results
        if r["exposure"] == exposure_name
    ]

    print(
        f"{exposure_name:<15} {len(group):7d} "
        f"{100 * min(group):18.4f}%"
    )

worst = min(
    results, key=lambda r: r["coverage"]
)

print()
print(f"Total exact coverage evaluations: {len(results)}")
print(
    f"Worst-case evaluated coverage: "
    f"{100 * worst['coverage']:.4f}%"
)
print(
    f"Worst-case configuration: "
    f"{worst['exposure']}, {worst['profile']}, "
    f"n={worst['n']}"
)

failures = [
    r for r in results
    if r["coverage"] < 0.95 - 1e-10
]

print(
    f"Configurations below nominal 95%: "
    f"{len(failures)}"
)

assert not failures, (
    "Exact conditional coverage below nominal level"
)

print()
print("PASS: All enumerated models satisfy nominal coverage")
print("NOTE: Finite tested models do not constitute a full proof.")
print("NOTE: Independent Poisson sampling remains an assumption.")
