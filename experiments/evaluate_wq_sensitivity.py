"""Multi-seed sensitivity of cumulative Hi-C width thresholds."""

import numpy as np

from hic_diagonal_narrowing.synthetic import simulate_contact_map

thresholds = [0.50, 0.70, 0.80, 0.81, 0.90, 0.95]
scenarios = ["baseline", "narrowing", "coverage_drop"]
region = (96, 160)
results = {}

for scenario in scenarios:
    widths = []

    for seed in range(100):
        matrix = simulate_contact_map(scenario, seed=seed)
        profile = []

        for d in range(1, 65):
            rows = np.arange(matrix.shape[0] - d)
            mid = rows + d / 2
            valid = (mid >= region[0]) & (mid < region[1])

            profile.append(
                matrix[rows[valid], rows[valid] + d].mean()
            )

        profile = np.asarray(profile)
        cumulative = np.cumsum(profile) / profile.sum()

        widths.append([
            np.searchsorted(cumulative, q) + 1
            for q in thresholds
        ])

    results[scenario] = np.asarray(widths)

print("\nWq MULTI-SEED SENSITIVITY — 100 REPLICATES")
print("=" * 68)

for i, q in enumerate(thresholds):
    print(f"\nW{int(q * 100)}")
    for scenario in scenarios:
        values = results[scenario][:, i]
        print(
            f"  {scenario:<15}"
            f"mean={values.mean():6.2f}  "
            f"std={values.std(ddof=1):5.2f}"
        )

agreement = np.mean(
    results["baseline"][:, 2] == results["baseline"][:, 3]
)

print("\nW80/W81 agreement in baseline simulations:", f"{agreement:.1%}")
print("PASS: Multi-seed sensitivity analysis completed")
