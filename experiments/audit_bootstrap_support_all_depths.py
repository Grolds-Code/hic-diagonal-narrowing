"""Audit distal-support explanations across all bootstrap trials."""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map

df = pd.read_csv(
    "research_artifacts/stage1/bootstrap_calibration_trials.csv"
)

assert len(df) == 1200

scenarios = [
    ("baseline", 28),
    ("narrowing", 9),
    ("coverage_drop", 28),
]

depths = (1.0, 0.01, 0.001, 0.0001)

rows, cols = np.triu_indices(256, k=1)
distances = cols - rows
midpoints = (rows + cols) / 2

local = (
    (midpoints >= 96)
    & (midpoints < 160)
    & (distances <= 64)
)

print("BOOTSTRAP SUPPORT — ALL DEPTHS")
print("=" * 85)
print(
    f"{'Scenario':<15} {'Depth':>8} {'Valid':>7} "
    f"{'Too low':>9} {'Missing':>9} "
    f"{'Explained':>10} {'Unexplained':>12}"
)

total_low = 0
total_explained = 0

for si, (scenario, truth) in enumerate(scenarios):
    for di, depth in enumerate(depths):
        subset = df[
            (df.scenario == scenario)
            & np.isclose(df.depth, depth)
        ]

        assert len(subset) == 100

        valid_count = 0
        too_low_count = 0
        missing_count = 0
        explained_count = 0

        for trial in subset.itertuples():
            rep = int(trial.replicate)

            matrix = simulate_contact_map(
                scenario,
                seed=42000 + si * 10000 + rep,
            )

            rng = np.random.default_rng(
                500000 + si * 100000 + di * 1000 + rep
            )

            sampled = rng.binomial(
                matrix[rows, cols].astype(np.int64),
                depth,
            )

            counts = sampled[local]
            assert int(counts.sum()) == int(trial.contact_count)

            observed = distances[local][counts > 0]
            max_distance = int(observed.max()) if observed.size else 0

            if not bool(trial.valid):
                continue

            valid_count += 1

            missing = max_distance < truth
            too_low = float(trial.upper_bins) < truth

            missing_count += int(missing)
            too_low_count += int(too_low)
            explained_count += int(missing and too_low)

            # If all observed contacts are closer than the
            # population width, empirical bootstrap samples
            # cannot extend beyond those observed separations.
            if missing:
                assert too_low

        unexplained = too_low_count - explained_count
        total_low += too_low_count
        total_explained += explained_count

        print(
            f"{scenario:<15} {depth:8.4f} "
            f"{valid_count:7d} {too_low_count:9d} "
            f"{missing_count:9d} "
            f"{explained_count:10d} {unexplained:12d}"
        )

print("\nTotal downward failures:", total_low)
print("Failures explained by missing distal support:", total_explained)
print("Failures despite distal support:", total_low - total_explained)
print("\nPASS: All 1,200 trials audited")
