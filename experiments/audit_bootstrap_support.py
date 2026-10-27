"""Diagnose distal-support limitations of the empirical W80 bootstrap.

Reconstruct the lowest-depth synthetic trials using original seeds.
No biological calibration is implied.
"""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map

trials = pd.read_csv(
    "research_artifacts/stage1/bootstrap_calibration_trials.csv"
)

rows, cols = np.triu_indices(256, k=1)
distances = cols - rows
midpoints = (rows + cols) / 2

local = (
    (midpoints >= 96)
    & (midpoints < 160)
    & (distances <= 64)
)

scenarios = [
    ("baseline", 0, 28),
    ("narrowing", 1, 9),
    ("coverage_drop", 2, 28),
]

print("BOOTSTRAP DISTAL-SUPPORT VALIDATION")
print("=" * 68)

for name, scenario_index, truth in scenarios:
    subset = trials[
        (trials.scenario == name)
        & np.isclose(trials.depth, 0.0001)
    ].sort_values("replicate")

    assert len(subset) == 100

    records = []

    for trial in subset.itertuples():
        rep = int(trial.replicate)

        matrix = simulate_contact_map(
            name,
            seed=42000 + scenario_index * 10000 + rep,
        )

        rng = np.random.default_rng(
            500000 + scenario_index * 100000 + 3000 + rep
        )

        sampled = rng.binomial(
            matrix[rows, cols].astype(np.int64),
            0.0001,
        )

        counts = sampled[local]
        d = distances[local]

        assert int(counts.sum()) == int(trial.contact_count)

        observed_distances = d[counts > 0]
        maximum = (
            int(observed_distances.max())
            if len(observed_distances)
            else 0
        )

        valid = bool(trial.valid)
        too_low = (
            valid and float(trial.upper_bins) < truth
        )

        records.append({
            "valid": valid,
            "too_low": too_low,
            "missing_distal": maximum < truth,
        })

    result = pd.DataFrame(records)
    valid = result[result.valid]
    failures = valid[valid.too_low]

    missing = int(valid.missing_distal.sum())
    low = len(failures)
    explained = int(failures.missing_distal.sum())

    print(f"\n{name}")
    print(f"  Valid intervals:                 {len(valid)}/100")
    print(f"  Missing distal support:          {missing}")
    print(f"  Intervals entirely below truth:  {low}")
    print(f"  Too-low with missing support:    {explained}")
    print(f"  Too-low despite distal support:  {low - explained}")

    assert (
        valid.loc[valid.missing_distal, "too_low"]
        .all()
    ), "Distal-support consistency failed"

print("\nPASS: All three scenarios reconstructed and audited")
