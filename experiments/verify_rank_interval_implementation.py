"""Verify reusable W80 rank intervals against prior calibration results."""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.uncertainty import (
    rank_interval_w80_equal_support,
)

calibration = pd.read_csv(
    "research_artifacts/stage1/bootstrap_calibration_trials.csv"
)

reference = pd.read_csv(
    "research_artifacts/stage1/w80_interval_method_comparison.csv"
)

assert len(calibration) == 1200
assert len(reference) == 1190

keys = ["scenario", "depth", "replicate"]

assert not calibration.duplicated(keys).any()
assert not reference.duplicated(keys).any()

calibration = calibration.set_index(keys)
reference = reference.set_index(keys)

scenarios = ("baseline", "narrowing", "coverage_drop")
depths = (1.0, 0.01, 0.001, 0.0001)

rows, cols = np.triu_indices(256, k=1)

verified = 0
invalid_verified = 0
mismatches = []

for scenario_index, scenario in enumerate(scenarios):
    for replicate in range(100):
        matrix = simulate_contact_map(
            scenario,
            seed=42000 + scenario_index * 10000 + replicate,
        )

        source_counts = matrix[rows, cols].astype(np.int64)

        for depth_index, depth in enumerate(depths):
            key = (scenario, depth, replicate)
            trial = calibration.loc[key]

            rng = np.random.default_rng(
                500000
                + scenario_index * 100000
                + depth_index * 1000
                + replicate
            )

            sampled = rng.binomial(source_counts, depth)

            thinned = np.zeros_like(matrix)
            thinned[rows, cols] = sampled
            thinned[cols, rows] = sampled

            if not bool(trial["valid"]):
                assert key not in reference.index
                assert int(trial["contact_count"]) == 0

                try:
                    rank_interval_w80_equal_support(
                        thinned,
                        (96, 160),
                        max_separation_bins=64,
                    )
                except ValueError:
                    invalid_verified += 1
                else:
                    raise AssertionError(
                        f"Expected zero-contact failure: {key}"
                    )
                continue

            expected = reference.loc[key]

            actual = rank_interval_w80_equal_support(
                thinned,
                (96, 160),
                max_separation_bins=64,
                confidence_level=0.95,
            )

            matches = (
                actual["rank_lower"] == int(expected["rank_lower"])
                and actual["rank_upper"] == int(expected["rank_upper"])
                and actual["total_contacts"] == int(trial["contact_count"])
                and actual["observed_w80"] == int(trial["observed_w80_bins"])
            )

            if not matches:
                mismatches.append({
                    "key": key,
                    "expected": (
                        int(expected["rank_lower"]),
                        int(expected["rank_upper"]),
                    ),
                    "actual": (
                        actual["rank_lower"],
                        actual["rank_upper"],
                    ),
                })

            verified += 1

assert verified == 1190
assert invalid_verified == 10

print("REUSABLE RANK-INTERVAL IMPLEMENTATION VERIFICATION")
print("=" * 72)
print(f"Valid trials checked: {verified}")
print(f"Zero-contact failures verified: {invalid_verified}")
print(f"Interval mismatches: {len(mismatches)}")

if mismatches:
    for mismatch in mismatches[:10]:
        print(mismatch)
    raise AssertionError("Rank-interval calibration mismatch")

print("PASS: All 1,190 interval bounds match previous calibration")
print("PASS: All 10 zero-contact cases rejected")
