"""Calibrate Poisson W80 intervals with unequal distance support.

Independent Poisson band counts with fixed eligible-pair exposures.
These are controlled statistical simulations, not biological data.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.poisson_uncertainty import (
    poisson_w80_interval_from_bands,
)

N_BINS = 256
REGION = (96, 160)
MAX_DISTANCE = 64
N_REPLICATES = 200

DEPTHS = (1.0, 0.01, 0.001, 0.0001)
SCENARIOS = ("baseline", "narrowing", "coverage_drop")

distances = np.arange(1, MAX_DISTANCE + 1)

# Construct a fixed, reproducible missing-bin pattern.
mask_rng = np.random.default_rng(20261028)

valid_bins = np.ones(N_BINS, dtype=bool)
invalid_indices = mask_rng.choice(
    N_BINS, size=58, replace=False
)
valid_bins[invalid_indices] = False

eligible = []

for distance in distances:
    rows = np.arange(N_BINS - distance)
    midpoints = rows + distance / 2

    selected = (
        (midpoints >= REGION[0])
        & (midpoints < REGION[1])
        & valid_bins[rows]
        & valid_bins[rows + distance]
    )

    eligible.append(int(selected.sum()))

eligible = np.asarray(eligible, dtype=np.int64)

assert np.all(eligible > 0)
assert len(np.unique(eligible)) > 1

baseline = 90.0 * np.exp(-distances / 18.0)
narrowing = 90.0 * np.exp(-distances / 5.0)

# Preserve expected regional contacts for narrowing.
narrowing *= (
    np.sum(eligible * baseline)
    / np.sum(eligible * narrowing)
)

intensities = {
    "baseline": baseline,
    "narrowing": narrowing,
    "coverage_drop": 0.25 * baseline,
}

def population_w80(intensity):
    cumulative = np.cumsum(intensity)
    return int(
        np.searchsorted(
            cumulative,
            0.80 * intensity.sum(),
        ) + 1
    )

truths = {
    scenario: population_w80(intensity)
    for scenario, intensity in intensities.items()
}

assert truths == {
    "baseline": 28,
    "narrowing": 9,
    "coverage_drop": 28,
}

rng = np.random.default_rng(810241)
records = []

for scenario in SCENARIOS:
    truth = truths[scenario]
    intensity = intensities[scenario]

    for depth in DEPTHS:
        expected_counts = eligible * intensity * depth

        for replicate in range(N_REPLICATES):
            counts = rng.poisson(expected_counts).astype(
                np.int64
            )

            total = int(counts.sum())

            record = {
                "scenario": scenario,
                "depth": depth,
                "replicate": replicate,
                "population_w80": truth,
                "contacts": total,
                "valid": False,
                "covered": False,
                "observed_w80": np.nan,
                "lower": np.nan,
                "upper": np.nan,
                "span": np.nan,
            }

            if total == 0:
                records.append(record)
                continue

            result = poisson_w80_interval_from_bands(
                counts,
                eligible,
                confidence_level=0.95,
            )

            lower = result["poisson_lower"]
            upper = result["poisson_upper"]

            record.update({
                "valid": True,
                "covered": lower <= truth <= upper,
                "observed_w80": result["observed_w80"],
                "lower": lower,
                "upper": upper,
                "span": upper - lower,
            })

            records.append(record)

df = pd.DataFrame(records)

expected_trials = (
    len(SCENARIOS) * len(DEPTHS) * N_REPLICATES
)

assert len(df) == expected_trials

output = Path(
    "research_artifacts/stage1/"
    "poisson_unequal_support_calibration.csv"
)
output.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output, index=False)

print("POISSON W80 — UNEQUAL-SUPPORT CALIBRATION")
print("=" * 99)
print(f"Genomic bins: {N_BINS}")
print(f"Invalid genomic bins: {(~valid_bins).sum()}")
print(f"Distance bands: {MAX_DISTANCE}")
print(
    f"Eligible pairs per band: "
    f"min={eligible.min()}, max={eligible.max()}"
)
print(f"Unique eligible-pair counts: {len(np.unique(eligible))}")
print(f"Replicates per condition: {N_REPLICATES}")
print("Nominal confidence level: 95%")
print()

print(
    f"{'Scenario':<15} {'Depth':>8} {'Truth':>6} "
    f"{'Valid':>9} {'Coverage*':>11} "
    f"{'All trials':>11} {'Contacts':>10} "
    f"{'Med span':>10}"
)
print("-" * 99)

for (scenario, depth), group in df.groupby(
    ["scenario", "depth"], sort=False
):
    valid = group[group["valid"]]

    n_valid = len(valid)
    covered = int(valid["covered"].sum())

    conditional_coverage = (
        f"{100 * covered / n_valid:.1f}%"
        if n_valid else "NA"
    )

    all_trial_coverage = (
        f"{100 * covered / len(group):.1f}%"
    )

    median_span = (
        f"{valid['span'].median():.2f}"
        if n_valid else "NA"
    )

    print(
        f"{scenario:<15} {depth:8.4f} "
        f"{truths[scenario]:6d} "
        f"{n_valid:3d}/{N_REPLICATES:<5d} "
        f"{conditional_coverage:>11} "
        f"{all_trial_coverage:>11} "
        f"{group['contacts'].median():10.1f} "
        f"{median_span:>10}"
    )

print()
print(f"Output: {output}")
print("*Coverage conditional on a nonzero contact count.")
print("All trials counts abstentions as non-covering.")
print("Interval spans are in genomic bins.")
print("Assumes independent Poisson counts and fixed exposures.")
print("No biological Hi-C calibration is established.")
print("PASS: Unequal-support calibration completed")
