"""Evaluate QC metrics against W80 reliability.

Synthetic calibration experiment. Full-depth widths serve as
simulation references, not biological ground truth.
"""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.width import estimate_w80
from hic_diagonal_narrowing.qc import summarize_local_support


REGION = (96, 160)
MAX_DISTANCE = 64
REPLICATES = 50
DEPTHS = (1.0, 0.1, 0.01, 0.001, 0.0001)
SCENARIOS = ("baseline", "narrowing", "coverage_drop")

ROWS, COLS = np.triu_indices(256, k=1)


def thin_contacts(matrix, depth, rng):
    """Independently retain contacts while preserving symmetry."""
    sampled = rng.binomial(
        matrix[ROWS, COLS].astype(np.int64),
        depth,
    )

    result = np.zeros_like(matrix)
    result[ROWS, COLS] = sampled
    result[COLS, ROWS] = sampled

    return result


def main():
    records = []

    for replicate in range(REPLICATES):
        seed = 5000 + replicate

        baseline = simulate_contact_map(
            "baseline", seed=seed
        )
        baseline_width = estimate_w80(
            baseline, REGION,
            max_separation_bins=MAX_DISTANCE,
        )

        for scenario_idx, scenario in enumerate(SCENARIOS):
            full_matrix = simulate_contact_map(
                scenario, seed=seed
            )

            full_width = estimate_w80(
                full_matrix, REGION,
                max_separation_bins=MAX_DISTANCE,
            )

            for depth_idx, depth in enumerate(DEPTHS):
                rng = np.random.default_rng(
                    100000 + replicate * 100
                    + scenario_idx * 10 + depth_idx
                )

                matrix = thin_contacts(
                    full_matrix, depth, rng
                )

                qc = summarize_local_support(
                    matrix,
                    REGION,
                    max_separation_bins=MAX_DISTANCE,
                )

                try:
                    width = estimate_w80(
                        matrix,
                        REGION,
                        max_separation_bins=MAX_DISTANCE,
                    )
                    error = abs(width - full_width)
                    score = np.log2(baseline_width / width)
                except ValueError:
                    width = np.nan
                    error = np.nan
                    score = np.nan

                records.append({
                    "scenario": scenario,
                    "depth": depth,
                    "contacts": qc["total_contact_signal"],
                    "occupancy": qc["pair_occupancy"],
                    "distance_support": qc[
                        "distance_support_fraction"
                    ],
                    "width": width,
                    "absolute_error": error,
                    "narrowing_call": (
                        float(score > 0.5)
                        if np.isfinite(score)
                        else np.nan
                    ),
                })

    df = pd.DataFrame(records)

    print("QC reliability experiment")
    print(f"Replicates: {REPLICATES}")
    print("W80 error reference: same-scenario full-depth width")
    print("Narrowing-call reference: full-depth baseline")
    print("Threshold: log2 narrowing score > 0.5")
    print()

    summary = (
        df.groupby(["scenario", "depth"], sort=False)
        .agg(
            median_contacts=("contacts", "median"),
            median_occupancy=("occupancy", "median"),
            median_distance_support=(
                "distance_support", "median"
            ),
            valid_widths=("width", "count"),
            mean_absolute_error=("absolute_error", "mean"),
            call_rate_valid=("narrowing_call", "mean"),
        )
        .reset_index()
    )

    print(summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}",
    ))

    print()
    print("QC reliability evaluation: COMPLETE")
    print("No QC acceptance threshold has been established.")


if __name__ == "__main__":
    main()
