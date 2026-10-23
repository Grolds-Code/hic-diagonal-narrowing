"""Evaluate narrowing detection across genomic positions.

Synthetic benchmark only.
Detection parameters are fixed before evaluation.
"""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.scan import scan_local_width
from hic_diagonal_narrowing.detect import detect_narrowing


RESOLUTION = 500
SEEDS = range(20)

REGIONS = [
    (32, 96),
    (64, 128),
    (96, 160),
    (128, 192),
    (160, 224),
]

SCENARIOS = (
    "baseline",
    "narrowing",
    "coverage_drop",
)


def interval_iou(start_a, end_a, start_b, end_b):
    """Intersection-over-union of two genomic intervals."""
    intersection = max(
        0,
        min(end_a, end_b) - max(start_a, start_b),
    )

    union = (
        (end_a - start_a)
        + (end_b - start_b)
        - intersection
    )

    return intersection / union if union > 0 else 0.0


def main():
    records = []

    for region in REGIONS:
        truth_start = region[0] * RESOLUTION
        truth_end = region[1] * RESOLUTION

        for replicate in SEEDS:
            seed = 10000 + replicate

            for scenario in SCENARIOS:
                matrix = simulate_contact_map(
                    scenario,
                    region=region,
                    seed=seed,
                )

                profile = scan_local_width(
                    matrix,
                    chromosome="chrSynthetic",
                    resolution_bp=RESOLUTION,
                    window_bins=64,
                    step_bins=16,
                    max_separation_bins=64,
                )

                _, candidates = detect_narrowing(
                    profile,
                    score_threshold=0.5,
                    min_reference_windows=3,
                )

                overlaps = []
                peak_overlaps = []

                for candidate in candidates.itertuples():
                    overlaps.append(
                        interval_iou(
                            candidate.start_bp,
                            candidate.end_bp,
                            truth_start,
                            truth_end,
                        )
                    )

                    peak_overlaps.append(
                        interval_iou(
                            candidate.peak_start_bp,
                            candidate.peak_end_bp,
                            truth_start,
                            truth_end,
                        )
                    )

                records.append({
                    "scenario": scenario,
                    "region": str(region),
                    "replicate": replicate,
                    "n_candidates": len(candidates),
                    "any_call": int(len(candidates) > 0),
                    "event_detected": int(
                        any(value > 0 for value in overlaps)
                    ),
                    "best_footprint_iou": max(
                        overlaps, default=0.0
                    ),
                    "best_peak_iou": max(
                        peak_overlaps, default=0.0
                    ),
                })

    df = pd.DataFrame(records)

    summary = (
        df.groupby(
            ["scenario", "region"],
            sort=False,
        )
        .agg(
            trials=("replicate", "size"),
            call_rate=("any_call", "mean"),
            event_detection_rate=("event_detected", "mean"),
            mean_candidates=("n_candidates", "mean"),
            mean_footprint_iou=(
                "best_footprint_iou", "mean"
            ),
            mean_peak_iou=("best_peak_iou", "mean"),
        )
        .reset_index()
    )

    print("MULTI-POSITION DETECTOR EVALUATION")
    print("=" * 70)
    print(f"Genomic positions: {len(REGIONS)}")
    print(f"Replicates per position: {len(SEEDS)}")
    print(f"Total simulations: {len(df)}")
    print("Fixed exploratory score threshold: 0.5")
    print()

    print(summary.to_string(
        index=False,
        float_format=lambda value: f"{value:.3f}",
    ))

    print()
    print("OVERALL SUMMARY")

    for scenario in SCENARIOS:
        subset = df[df["scenario"] == scenario]

        print(
            f"{scenario:15s} | "
            f"any-call rate={subset['any_call'].mean():.3f} | "
            f"event detection={subset['event_detected'].mean():.3f} | "
            f"mean footprint IoU="
            f"{subset['best_footprint_iou'].mean():.3f} | "
            f"mean peak IoU="
            f"{subset['best_peak_iou'].mean():.3f}"
        )

    assert len(df) == 300

    print()
    print("Multi-position evaluation: COMPLETE")
    print("No detector parameters were optimized in this test.")
    print("All results are synthetic, not biological validation.")


if __name__ == "__main__":
    main()
