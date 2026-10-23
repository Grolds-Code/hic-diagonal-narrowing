"""Evaluate W80 detection on off-grid, variable-width events.

The detection configuration is fixed before evaluation.
All results are synthetic and exploratory.
"""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.synthetic import simulate_contact_map
from hic_diagonal_narrowing.scan import scan_local_width
from hic_diagonal_narrowing.detect import detect_narrowing


RESOLUTION_BP = 500
N_REPLICATES = 12

# Every start is deliberately off the 16-bin scanning grid.
REGIONS = [
    (19, 43),     # 24 bins
    (45, 77),     # 32 bins
    (91, 131),    # 40 bins
    (139, 187),   # 48 bins
    (73, 129),    # 56 bins
    (105, 169),   # 64 bins
    (27, 99),     # 72 bins
    (119, 207),   # 88 bins
    (151, 247),   # 96 bins
]

SCENARIOS = ("baseline", "narrowing", "coverage_drop")


def interval_iou(a_start, a_end, b_start, b_end):
    intersection = max(
        0,
        min(a_end, b_end) - max(a_start, b_start),
    )

    union = (
        (a_end - a_start)
        + (b_end - b_start)
        - intersection
    )

    return intersection / union if union > 0 else 0.0


def main():
    records = []

    assert all(
        start % 16 != 0
        for start, _ in REGIONS
    )

    for region in REGIONS:
        truth_start = region[0] * RESOLUTION_BP
        truth_end = region[1] * RESOLUTION_BP
        width_bins = region[1] - region[0]

        for replicate in range(N_REPLICATES):
            seed = 20000 + replicate

            for scenario in SCENARIOS:
                matrix = simulate_contact_map(
                    scenario,
                    region=region,
                    seed=seed,
                )

                profile = scan_local_width(
                    matrix,
                    chromosome="chrSynthetic",
                    resolution_bp=RESOLUTION_BP,
                    window_bins=64,
                    step_bins=16,
                    max_separation_bins=64,
                )

                scored, candidates = detect_narrowing(
                    profile,
                    score_threshold=0.5,
                    min_reference_windows=3,
                )

                footprint_ious = []
                peak_ious = []

                for candidate in candidates.itertuples():
                    footprint_ious.append(
                        interval_iou(
                            candidate.start_bp,
                            candidate.end_bp,
                            truth_start,
                            truth_end,
                        )
                    )

                    peak_ious.append(
                        interval_iou(
                            candidate.peak_start_bp,
                            candidate.peak_end_bp,
                            truth_start,
                            truth_end,
                        )
                    )

                best_footprint = max(
                    footprint_ious, default=0.0
                )
                best_peak = max(
                    peak_ious, default=0.0
                )

                records.append({
                    "scenario": scenario,
                    "region": str(region),
                    "width_bins": width_bins,
                    "replicate": replicate,
                    "n_candidates": len(candidates),
                    "any_call": int(len(candidates) > 0),
                    "any_overlap": int(
                        best_footprint > 0
                    ),
                    "peak_iou_ge_025": int(
                        best_peak >= 0.25
                    ),
                    "best_footprint_iou": best_footprint,
                    "best_peak_iou": best_peak,
                    "scored_windows": int(
                        scored["narrowing_score"]
                        .notna().sum()
                    ),
                })

    df = pd.DataFrame(records)

    narrowing = df[
        df["scenario"] == "narrowing"
    ]

    summary = (
        narrowing.groupby(
            ["region", "width_bins"],
            sort=False,
        )
        .agg(
            trials=("replicate", "size"),
            call_rate=("any_call", "mean"),
            overlap_rate=("any_overlap", "mean"),
            peak_iou_ge_025=(
                "peak_iou_ge_025", "mean"
            ),
            mean_footprint_iou=(
                "best_footprint_iou", "mean"
            ),
            mean_peak_iou=(
                "best_peak_iou", "mean"
            ),
            mean_candidates=("n_candidates", "mean"),
        )
        .reset_index()
    )

    print("OFF-GRID NARROWING BENCHMARK")
    print("=" * 72)
    print(f"Regions: {len(REGIONS)}")
    print(f"Replicates per scenario/region: {N_REPLICATES}")
    print(f"Total simulations: {len(df)}")
    print("Scanning window: 64 bins")
    print("Scanning step: 16 bins")
    print("Exploratory score threshold: 0.5")
    print()

    print("NARROWING SCENARIOS")
    print(summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    ))

    print()
    print("OVERALL RESULTS")

    for scenario in SCENARIOS:
        subset = df[
            df["scenario"] == scenario
        ]

        print(
            f"{scenario:15s} | "
            f"call rate={subset['any_call'].mean():.3f} | "
            f"overlap={subset['any_overlap'].mean():.3f} | "
            f"peak IoU>=0.25="
            f"{subset['peak_iou_ge_025'].mean():.3f} | "
            f"mean peak IoU="
            f"{subset['best_peak_iou'].mean():.3f} | "
            f"mean candidates="
            f"{subset['n_candidates'].mean():.3f}"
        )

    assert len(df) == (
        len(REGIONS) * N_REPLICATES * len(SCENARIOS)
    )

    print()
    print("Off-grid benchmark: COMPLETE")
    print("No parameters were optimized using these results.")
    print("Synthetic performance is not biological validation.")


if __name__ == "__main__":
    main()
