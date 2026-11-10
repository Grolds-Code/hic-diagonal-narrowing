"""Provisional coverage-aware Stage 1 detection comparison."""

from pathlib import Path

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.detect import detect_narrowing

root = Path("research_artifacts/stage1")

profile = pd.read_csv(
    root / "hct116_stage1_500bp_step_profile.csv"
)
original = pd.read_csv(
    root / "hct116_stage1_500bp_step_scored.csv"
)
coverage = pd.read_csv(
    root / "hct116_stage1_low_support_full_scan_v1.csv"
)
regions = pd.read_csv(
    root / "hct116_stage1_candidate_deliverable_v1.csv"
)

output_windows = (
    root / "hct116_stage1_coverage_triage_v1.csv"
)
output_regions = (
    root / "hct116_stage1_coverage_filtered_preview_regions_v1.csv"
)

if output_windows.exists() or output_regions.exists():
    raise RuntimeError("Refusing to overwrite previous output")

if (len(profile), len(original), len(coverage), len(regions)) != (
    1937, 1937, 1809, 5
):
    raise RuntimeError("Unexpected source dimensions")

keys = ["chromosome", "start_bp", "end_bp"]

if not profile[keys].equals(original[keys]):
    raise RuntimeError("Original profiles do not align")

metrics = [
    "start_bp",
    "end_bp",
    "endpoint_longest_zero_run",
    "endpoint_low10_longest_run",
    "endpoint_low25_longest_run",
    "endpoint_p10_to_median_ratio",
    "endpoint_zero_bins",
]

joined = original.merge(
    coverage[metrics],
    on=["start_bp", "end_bp"],
    how="left",
    validate="one_to_one",
    sort=False,
)

if len(joined) != 1937:
    raise RuntimeError("Incomplete coverage integration")

interior = joined["endpoint_longest_zero_run"].notna()

if int(interior.sum()) != 1809:
    raise RuntimeError("Interior geometry mismatch")

zero = joined["endpoint_longest_zero_run"].ge(2)
low10 = joined["endpoint_low10_longest_run"].ge(2)
low25 = joined["endpoint_low25_longest_run"].ge(2)

# Ordered classifications: the strongest applicable
# warning takes precedence over weaker warnings.
flags = np.select(
    [
        ~interior,
        zero,
        low10,
        low25,
    ],
    [
        "boundary_unassessed",
        "consecutive_zero_support",
        "severe_low_support",
        "moderate_low_support",
    ],
    default="no_coverage_warning",
)

eligible = flags == "no_coverage_warning"

screened = profile.copy()

screened["qc_status"] = np.where(
    eligible,
    "unassessed",
    "coverage_warning",
)

rescored, preview_regions = detect_narrowing(
    screened,
    score_threshold=0.20,
    min_reference_windows=3,
)

original_calls = (
    original["candidate_call"]
    .astype(str)
    .str.lower()
    .eq("true")
    .to_numpy()
)

new_calls = (
    rescored["candidate_call"].to_numpy(dtype=bool)
)

if (
    int(original_calls.sum()) != 93
    or int(new_calls.sum()) != 25
    or len(preview_regions) != 2
    or np.any(new_calls & ~eligible)
):
    raise RuntimeError(
        "Filtered preview disagrees with Step 381"
    )

triage = joined.copy()

triage["coverage_flag"] = flags
triage["provisionally_eligible"] = eligible
triage["original_candidate_call"] = original_calls
triage["preview_candidate_call"] = new_calls
triage["preview_reference_w80_bins"] = (
    rescored["reference_w80_bins"].to_numpy()
)
triage["preview_narrowing_score"] = (
    rescored["narrowing_score"].to_numpy()
)

# A lack of a warning is not a verified QC pass.
triage["formal_qc_status"] = "unassessed"

for region in regions.itertuples(index=False):
    subset = triage.loc[
        original_calls
        & (
            triage["start_bp"]
            >= region.footprint_start_bp
        )
        & (
            triage["end_bp"]
            <= region.footprint_end_bp
        )
    ]

    if len(subset) != int(region.supporting_windows):
        raise RuntimeError(
            f"{region.candidate_id}: membership mismatch"
        )

    print(f"\n{region.candidate_id}")
    print("-" * 60)
    print("Original calls:", len(subset))
    print(
        "Coverage classifications:",
        subset["coverage_flag"].value_counts().to_dict()
    )
    print(
        "Retained provisional calls:",
        int(subset["preview_candidate_call"].sum())
    )

triage.to_csv(output_windows, index=False)

preview_regions = preview_regions.copy()
preview_regions["analysis_status"] = (
    "provisional_coverage_screen"
)
preview_regions["formal_qc_status"] = "unassessed"
preview_regions.to_csv(output_regions, index=False)

print("\nSTAGE 1 — SUPERVISOR COVERAGE-SCREENING PREVIEW")
print("=" * 80)
print("Original scanned windows:", len(triage))
print("Interior assessed windows:", int(interior.sum()))
print("Boundary-unassessed windows:", int((~interior).sum()))
print("Original called windows:", int(original_calls.sum()))
print("Provisional called windows:", int(new_calls.sum()))
print("Provisional merged regions:", len(preview_regions))

print("\nALL-WINDOW COVERAGE CLASSIFICATIONS")
print(
    triage["coverage_flag"]
    .value_counts()
    .to_string()
)

print("\nPROVISIONAL DETECTOR REGIONS")
print(
    preview_regions[
        [
            "chromosome",
            "start_bp",
            "end_bp",
            "peak_w80_bins",
            "peak_score",
            "supporting_windows",
        ]
    ].to_string(index=False)
)

print("\nVERIFICATION")
print("Original measurements preserved: PASS")
print("Reference widths recalculated: PASS")
print("Filtered calls match Step 381: PASS")
print("Warnings distinguished from formal QC: PASS")
print("Mappability-masked W80: NOT DONE")
print("Independent biological validation: NOT DONE")
print("Threshold calibration: INCOMPLETE")

print("\nSaved:", output_windows)
print("Saved:", output_regions)
print("\nStep 383: PASS")
