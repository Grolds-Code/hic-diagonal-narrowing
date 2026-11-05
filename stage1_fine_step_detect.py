"""Compare original and 500-bp-step HCT116 narrowing calls."""

from pathlib import Path

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.detect import detect_narrowing

ROOT = Path("research_artifacts/stage1")

BASELINE = ROOT / "w80_profile_v2_offdiag.csv"
FINE = ROOT / "hct116_stage1_500bp_step_profile.csv"

baseline = pd.read_csv(BASELINE)
fine = pd.read_csv(FINE)

if len(baseline) != 122 or len(fine) != 1937:
    raise RuntimeError("Unexpected input profile lengths")

if not fine["start_bp"].diff().iloc[1:].eq(500).all():
    raise RuntimeError("Fine scan does not advance by 500 bp")

if not (fine["end_bp"] - fine["start_bp"]).eq(32000).all():
    raise RuntimeError("Unexpected detection window size")

if not fine["qc_status"].eq("unassessed").all():
    raise RuntimeError("Unexpected fine-profile QC status")

pd.testing.assert_frame_equal(
    fine.iloc[::16].reset_index(drop=True),
    baseline.reset_index(drop=True),
    check_dtype=False,
    atol=1e-12,
    rtol=0,
)

kwargs = {
    "score_threshold": 0.2,
    "min_reference_windows": 3,
}

old_scored, old_regions = detect_narrowing(
    baseline, **kwargs
)

new_scored, new_regions = detect_narrowing(
    fine, **kwargs
)

old_calls = old_scored.loc[
    old_scored["candidate_call"]
].copy()

new_calls = new_scored.loc[
    new_scored["candidate_call"]
].copy()

if len(old_calls) != 5 or len(old_regions) != 2:
    raise RuntimeError("Original detector calls changed")

if len(new_scored) != 1937:
    raise RuntimeError("Fine detector output incomplete")

if not set(new_calls.index).issubset(new_scored.index):
    raise RuntimeError("Invalid candidate indexing")

outputs = {
    ROOT / "hct116_stage1_500bp_step_scored.csv":
        new_scored,
    ROOT / "hct116_stage1_500bp_step_calls.csv":
        new_calls,
    ROOT / "hct116_stage1_500bp_step_regions.csv":
        new_regions,
}

for path in outputs:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite: {path}")

print("\nSTAGE 1: COARSE VERSUS FINE DETECTION")
print("=" * 75)
print(f"Original windows:         {len(baseline)}")
print(f"Original candidate calls: {len(old_calls)}")
print(f"Original merged regions:  {len(old_regions)}")
print(f"Fine-grid windows:        {len(fine)}")
print(f"Fine-grid candidate calls:{len(new_calls):>5}")
print(f"Fine-grid merged regions: {len(new_regions)}")

print("\nORIGINAL FIVE CALLS UNDER BOTH REFERENCE SETS")
print("=" * 75)

old_subset = old_calls[
    ["start_bp", "end_bp", "w80_bins", "narrowing_score"]
].rename(columns={
    "narrowing_score": "original_score"
})

fine_subset = new_scored[
    ["start_bp", "end_bp", "narrowing_score", "candidate_call"]
].rename(columns={
    "narrowing_score": "fine_reference_score",
    "candidate_call": "fine_reference_call",
})

comparison = old_subset.merge(
    fine_subset,
    on=["start_bp", "end_bp"],
    validate="one_to_one",
)

if len(comparison) != 5:
    raise RuntimeError("Failed to compare all original calls")

print(comparison.to_string(index=False))

print("\nFINE-GRID RESULTS AROUND ORIGINAL CANDIDATES")
print("=" * 75)

for label, region in zip(("A", "B"), old_regions.itertuples()):
    overlapping = new_calls.loc[
        (new_calls["start_bp"] < region.end_bp)
        & (new_calls["end_bp"] > region.start_bp)
    ]

    print(f"\nRegion {label}")
    print(
        f"Original footprint: chr10:"
        f"{region.start_bp:,}-{region.end_bp:,}"
    )
    print("Overlapping fine-grid calls:", len(overlapping))

    if not overlapping.empty:
        peak = overlapping.sort_values(
            "narrowing_score", ascending=False
        ).iloc[0]

        print(
            "Highest-scoring overlapping window:",
            f"chr10:{int(peak.start_bp):,}-"
            f"{int(peak.end_bp):,}"
        )
        print("W80 (bp):", int(peak.w80_bp))
        print(
            "Narrowing score:",
            round(float(peak.narrowing_score), 6)
        )

print("\nFINE-GRID MERGED REGIONS")
print("=" * 75)
if new_regions.empty:
    print("No regions called at the specified threshold.")
else:
    print(new_regions.to_string(index=False, max_rows=25))

for path, data in outputs.items():
    data.to_csv(path, index=False)
    print("\nSaved:", path)

print("\nInterpretation:")
print("- Input resolution: 500 bp")
print("- Measurement window: 32 kb")
print("- Fine scan step: 500 bp")
print("- Overlapping windows are not independent")
print("- QC remains unassessed")
print("- Coordinates remain exploratory")
print("\nStep 354: PASS")
