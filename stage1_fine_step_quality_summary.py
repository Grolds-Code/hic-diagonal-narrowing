"""Descriptive quality and peak-ambiguity audit of fine-grid calls."""

from pathlib import Path

import numpy as np
import pandas as pd

root = Path("research_artifacts/stage1")

calls = pd.read_csv(
    root / "hct116_stage1_500bp_step_calls.csv"
)
regions = pd.read_csv(
    root / "hct116_stage1_500bp_step_regions.csv"
)

if len(calls) != 93 or len(regions) != 5:
    raise RuntimeError("Unexpected fine-grid results")

if not calls["qc_status"].eq("unassessed").all():
    raise RuntimeError("Unexpected QC status")

membership = np.zeros(len(calls), dtype=int)
records = []

for i, region in enumerate(regions.itertuples(), 1):
    inside = (
        (calls["start_bp"] >= region.start_bp)
        & (calls["end_bp"] <= region.end_bp)
    )

    membership += inside.to_numpy().astype(int)

    group = calls.loc[inside].copy()

    if len(group) != region.supporting_windows:
        raise RuntimeError(
            f"Supporting-window mismatch in region {i}"
        )

    best = group["narrowing_score"].max()

    tied = group.loc[
        np.isclose(
            group["narrowing_score"].to_numpy(),
            best,
            atol=1e-12,
            rtol=0,
        )
    ].sort_values("start_bp")

    if not np.isclose(
        best, region.peak_score,
        atol=1e-12,
        rtol=0,
    ):
        raise RuntimeError(f"Peak-score mismatch: {i}")

    records.append({
        "region": i,
        "chromosome": region.chromosome,
        "start_bp": region.start_bp,
        "end_bp": region.end_bp,
        "supporting_windows": len(group),
        "peak_score": best,
        "minimum_w80_bp": int(group["w80_bp"].min()),
        "maximum_w80_bp": int(group["w80_bp"].max()),
        "peak_score_ties": len(tied),
        "earliest_tied_peak_start_bp": int(
            tied["start_bp"].min()
        ),
        "latest_tied_peak_start_bp": int(
            tied["start_bp"].max()
        ),
        "minimum_total_contacts": float(
            group["total_contacts"].min()
        ),
        "maximum_total_contacts": float(
            group["total_contacts"].max()
        ),
        "minimum_pair_occupancy": float(
            group["pair_occupancy"].min()
        ),
        "maximum_pair_occupancy": float(
            group["pair_occupancy"].max()
        ),
        "minimum_observed_distance_bands": int(
            group["supported_distance_bands"].min()
        ),
        "qc_status": "unassessed",
    })

if not np.all(membership == 1):
    raise RuntimeError(
        "Each called window must belong to exactly one region"
    )

summary = pd.DataFrame(records)

destination = (
    root / "hct116_stage1_fine_step_quality_summary.csv"
)

if destination.exists():
    raise RuntimeError(f"Refusing to overwrite: {destination}")

summary.to_csv(destination, index=False)

print("\nSTAGE 1 — FINE-GRID CANDIDATE QUALITY AUDIT")
print("=" * 85)

columns = [
    "region",
    "start_bp",
    "end_bp",
    "supporting_windows",
    "minimum_w80_bp",
    "peak_score",
    "peak_score_ties",
    "earliest_tied_peak_start_bp",
    "latest_tied_peak_start_bp",
    "minimum_total_contacts",
    "minimum_pair_occupancy",
    "minimum_observed_distance_bands",
]

print(summary[columns].to_string(index=False))

print("\nINTERPRETATION")
print("- Every called window assigned to one merged region: PASS")
print("- Coverage values are descriptive, not calibrated QC.")
print("- Observed distance bands do not establish bin validity.")
print("- Equal-score windows imply ambiguous peak localization.")
print("- Formal missing-bin QC remains unassessed.")
print("\nSaved:", destination)
print("Step 355: PASS")
