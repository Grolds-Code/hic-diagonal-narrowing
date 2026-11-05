"""Integrate verified exploratory HCT116 Stage 1 calls and QC audits."""

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
quality = pd.read_csv(
    root / "hct116_stage1_fine_step_quality_summary.csv"
)
support = pd.read_csv(
    root / "hct116_stage1_endpoint_support_coordinates.csv"
).fillna({"zero_contact_bin_starts_bp": ""})
sensitivity = pd.read_csv(
    root / "hct116_stage1_zero_bin_sensitivity.csv"
)

if (
    len(calls) != 93
    or len(regions) != 5
    or len(quality) != 5
    or len(sensitivity) != 5
):
    raise RuntimeError("Unexpected Stage 1 source dimensions")

candidate_support = support.loc[
    support["kind"] == "candidate"
]

if len(candidate_support) != 5:
    raise RuntimeError("Expected five candidate support records")

records = []

for i, region in enumerate(regions.itertuples(), 1):
    label = f"R{i}"

    q = quality.loc[quality["region"] == i]
    s = candidate_support.loc[
        candidate_support["label"] == label
    ]
    m = sensitivity.loc[
        sensitivity["region"] == label
    ]

    if any(len(frame) != 1 for frame in (q, s, m)):
        raise RuntimeError(f"Missing or duplicate audit for {label}")

    q = q.iloc[0]
    s = s.iloc[0]
    m = m.iloc[0]

    group = calls.loc[
        (calls["start_bp"] >= region.start_bp)
        & (calls["end_bp"] <= region.end_bp)
    ]

    if len(group) != int(region.supporting_windows):
        raise RuntimeError(f"Window membership mismatch: {label}")

    peak_score = float(group["narrowing_score"].max())

    ties = group.loc[
        np.isclose(
            group["narrowing_score"],
            peak_score,
            rtol=0,
            atol=1e-12,
        )
    ].sort_values("start_bp")

    representative = ties.iloc[len(ties) // 2]

    if (
        len(ties) != int(q["peak_score_ties"])
        or not np.isclose(
            peak_score, region.peak_score,
            rtol=0, atol=1e-12,
        )
        or int(representative["start_bp"]) != int(s["start_bp"])
        or int(representative["start_bp"]) != int(m["start_bp"])
        or int(representative["end_bp"]) != int(m["end_bp"])
        or int(representative["w80_bp"]) != int(m["original_w80_bp"])
        or int(representative["total_contacts"]) != int(
            m["original_total_contacts"]
        )
        or int(representative["observed_pairs"]) != int(
            m["original_observed_pairs"]
        )
        or int(q["start_bp"]) != int(region.start_bp)
        or int(q["end_bp"]) != int(region.end_bp)
    ):
        raise RuntimeError(f"Cross-artifact mismatch: {label}")

    zero_text = str(s["zero_contact_bin_starts_bp"])
    zero_positions = (
        [int(v) for v in zero_text.split(";")]
        if zero_text else []
    )

    if len(zero_positions) != int(s["zero_contact_bins"]):
        raise RuntimeError(f"Zero-contact-bin mismatch: {label}")

    records.append({
        "candidate_id": label,
        "chromosome": region.chromosome,
        "footprint_start_bp": int(region.start_bp),
        "footprint_end_bp": int(region.end_bp),
        "coordinate_convention": "0-based_half-open",
        "input_resolution_bp": 500,
        "measurement_window_bp": 32000,
        "scan_step_bp": 500,
        "supporting_windows": len(group),
        "peak_score": peak_score,
        "peak_w80_bp": int(
            float(region.peak_w80_bins) * 500
        ),
        "equal_peak_score_windows": len(ties),
        "earliest_tied_start_bp": int(
            ties["start_bp"].min()
        ),
        "latest_tied_start_bp": int(
            ties["start_bp"].max()
        ),
        "representative_start_bp": int(
            representative["start_bp"]
        ),
        "representative_end_bp": int(
            representative["end_bp"]
        ),
        "representative_total_contacts": int(
            representative["total_contacts"]
        ),
        "representative_pair_occupancy": float(
            representative["pair_occupancy"]
        ),
        "minimum_called_window_occupancy": float(
            q["minimum_pair_occupancy"]
        ),
        "minimum_observed_distance_bands": int(
            q["minimum_observed_distance_bands"]
        ),
        "representative_zero_contact_bins": len(zero_positions),
        "longest_zero_contact_run_bins": int(
            s["longest_zero_run_bins"]
        ),
        "zero_contact_bin_starts_bp": zero_text,
        "zero_bin_exclusion_w80_bp": int(
            m["masked_w80_bp"]
        ),
        "zero_bin_exclusion_w80_change_bp": int(
            m["w80_change_bp"]
        ),
        "bin_support_warning": (
            "zero_contact_bins_observed"
            if zero_positions
            else "none_observed_in_representative_window"
        ),
        "missing_bin_status": "not_assessed",
        "formal_qc_status": "unassessed",
        "interpretation": "exploratory_candidate",
    })

table = pd.DataFrame(records)

if (
    len(table) != 5
    or int(table["supporting_windows"].sum()) != 93
    or not table["formal_qc_status"].eq("unassessed").all()
    or not table["zero_bin_exclusion_w80_change_bp"].eq(0).all()
):
    raise RuntimeError("Integrated table failed validation")

output = root / "hct116_stage1_candidate_deliverable_v1.csv"

if output.exists():
    raise RuntimeError(f"Refusing to overwrite: {output}")

table.to_csv(output, index=False)

print("\nSTAGE 1 — INTEGRATED CANDIDATE DELIVERABLE")
print("=" * 105)

display_columns = [
    "candidate_id",
    "footprint_start_bp",
    "footprint_end_bp",
    "supporting_windows",
    "peak_w80_bp",
    "peak_score",
    "equal_peak_score_windows",
    "representative_total_contacts",
    "representative_pair_occupancy",
    "representative_zero_contact_bins",
    "zero_bin_exclusion_w80_change_bp",
]

print(
    table[display_columns].to_string(
        index=False,
        float_format=lambda v: f"{v:.4f}",
    )
)

print("\nVERIFICATION")
print("Five merged genomic regions: PASS")
print("Ninety-three supporting windows: PASS")
print("Peak scores and coverage reconciled: PASS")
print("Zero-contact-bin audits reconciled: PASS")
print("Sensitivity measurements reconciled: PASS")
print("Formal missing-bin QC: NOT ASSESSED")
print("Biological significance: NOT ESTABLISHED")
print("Saved:", output)
print("\nStep 362: PASS")
