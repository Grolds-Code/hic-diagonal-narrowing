"""Read-only audit of exploratory HCT116 narrowing calls."""

import hashlib
import importlib.util
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-detector-audit")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy==2.5.4",
        "pandas==3.1.0",
    )
    .add_local_file(
        "src/hic_diagonal_narrowing/detect.py",
        remote_path="/root/detect_reference.py",
    )
)

SOURCE = (
    "hct116_chr10_20000000_21000000_500bp_"
    "midpoint_32kb_w80_profile_v2_offdiag"
)


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def audit_detector():
    import numpy as np
    import pandas as pd

    spec = importlib.util.spec_from_file_location(
        "detect_reference",
        "/root/detect_reference.py",
    )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    profile_path = Path(
        f"/research/results/{SOURCE}.csv"
    )

    metadata_path = Path(
        f"/research/metadata/{SOURCE}.json"
    )

    metadata = json.loads(metadata_path.read_text())
    profile_bytes = profile_path.read_bytes()

    digest = hashlib.sha256(
        profile_bytes
    ).hexdigest()

    if digest != metadata["output_sha256"]:
        raise RuntimeError("Profile checksum mismatch")

    profile = pd.read_csv(profile_path)

    if len(profile) != 122:
        raise RuntimeError("Expected 122 windows")

    if not np.all(
        np.isfinite(profile["w80_bins"])
    ):
        raise RuntimeError("Invalid W80 values")

    if not profile["qc_status"].eq("unassessed").all():
        raise RuntimeError("Unexpected QC status")

    thresholds = [0.0, 0.1, 0.2, 0.3, 0.5]

    summaries = []
    baseline_scored = None

    for threshold in thresholds:
        scored, candidates = module.detect_narrowing(
            profile,
            score_threshold=threshold,
            min_reference_windows=3,
        )

        if baseline_scored is None:
            baseline_scored = scored

        summaries.append({
            "threshold": threshold,
            "called_windows": int(
                scored["candidate_call"].sum()
            ),
            "merged_candidates": len(candidates),
            "largest_peak_score": (
                float(candidates["peak_score"].max())
                if len(candidates) else None
            ),
        })

    valid = baseline_scored[
        np.isfinite(
            baseline_scored["narrowing_score"]
        )
    ]

    if len(valid) != 122:
        raise RuntimeError(
            "Not all windows received reference scores"
        )

    top = valid.sort_values(
        "narrowing_score",
        ascending=False,
    ).head(10)

    top_windows = []

    for _, row in top.iterrows():
        top_windows.append({
            "start_bp": int(row["start_bp"]),
            "end_bp": int(row["end_bp"]),
            "w80_bins": float(row["w80_bins"]),
            "reference_w80_bins": float(
                row["reference_w80_bins"]
            ),
            "score": float(row["narrowing_score"]),
            "contact_signal": float(
                row["total_contacts"]
            ),
            "pair_occupancy": float(
                row["pair_occupancy"]
            ),
            "reference_windows": int(
                row["reference_windows"]
            ),
        })

    return {
        "checksum": "PASS",
        "profile_windows": len(profile),
        "scored_windows": len(valid),
        "score_min": float(
            valid["narrowing_score"].min()
        ),
        "score_median": float(
            valid["narrowing_score"].median()
        ),
        "score_max": float(
            valid["narrowing_score"].max()
        ),
        "threshold_results": summaries,
        "top_windows": top_windows,
        "qc_assessment": "NOT_PERFORMED",
        "significance_testing": "NOT_PERFORMED",
    }


@app.local_entrypoint()
def main():
    result = audit_detector.remote()

    print("\nHCT116 EXPLORATORY DETECTOR AUDIT")
    print("=" * 65)

    print("Checksum:", result["checksum"])
    print("Profile windows:", result["profile_windows"])
    print("Scored windows:", result["scored_windows"])

    print("\nNARROWING SCORE DISTRIBUTION")
    print("Minimum:", result["score_min"])
    print("Median:", result["score_median"])
    print("Maximum:", result["score_max"])

    print("\nTHRESHOLD SENSITIVITY")
    print(
        f"{'Threshold':>12}"
        f"{'Called windows':>17}"
        f"{'Merged regions':>17}"
    )

    for row in result["threshold_results"]:
        print(
            f"{row['threshold']:>12.2f}"
            f"{row['called_windows']:>17}"
            f"{row['merged_candidates']:>17}"
        )

    print("\nTEN HIGHEST-SCORING WINDOWS")

    for row in result["top_windows"]:
        print(
            f"chr10:{row['start_bp']}-{row['end_bp']}"
            f" | W80={row['w80_bins']:.0f}"
            f" | ref={row['reference_w80_bins']:.1f}"
            f" | score={row['score']:.4f}"
            f" | signal={row['contact_signal']:.0f}"
            f" | occupancy={row['pair_occupancy']:.3f}"
        )

    print("\nQC:", result["qc_assessment"])
    print(
        "Significance:",
        result["significance_testing"],
    )

    print("\nRead-only detector audit: PASS")
