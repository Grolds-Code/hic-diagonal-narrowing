"""Exploratory Hi-C diagonal narrowing candidate detection.

Uses non-overlapping reference windows and merges overlapping
candidate windows. Scores are descriptive, not p-values.
"""

import numpy as np
import pandas as pd


def detect_narrowing(
    profile: pd.DataFrame,
    *,
    score_threshold: float = 0.5,
    min_reference_windows: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return scored windows and merged candidate intervals.

    The reference for each window is the median W80 of
    non-overlapping, informative windows on the same chromosome.

    Overlapping candidate windows are merged. Reported spans
    are detection footprints, not precise biological boundaries.

    QC acceptance thresholds are not yet implemented.
    """

    required = {
        "chromosome", "start_bp", "end_bp",
        "w80_bins", "total_contacts", "qc_status",
    }

    missing = required - set(profile.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    if not np.isfinite(score_threshold) or score_threshold < 0:
        raise ValueError("Invalid score threshold")

    if min_reference_windows < 1:
        raise ValueError("At least one reference window required")

    scored = profile.copy().sort_values(
        ["chromosome", "start_bp", "end_bp"]
    ).reset_index(drop=True)

    if scored.empty:
        raise ValueError("Empty width profile")

    if (
        (scored["start_bp"] < 0).any()
        or (scored["end_bp"] <= scored["start_bp"]).any()
    ):
        raise ValueError("Invalid genomic intervals")

    scored["reference_w80_bins"] = np.nan
    scored["reference_windows"] = 0
    scored["narrowing_score"] = np.nan
    scored["candidate_call"] = False

    for chromosome, group in scored.groupby(
        "chromosome", sort=False
    ):
        indices = group.index.to_numpy()

        starts = group["start_bp"].to_numpy()
        ends = group["end_bp"].to_numpy()
        widths = group["w80_bins"].to_numpy(dtype=float)
        contacts = group["total_contacts"].to_numpy(dtype=float)

        eligible = (
            np.isfinite(widths)
            & (widths > 0)
            & np.isfinite(contacts)
            & (contacts > 0)
            & group["qc_status"].eq("unassessed").to_numpy()
        )

        for local_idx, global_idx in enumerate(indices):
            if not eligible[local_idx]:
                continue

            # Reference windows must not overlap this window.
            disjoint = (
                (ends <= starts[local_idx])
                | (starts >= ends[local_idx])
            )

            reference_mask = eligible & disjoint
            n_reference = int(reference_mask.sum())

            scored.at[global_idx, "reference_windows"] = n_reference

            if n_reference < min_reference_windows:
                continue

            reference_width = float(
                np.median(widths[reference_mask])
            )

            score = float(
                np.log2(reference_width / widths[local_idx])
            )

            scored.at[
                global_idx, "reference_w80_bins"
            ] = reference_width

            scored.at[
                global_idx, "narrowing_score"
            ] = score

            scored.at[
                global_idx, "candidate_call"
            ] = score > score_threshold

    called = scored.loc[scored["candidate_call"]]

    candidates = []
    current = None

    for _, row in called.iterrows():
        chromosome = row["chromosome"]
        start = int(row["start_bp"])
        end = int(row["end_bp"])
        score = float(row["narrowing_score"])

        if (
            current is not None
            and current["chromosome"] == chromosome
            and start <= current["end_bp"]
        ):
            current["end_bp"] = max(
                current["end_bp"], end
            )
            current["supporting_windows"] += 1

            if score > current["peak_score"]:
                current["peak_score"] = score
                current["peak_w80_bins"] = float(
                    row["w80_bins"]
                )
                current["peak_start_bp"] = start
                current["peak_end_bp"] = end

        else:
            current = {
                "chromosome": chromosome,
                "start_bp": start,
                "end_bp": end,
                "peak_start_bp": start,
                "peak_end_bp": end,
                "peak_w80_bins": float(row["w80_bins"]),
                "peak_score": score,
                "supporting_windows": 1,
                "qc_status": "unassessed",
            }
            candidates.append(current)

    return scored, pd.DataFrame(
        candidates,
        columns=[
            "chromosome", "start_bp", "end_bp",
            "peak_start_bp", "peak_end_bp",
            "peak_w80_bins", "peak_score",
            "supporting_windows", "qc_status",
        ],
    )
