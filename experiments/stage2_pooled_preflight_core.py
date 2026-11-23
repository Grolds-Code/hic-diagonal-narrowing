"""Read-only checks for the pooled HCT116 Stage 1 extraction used by Stage 2.

This reproduces the original 500-bp scan from midpoint-assigned pairs.
It does not establish coverage validity or biological replication.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np

START = 20_000_000
END = 21_000_000
RES = 500
WINDOW = 64
MAX_DISTANCE = 64
EXPECTED_SHA256 = (
    "fb01bc848b4fecb342a450b8d87165393"
    "a89de019fc27c82e72d8519c07f5de5"
)
SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"

# These are Stage 1 *footprints*, not inferred genomic boundaries.
FOOTPRINTS = {
    "R1": (20_258_500, 20_336_000),
    "R2": (20_427_000, 20_459_500),
}


def midpoint_w80_profile(
    x: np.ndarray,
    y: np.ndarray,
    counts: np.ndarray,
    *,
    start: int,
    end: int,
    resolution: int,
    window_bins: int,
    max_distance_bins: int,
):
    """Compute original unmasked W80 across all one-bin-step windows.

    Midpoint selection and distance aggregation follow the Stage 1 source
    estimator; no per-bin coverage or mappability masking is applied.
    """
    x, y, counts = (
        np.asarray(x), np.asarray(y), np.asarray(counts)
    )
    if x.ndim != 1 or y.shape != x.shape or counts.shape != x.shape:
        raise ValueError("Contact arrays must be 1D and equal-length")
    if (
        resolution <= 0 or window_bins <= 0 or max_distance_bins <= 0
        or start >= end or (end - start) % resolution
    ):
        raise ValueError("Invalid profile geometry")
    if (
        not np.all(np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
        or np.any(x > y)
        or np.any((x - start) % resolution)
        or np.any((y - start) % resolution)
    ):
        raise ValueError("Invalid contact records")

    n_bins = (end - start) // resolution
    if window_bins > n_bins:
        raise ValueError("Measurement window exceeds target region")

    sep = (y - x) // resolution
    selected = (sep >= 1) & (sep <= max_distance_bins)
    x, y, sep = x[selected], y[selected], sep[selected]
    weights = counts[selected].astype(np.int64)

    midpoints = (x + y - 2 * start) // (2 * resolution)
    if np.any((midpoints < 0) | (midpoints >= n_bins)):
        raise ValueError("Off-diagonal record midpoint outside target")
    signal = np.zeros((n_bins, max_distance_bins), dtype=np.int64)
    np.add.at(signal, (midpoints, sep.astype(int) - 1), weights)
    prefix = np.vstack((
        np.zeros((1, max_distance_bins), dtype=np.int64),
        np.cumsum(signal, axis=0),
    ))

    starts = np.arange(n_bins - window_bins + 1, dtype=np.int64)
    bands = prefix[starts + window_bins] - prefix[starts]
    totals = bands.sum(axis=1)
    widths = np.full(len(starts), -1, dtype=np.int64)
    good = totals > 0
    if np.any(good):
        widths[good] = 1 + np.argmax(
            np.cumsum(bands[good], axis=1)
            >= 0.8 * totals[good, None],
            axis=1,
        )
    return widths, totals, int(selected.sum())


def audit_pooled_source(
    source_path: Path,
    metadata_path: Path,
    triage_csv_path: Path,
):
    """Reproduce every Stage 1 W80 and contact-total row; fail closed."""
    source_bytes = source_path.read_bytes()
    digest = hashlib.sha256(source_bytes).hexdigest()
    if digest != EXPECTED_SHA256:
        raise RuntimeError("Pooled contact NPZ SHA256 mismatch")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if (
        metadata.get("sha256") != EXPECTED_SHA256
        or metadata.get("start_bp") != START
        or metadata.get("end_bp_exclusive") != END
        or metadata.get("resolution_bp") != RES
        or metadata.get("max_separation_bp") != MAX_DISTANCE * RES
    ):
        raise RuntimeError("Pooled contact provenance metadata mismatch")

    with np.load(io.BytesIO(source_bytes), allow_pickle=False) as z:
        x = z["bin_x_bp"].astype(np.int64)
        y = z["bin_y_bp"].astype(np.int64)
        counts = z["counts"]

    if len(x) != 45_831:
        raise RuntimeError("Pooled source record count changed")

    widths, totals, offdiag_count = midpoint_w80_profile(
        x, y, counts,
        start=START, end=END, resolution=RES,
        window_bins=WINDOW, max_distance_bins=MAX_DISTANCE,
    )
    if len(widths) != 1937 or offdiag_count != 43_918:
        raise RuntimeError("Pooled source geometry or contact count changed")

    with triage_csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1937:
        raise RuntimeError("Stage 1 QC triage scan length changed")

    original_calls = []
    screened_calls = []
    for i, row in enumerate(rows):
        expected_start = START + i * RES
        if (
            row.get("chromosome") != "chr10"
            or int(row["start_bp"]) != expected_start
            or int(row["end_bp"]) != expected_start + WINDOW * RES
        ):
            raise RuntimeError(f"Stage 1 coordinate mismatch at row {i}")
        if (
            int(float(row["w80_bins"])) != int(widths[i])
            or int(float(row["total_contacts"])) != int(totals[i])
        ):
            raise RuntimeError(f"W80/contact-total mismatch at {expected_start}")
        if row["original_candidate_call"].strip().lower() == "true":
            original_calls.append(row)
        if row["preview_candidate_call"].strip().lower() == "true":
            if (
                row["coverage_flag"] != "no_coverage_warning"
                or row["provisionally_eligible"].strip().lower() != "true"
            ):
                raise RuntimeError("Provisional screen contains flagged window")
            screened_calls.append(row)

    if len(original_calls) != 93 or len(screened_calls) != 25:
        raise RuntimeError("Original or provisional Stage 1 calls changed")

    candidates = {}
    for name, (foot_start, foot_end) in FOOTPRINTS.items():
        subset = [
            row for row in screened_calls
            if (
                int(row["start_bp"]) >= foot_start
                and int(row["end_bp"]) <= foot_end
            )
        ]
        candidates[name] = {
            "footprint_start_bp": foot_start,
            "footprint_end_bp": foot_end,
            "provisional_called_windows": len(subset),
            "w80_bp_observed": sorted({
                int(float(row["w80_bp"])) for row in subset
            }),
        }
    if (
        candidates["R1"]["provisional_called_windows"] != 23
        or candidates["R2"]["provisional_called_windows"] != 2
    ):
        raise RuntimeError("Provisional candidate footprint accounting changed")

    return {
        "source_sha256": digest,
        "source_records": len(x),
        "offdiagonal_records": offdiag_count,
        "scan_windows": len(widths),
        "stage1_rows_reproduced": len(rows),
        "original_stage1_called_windows": len(original_calls),
        "provisional_coverage_screen_calls": len(screened_calls),
        "candidates": candidates,
        "formal_qc_status": "unassessed",
        "interpretation": (
            "Exact pooled-data regression only; overlapping windows are "
            "not independent; no biological replicates tested."
        ),
    }
