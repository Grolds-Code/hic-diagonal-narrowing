"""Genome-coordinate scanning of local Hi-C diagonal width.

Research prototype v0.1.
Window-level scanning only; candidate calling is not yet
implemented or validated.
"""

import numpy as np
import pandas as pd

from hic_diagonal_narrowing.width import estimate_w80
from hic_diagonal_narrowing.qc import summarize_local_support


def scan_local_width(
    matrix: np.ndarray,
    chromosome: str,
    *,
    region_start_bp: int = 0,
    resolution_bp: int = 500,
    window_bins: int = 64,
    step_bins: int = 16,
    max_separation_bins: int = 64,
) -> pd.DataFrame:
    """Scan local W80 using overlapping genomic windows.

    Coordinates are zero-based, half-open genomic intervals.
    The input matrix represents consecutive bins starting at
    region_start_bp.

    No statistical significance or biological confidence
    is implied by the returned widths.
    """
    matrix = np.asarray(matrix)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square")

    if not np.all(np.isfinite(matrix)):
        raise ValueError("Matrix contains non-finite values")

    if np.any(matrix < 0):
        raise ValueError("Matrix contains negative values")

    if not np.allclose(matrix, matrix.T):
        raise ValueError("Matrix must be symmetric")

    n_bins = matrix.shape[0]

    if not chromosome:
        raise ValueError("Chromosome must be specified")

    if region_start_bp < 0 or resolution_bp <= 0:
        raise ValueError("Invalid genomic coordinates or resolution")

    if not (1 <= window_bins <= n_bins):
        raise ValueError("Invalid window size")

    if step_bins < 1:
        raise ValueError("Step size must be positive")

    # Conservative restriction for this initial implementation.
    if not (1 <= max_separation_bins <= window_bins):
        raise ValueError(
            "Maximum separation must not exceed window size"
        )

    records = []

    for start_bin in range(
        0,
        n_bins - window_bins + 1,
        step_bins,
    ):
        end_bin = start_bin + window_bins
        region = (start_bin, end_bin)

        qc = summarize_local_support(
            matrix,
            region,
            max_separation_bins=max_separation_bins,
        )

        if qc["total_contact_signal"] == 0:
            width = np.nan
            status = "no_signal"
        else:
            width = estimate_w80(
                matrix,
                region,
                max_separation_bins=max_separation_bins,
            )
            status = "unassessed"

        records.append({
            "chromosome": chromosome,
            "start_bp": region_start_bp
                        + start_bin * resolution_bp,
            "end_bp": region_start_bp
                      + end_bin * resolution_bp,
            "start_bin": start_bin,
            "end_bin": end_bin,
            "w80_bins": width,
            "w80_bp": (
                width * resolution_bp
                if np.isfinite(width)
                else np.nan
            ),
            "total_contacts": qc["total_contact_signal"],
            "pair_occupancy": qc["pair_occupancy"],
            "distance_support": qc[
                "distance_support_fraction"
            ],
            "qc_status": status,
        })

    return pd.DataFrame.from_records(records)
