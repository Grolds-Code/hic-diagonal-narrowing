"""500-bp scan-step HCT116 W80 profile with coarse-grid validation."""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-stage1-fine-step-profile")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
START = 20_000_000
END = 21_000_000
RES = 500
WINDOW = 64
DMAX = 64

COLUMNS = [
    "chromosome", "start_bp", "end_bp",
    "w80_bins", "w80_bp", "total_contacts",
    "observed_pairs", "possible_pairs",
    "pair_occupancy", "supported_distance_bands",
    "distance_support", "qc_status",
]


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def compute_profile():
    import numpy as np

    root = Path("/research")

    source_path = root / "matrices" / f"{SOURCE}.npz"
    metadata_path = root / "metadata" / f"{SOURCE}.json"

    source_bytes = source_path.read_bytes()
    metadata = json.loads(metadata_path.read_text())
    digest = hashlib.sha256(source_bytes).hexdigest()

    if digest != metadata["sha256"]:
        raise RuntimeError("Source checksum mismatch")

    if (
        metadata["start_bp"] != START
        or metadata["end_bp_exclusive"] != END
        or metadata["resolution_bp"] != RES
        or metadata["max_separation_bp"] != RES * DMAX
    ):
        raise RuntimeError("Unexpected source geometry")

    with np.load(io.BytesIO(source_bytes), allow_pickle=False) as data:
        x = data["bin_x_bp"].astype(np.int64)
        y = data["bin_y_bp"].astype(np.int64)
        raw_counts = data["counts"]

    if (
        len(x) != 45_831
        or len(y) != len(x)
        or len(raw_counts) != len(x)
        or np.any(~np.isfinite(raw_counts))
        or np.any(raw_counts < 0)
        or np.any(raw_counts != np.rint(raw_counts))
        or np.any(x > y)
        or np.any((y - x) % RES != 0)
    ):
        raise RuntimeError("Invalid original contact records")

    counts = raw_counts.astype(np.int64)
    distances = (y - x) // RES

    # Exclude self-contacts, exactly as in the baseline profile.
    selected = (distances >= 1) & (distances <= DMAX)

    if int(selected.sum()) != 43_918:
        raise RuntimeError("Unexpected off-diagonal record count")

    x = x[selected]
    y = y[selected]
    distances = distances[selected].astype(int)
    counts = counts[selected]

    n_bins = (END - START) // RES

    # Assign each contact by its genomic midpoint.
    # Retain external endpoints when the midpoint is inside.
    midpoint_bin = (x + y - 2 * START) // (2 * RES)

    if np.any(midpoint_bin < 0) or np.any(midpoint_bin >= n_bins):
        raise RuntimeError("Contact midpoint outside target region")

    # Integer-valued contact totals and occupied pairs,
    # grouped by midpoint bin and genomic separation.
    signals = np.zeros((n_bins, DMAX), dtype=np.int64)
    occupied = np.zeros((n_bins, DMAX), dtype=np.int64)

    index = (midpoint_bin.astype(int), distances - 1)

    np.add.at(signals, index, counts)
    np.add.at(occupied, index, 1)

    signal_prefix = np.vstack([
        np.zeros((1, DMAX), dtype=np.int64),
        np.cumsum(signals, axis=0),
    ])

    occupied_prefix = np.vstack([
        np.zeros((1, DMAX), dtype=np.int64),
        np.cumsum(occupied, axis=0),
    ])

    starts = np.arange(n_bins - WINDOW + 1)
    ends = starts + WINDOW

    band_signal = signal_prefix[ends] - signal_prefix[starts]
    band_occupied = (
        occupied_prefix[ends] - occupied_prefix[starts]
    )

    totals = band_signal.sum(axis=1)
    observed = band_occupied.sum(axis=1)
    supported = np.count_nonzero(band_occupied, axis=1)

    widths = np.full(len(starts), np.nan)

    positive = totals > 0
    cumulative = np.cumsum(band_signal[positive], axis=1)

    widths[positive] = (
        np.argmax(
            cumulative >= 0.8 * totals[positive, None],
            axis=1,
        ) + 1
    )

    # Under the original unmasked, complete midpoint geometry,
    # a 64-bin window has 64 eligible pairs per separation.
    # This is structural support, not independent missing-bin QC.
    possible = WINDOW * DMAX

    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=COLUMNS,
        lineterminator="\n",
    )

    writer.writeheader()

    for i, start_bin in enumerate(starts):
        width = widths[i]

        writer.writerow({
            "chromosome": "chr10",
            "start_bp": START + int(start_bin) * RES,
            "end_bp": START + int(ends[i]) * RES,
            "w80_bins": width,
            "w80_bp": width * RES,
            "total_contacts": float(totals[i]),
            "observed_pairs": int(observed[i]),
            "possible_pairs": possible,
            "pair_occupancy": float(observed[i] / possible),
            "supported_distance_bands": int(supported[i]),
            "distance_support": float(supported[i] / DMAX),
            "qc_status": (
                "unassessed" if positive[i] else "no_signal"
            ),
        })

    return {
        "source_sha256": digest,
        "offdiagonal_records": int(selected.sum()),
        "csv_bytes": output.getvalue().encode("utf-8"),
    }


@app.local_entrypoint()
def main():
    import pandas as pd

    result = compute_profile.remote()

    expected_sha = (
        "fb01bc848b4fecb342a450b8d87165393"
        "a89de019fc27c82e72d8519c07f5de5"
    )

    if result["source_sha256"] != expected_sha:
        raise RuntimeError("Unexpected HCT116 source version")

    fine = pd.read_csv(io.BytesIO(result["csv_bytes"]))

    baseline_path = Path(
        "research_artifacts/stage1/w80_profile_v2_offdiag.csv"
    )

    baseline = pd.read_csv(baseline_path)

    if len(fine) != 1937 or len(baseline) != 122:
        raise RuntimeError("Unexpected window counts")

    # Every 16th fine-grid window must match the original
    # 8-kb-step profile, including all coverage columns.
    coarse = fine.iloc[::16].reset_index(drop=True)

    pd.testing.assert_frame_equal(
        coarse[COLUMNS],
        baseline[COLUMNS].reset_index(drop=True),
        check_dtype=False,
        check_exact=False,
        atol=1e-12,
        rtol=0,
    )

    destination = Path(
        "research_artifacts/stage1/"
        "hct116_stage1_500bp_step_profile.csv"
    )

    if destination.exists():
        raise RuntimeError(
            f"Refusing to overwrite existing file: {destination}"
        )

    destination.write_bytes(result["csv_bytes"])

    print("\nSTAGE 1 — 500-BP SCAN-STEP VERIFICATION")
    print("=" * 65)
    print("Source checksum: PASS")
    print("Off-diagonal records:", result["offdiagonal_records"])
    print("Input data resolution: 500 bp")
    print("Measurement window: 32 kb")
    print("Maximum separation: 32 kb")
    print("New scan step: 500 bp")
    print("New profile windows:", len(fine))
    print("Original windows reproduced:", len(coarse))
    print("W80 and coverage baseline match: PASS")
    print("Missing-bin QC: NOT ASSESSED")
    print("Output:", destination)

    print("\nFIRST FIVE WINDOWS")
    print(
        fine[
            [
                "start_bp", "end_bp", "w80_bins",
                "total_contacts", "pair_occupancy",
            ]
        ].head().to_string(index=False)
    )

    print("\nStep 353: PASS")
