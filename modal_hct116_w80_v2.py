"""Exploratory W80 profile of the real HCT116 1 Mb contact band."""

import csv
import hashlib
import io
import json
import os
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-w80-v2")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
OUTPUT = SOURCE + "_w80_profile_v2_offdiag"

START_BP = 20_000_000
END_BP = 21_000_000
RESOLUTION = 500
WINDOW_BINS = 64
STEP_BINS = 16
MAX_SEPARATION_BINS = 64


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def compute_w80_profile():
    import numpy as np

    source_path = Path(f"/research/matrices/{SOURCE}.npz")
    source_meta_path = Path(f"/research/metadata/{SOURCE}.json")

    output_path = Path(f"/research/results/{OUTPUT}.csv")
    output_meta_path = Path(f"/research/metadata/{OUTPUT}.json")

    if output_path.exists() or output_meta_path.exists():
        raise RuntimeError(
            "Profile output already exists; refusing overwrite"
        )

    source_meta = json.loads(source_meta_path.read_text())

    source_hash = hashlib.sha256(
        source_path.read_bytes()
    ).hexdigest()

    if source_hash != source_meta["sha256"]:
        raise RuntimeError("Source checksum mismatch")

    with np.load(source_path, allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    expected_records = source_meta["selected_records"]

    if not all(
        array.shape == (expected_records,)
        for array in (x, y, counts)
    ):
        raise RuntimeError("Unexpected source dimensions")

    if not np.all(np.isfinite(counts)) or np.any(counts < 0):
        raise RuntimeError("Invalid contact counts")

    if not np.all(
        (x % RESOLUTION == 0)
        & (y % RESOLUTION == 0)
        & (x <= y)
        & (y - x <= MAX_SEPARATION_BINS * RESOLUTION)
        & (x + y >= 2 * START_BP)
        & (x + y < 2 * END_BP)
    ):
        raise RuntimeError("Invalid contact geometry")

    if not np.isclose(
        counts.sum(),
        source_meta["selected_contact_signal"],
    ):
        raise RuntimeError("Source signal mismatch")

    midpoint_twice = x + y
    separations = (y - x) // RESOLUTION

    rows = []

    total_bins = (END_BP - START_BP) // RESOLUTION

    for start_bin in range(
        0,
        total_bins - WINDOW_BINS + 1,
        STEP_BINS,
    ):
        window_start = START_BP + start_bin * RESOLUTION
        window_end = window_start + WINDOW_BINS * RESOLUTION

        selected = (
            (midpoint_twice >= 2 * window_start)
            & (midpoint_twice < 2 * window_end)
            & (separations >= 1)
        )

        local_counts = counts[selected]
        local_separations = separations[selected]

        distance_signal = np.bincount(
            local_separations,
            weights=local_counts,
            minlength=MAX_SEPARATION_BINS + 1,
        )[1:]

        if len(distance_signal) != MAX_SEPARATION_BINS:
            raise RuntimeError("Unexpected distance-band count")

        # With midpoint selection and two-sided buffering,
        # each separation has exactly WINDOW_BINS possible
        # genomic pairs in a full, unmasked window.
        mean_intensity = distance_signal / WINDOW_BINS
        total_intensity = float(mean_intensity.sum())

        positive_pairs = int(np.count_nonzero(local_counts > 0))
        possible_pairs = WINDOW_BINS * MAX_SEPARATION_BINS

        support_bands = int(
            np.count_nonzero(distance_signal > 0)
        )

        if total_intensity > 0:
            cumulative = np.cumsum(mean_intensity)

            w80_bins = int(
                np.searchsorted(
                    cumulative,
                    0.8 * total_intensity,
                    side="left",
                ) + 1
            )

            w80_bp = w80_bins * RESOLUTION
        else:
            w80_bins = ""
            w80_bp = ""

        rows.append({
            "chromosome": "chr10",
            "start_bp": window_start,
            "end_bp": window_end,
            "w80_bins": w80_bins,
            "w80_bp": w80_bp,
            "total_contacts": float(local_counts.sum()),
            "observed_pairs": positive_pairs,
            "possible_pairs": possible_pairs,
            "pair_occupancy": positive_pairs / possible_pairs,
            "supported_distance_bands": support_bands,
            "distance_support": (
                support_bands / MAX_SEPARATION_BINS
            ),
            "qc_status": "unassessed",
        })

    if len(rows) != 122:
        raise RuntimeError(
            f"Expected 122 windows, obtained {len(rows)}"
        )

    fieldnames = list(rows[0])

    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

    csv_bytes = buffer.getvalue().encode("utf-8")
    csv_hash = hashlib.sha256(csv_bytes).hexdigest()

    widths = [
        row["w80_bins"]
        for row in rows
        if isinstance(row["w80_bins"], int)
    ]

    metadata = {
        "source_dataset": SOURCE,
        "source_sha256": source_hash,
        "output_sha256": csv_hash,
        "chromosome": "chr10",
        "start_bp": START_BP,
        "end_bp_exclusive": END_BP,
        "resolution_bp": RESOLUTION,
        "window_bins": WINDOW_BINS,
        "step_bins": STEP_BINS,
        "max_separation_bins": MAX_SEPARATION_BINS,
        "w80_definition": (
            "Smallest positive separation bin (1-64) "
            "accumulating >=80% of off-diagonal "
            "per-separation mean observed contact intensity"
        ),
        "normalization": "NONE",
        "analysis_version": "v2_offdiag",
        "diagonal_policy": "excluded",
        "reference_equivalence": "122_of_122_verified",
        "missing_bin_status": "not_assessed",
        "qc_status": "unassessed",
        "interpretation": "exploratory_not_significance_tested",
        "window_count": len(rows),
        "valid_w80_windows": len(widths),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_meta_path.parent.mkdir(parents=True, exist_ok=True)

    staging_csv = output_path.with_suffix(".partial")
    staging_meta = output_meta_path.with_suffix(".partial")

    staging_csv.write_bytes(csv_bytes)
    staging_meta.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    os.replace(staging_csv, output_path)
    os.replace(staging_meta, output_meta_path)

    volume.commit()

    lowest = sorted(
        (r for r in rows if isinstance(r["w80_bins"], int)),
        key=lambda r: (r["w80_bins"], r["start_bp"]),
    )[:5]

    return {
        "result_path": str(output_path),
        "metadata_path": str(output_meta_path),
        "source_checksum": "PASS",
        "window_count": len(rows),
        "valid_w80_windows": len(widths),
        "minimum_w80_bins": min(widths) if widths else None,
        "median_w80_bins": (
            float(np.median(widths)) if widths else None
        ),
        "maximum_w80_bins": max(widths) if widths else None,
        "median_pair_occupancy": float(
            np.median([r["pair_occupancy"] for r in rows])
        ),
        "median_distance_support": float(
            np.median([r["distance_support"] for r in rows])
        ),
        "lowest_w80_windows": lowest,
    }


@app.local_entrypoint()
def main():
    result = compute_w80_profile.remote()

    print("\nHCT116 REAL-DATA W80 PROFILE")
    print("=" * 55)

    for key, value in result.items():
        if key == "lowest_w80_windows":
            continue
        print(f"{key}: {value}")

    print("\nLOWEST W80 WINDOWS (UNADJUSTED)")
    for row in result["lowest_w80_windows"]:
        print(
            f"{row['start_bp']}-{row['end_bp']}: "
            f"W80={row['w80_bins']} bins, "
            f"signal={row['total_contacts']:.1f}, "
            f"occupancy={row['pair_occupancy']:.3f}"
        )

    print("\nExploratory W80 computation: PASS")
