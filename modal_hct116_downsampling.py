"""Controlled contact-depth sensitivity for real HCT116 W80."""

import csv
import hashlib
import io
import json
import os
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-downsampling")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"
OUTPUT = SOURCE + "_w80_downsampling_v1"

RES = 500
MAX_SEP = 64
REPLICATES = 100
SEED = 20261023
RETENTION = (0.50, 0.25, 0.10, 0.05, 0.01)


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=2048,
    timeout=600,
)
def test_downsampling():
    import numpy as np

    source_path = Path(f"/research/matrices/{SOURCE}.npz")
    source_meta_path = Path(f"/research/metadata/{SOURCE}.json")

    profile_path = Path(f"/research/results/{PROFILE}.csv")
    profile_meta_path = Path(f"/research/metadata/{PROFILE}.json")

    result_path = Path(f"/research/results/{OUTPUT}.csv")
    result_meta_path = Path(f"/research/metadata/{OUTPUT}.json")

    if result_path.exists() or result_meta_path.exists():
        raise RuntimeError(
            "Output already exists; refusing to overwrite"
        )

    source_meta = json.loads(source_meta_path.read_text())
    profile_meta = json.loads(profile_meta_path.read_text())

    source_bytes = source_path.read_bytes()
    profile_bytes = profile_path.read_bytes()

    source_hash = hashlib.sha256(source_bytes).hexdigest()
    profile_hash = hashlib.sha256(profile_bytes).hexdigest()

    if source_hash != source_meta["sha256"]:
        raise RuntimeError("Source checksum mismatch")

    if profile_hash != profile_meta["output_sha256"]:
        raise RuntimeError("W80 profile checksum mismatch")

    with np.load(source_path, allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    rows = list(csv.DictReader(
        io.StringIO(profile_bytes.decode("utf-8"))
    ))

    if len(rows) != 122:
        raise RuntimeError("Unexpected window count")

    if (
        len(counts) != source_meta["selected_records"]
        or x.shape != counts.shape
        or y.shape != counts.shape
    ):
        raise RuntimeError("Source dimensions mismatch")

    if not np.all(np.isfinite(counts)):
        raise RuntimeError("Nonfinite contact counts")

    if np.any(counts < 0):
        raise RuntimeError("Negative contact counts")

    # Binomial thinning requires actual integer contact counts.
    # Do not silently round a non-integer normalized matrix.
    if not np.all(counts == np.rint(counts)):
        raise RuntimeError(
            "Noninteger counts: binomial thinning is invalid"
        )

    integer_counts = counts.astype(np.int64)

    if not np.array_equal(
        integer_counts.astype(float), counts
    ):
        raise RuntimeError("Contact-count conversion mismatch")

    midpoint_twice = x + y
    separation = (y - x) // RES

    windows = []
    distances = []
    baseline = []
    baseline_signal = []

    for row in rows:
        start = int(row["start_bp"])
        end = int(row["end_bp"])

        if end - start != 64 * RES:
            raise RuntimeError("Unexpected window width")

        mask = (
            (midpoint_twice >= 2 * start)
            & (midpoint_twice < 2 * end)
            & (separation >= 1)
            & (separation <= MAX_SEP)
        )

        indices = np.flatnonzero(mask)
        local_sep = separation[indices]

        reference_signal = np.bincount(
            local_sep,
            weights=integer_counts[indices],
            minlength=MAX_SEP + 1,
        )[1:]

        total = reference_signal.sum()

        if total <= 0:
            raise RuntimeError("No baseline window signal")

        reference_w80 = int(
            np.searchsorted(
                np.cumsum(reference_signal),
                0.8 * total,
            ) + 1
        )

        if reference_w80 != int(row["w80_bins"]):
            raise RuntimeError(
                f"Baseline W80 mismatch at {start}"
            )

        if not np.isclose(
            total, float(row["total_contacts"])
        ):
            raise RuntimeError(
                f"Baseline signal mismatch at {start}"
            )

        windows.append(indices)
        distances.append(local_sep)
        baseline.append(reference_w80)
        baseline_signal.append(float(total))

    baseline = np.asarray(baseline, dtype=int)

    rng = np.random.default_rng(SEED)

    detailed_rows = []
    summaries = []

    for retention in RETENTION:
        simulated = np.full(
            (REPLICATES, len(rows)),
            np.nan,
            dtype=float,
        )

        for replicate in range(REPLICATES):
            # One draw across the entire genomic region retains
            # shared randomness between overlapping windows.
            thinned = rng.binomial(
                integer_counts,
                retention,
            )

            for index, (positions, local_sep) in enumerate(
                zip(windows, distances)
            ):
                distance_signal = np.bincount(
                    local_sep,
                    weights=thinned[positions],
                    minlength=MAX_SEP + 1,
                )[1:]

                total = float(distance_signal.sum())

                if total > 0:
                    simulated[replicate, index] = int(
                        np.searchsorted(
                            np.cumsum(distance_signal),
                            0.8 * total,
                        ) + 1
                    )

        delta = simulated - baseline[np.newaxis, :]
        valid = np.isfinite(delta)
        valid_delta = delta[valid]

        if len(valid_delta) == 0:
            raise RuntimeError(
                "No valid downsampling estimates"
            )

        summaries.append({
            "retention": retention,
            "replicates": REPLICATES,
            "total_window_evaluations": int(delta.size),
            "valid_evaluations": int(valid.sum()),
            "invalid_evaluations": int((~valid).sum()),
            "median_absolute_error_bins": float(
                np.median(np.abs(valid_delta))
            ),
            "fraction_exact": float(
                np.mean(valid_delta == 0)
            ),
            "fraction_within_one_bin": float(
                np.mean(np.abs(valid_delta) <= 1)
            ),
            "fraction_narrower_by_two_or_more_bins": float(
                np.mean(valid_delta <= -2)
            ),
            "fraction_wider_by_two_or_more_bins": float(
                np.mean(valid_delta >= 2)
            ),
        })

        for index, row in enumerate(rows):
            values = simulated[:, index]
            valid_values = values[np.isfinite(values)]
            changes = valid_values - baseline[index]

            if len(valid_values):
                median = float(np.median(valid_values))
                q05 = float(np.percentile(valid_values, 5))
                q95 = float(np.percentile(valid_values, 95))
                fraction_exact = float(
                    np.mean(changes == 0)
                )
                fraction_narrower = float(
                    np.mean(changes <= -2)
                )
            else:
                median = ""
                q05 = ""
                q95 = ""
                fraction_exact = ""
                fraction_narrower = ""

            detailed_rows.append({
                "chromosome": row["chromosome"],
                "start_bp": int(row["start_bp"]),
                "end_bp": int(row["end_bp"]),
                "retention": retention,
                "baseline_w80_bins": int(baseline[index]),
                "baseline_contact_signal": baseline_signal[index],
                "replicates": REPLICATES,
                "valid_replicates": int(len(valid_values)),
                "median_thinned_w80_bins": median,
                "thinned_w80_q05_bins": q05,
                "thinned_w80_q95_bins": q95,
                "fraction_exact": fraction_exact,
                "fraction_narrower_by_two_or_more_bins": (
                    fraction_narrower
                ),
            })

    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(detailed_rows[0]),
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(detailed_rows)

    csv_bytes = buffer.getvalue().encode("utf-8")
    result_hash = hashlib.sha256(csv_bytes).hexdigest()

    metadata = {
        "experiment": "conditional_binomial_contact_thinning",
        "source": SOURCE,
        "source_sha256": source_hash,
        "reference_profile": PROFILE,
        "reference_sha256": profile_hash,
        "output_sha256": result_hash,
        "seed": SEED,
        "replicates_per_depth": REPLICATES,
        "retention_fractions": list(RETENTION),
        "genomic_windows": len(rows),
        "w80_definition": "off_diagonal_separations_1_to_64",
        "units": "500_bp_per_bin",
        "missing_bin_status": "not_assessed",
        "interpretation": (
            "Technical sensitivity conditional on observed counts; "
            "not biological replication or calibrated confidence intervals"
        ),
        "summaries": summaries,
    }

    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_meta_path.parent.mkdir(parents=True, exist_ok=True)

    staging_csv = result_path.with_suffix(".partial")
    staging_meta = result_meta_path.with_suffix(".partial")

    staging_csv.write_bytes(csv_bytes)
    staging_meta.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    os.replace(staging_csv, result_path)
    os.replace(staging_meta, result_meta_path)

    volume.commit()

    return {
        "result_path": str(result_path),
        "metadata_path": str(result_meta_path),
        "baseline_validation": "PASS",
        "windows": len(rows),
        "summary": summaries,
    }


@app.local_entrypoint()
def main():
    result = test_downsampling.remote()

    print("\nHCT116 W80 CONTACT-DEPTH SENSITIVITY")
    print("=" * 75)

    print("Baseline validation:", result["baseline_validation"])
    print("Windows:", result["windows"])

    print("\n")
    print(
        f"{'Retention':>10}"
        f"{'Valid':>9}"
        f"{'Invalid':>10}"
        f"{'Exact':>10}"
        f"{'Within 1':>11}"
        f"{'Narrow >=2':>13}"
        f"{'Wide >=2':>11}"
    )

    for item in result["summary"]:
        print(
            f"{item['retention']:>10.0%}"
            f"{item['valid_evaluations']:>9}"
            f"{item['invalid_evaluations']:>10}"
            f"{item['fraction_exact']:>10.1%}"
            f"{item['fraction_within_one_bin']:>11.1%}"
            f"{item['fraction_narrower_by_two_or_more_bins']:>13.1%}"
            f"{item['fraction_wider_by_two_or_more_bins']:>11.1%}"
        )

    print("\nDetailed results:", result["result_path"])
    print("Metadata:", result["metadata_path"])
    print("\nControlled downsampling experiment: PASS")
