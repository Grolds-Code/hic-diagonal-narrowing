"""Test HCT116 narrowing-call persistence under contact thinning."""

import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-candidate-stability")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4", "pandas==3.1.0")
    .add_local_file(
        "src/hic_diagonal_narrowing/detect.py",
        remote_path="/root/detect_reference.py",
    )
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"
OUTPUT = SOURCE + "_candidate_stability_v1"

RETENTIONS = (0.25, 0.10)
REPLICATES = 100
SEED = 20261024
THRESHOLD = 0.20


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=900,
)
def evaluate_candidate_stability():
    import numpy as np
    import pandas as pd

    spec = importlib.util.spec_from_file_location(
        "detect_reference",
        "/root/detect_reference.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    root = Path("/research")

    source_path = root / "matrices" / f"{SOURCE}.npz"
    profile_path = root / "results" / f"{PROFILE}.csv"

    source_meta = json.loads(
        (root / "metadata" / f"{SOURCE}.json").read_text()
    )
    profile_meta = json.loads(
        (root / "metadata" / f"{PROFILE}.json").read_text()
    )

    source_bytes = source_path.read_bytes()
    profile_bytes = profile_path.read_bytes()

    source_hash = hashlib.sha256(source_bytes).hexdigest()
    profile_hash = hashlib.sha256(profile_bytes).hexdigest()

    if source_hash != source_meta["sha256"]:
        raise RuntimeError("Source checksum mismatch")

    if profile_hash != profile_meta["output_sha256"]:
        raise RuntimeError("Profile checksum mismatch")

    with np.load(source_path, allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    if (
        not np.all(np.isfinite(counts))
        or np.any(counts < 0)
        or not np.all(counts == np.rint(counts))
    ):
        raise RuntimeError("Expected nonnegative integer counts")

    integer_counts = counts.astype(np.int64)

    profile = pd.read_csv(
        io.BytesIO(profile_bytes)
    ).sort_values(
        ["chromosome", "start_bp", "end_bp"]
    ).reset_index(drop=True)

    if len(profile) != 122:
        raise RuntimeError("Unexpected profile length")

    midpoints_twice = x + y
    separations = (y - x) // 500

    window_indices = []
    window_distances = []

    for row in profile.itertuples():
        selected = (
            (midpoints_twice >= 2 * int(row.start_bp))
            & (midpoints_twice < 2 * int(row.end_bp))
            & (separations >= 1)
            & (separations <= 64)
        )

        indices = np.flatnonzero(selected)

        window_indices.append(indices)
        window_distances.append(separations[indices])

    def make_profile(observed_counts):
        widths = []
        totals = []

        for indices, distances in zip(
            window_indices, window_distances
        ):
            signal = np.bincount(
                distances,
                weights=observed_counts[indices],
                minlength=65,
            )[1:]

            total = float(signal.sum())
            totals.append(total)

            if total > 0:
                width = int(
                    np.searchsorted(
                        np.cumsum(signal),
                        0.8 * total,
                    ) + 1
                )
            else:
                width = np.nan

            widths.append(width)

        result = profile.copy()
        result["w80_bins"] = widths
        result["total_contacts"] = totals

        return result

    baseline_profile = make_profile(integer_counts)

    if not np.array_equal(
        baseline_profile["w80_bins"].to_numpy(),
        profile["w80_bins"].to_numpy(),
    ):
        raise RuntimeError("Baseline W80 mismatch")

    if not np.allclose(
        baseline_profile["total_contacts"].to_numpy(),
        profile["total_contacts"].to_numpy(),
    ):
        raise RuntimeError("Baseline contact signal mismatch")

    baseline_scored, baseline_regions = module.detect_narrowing(
        baseline_profile,
        score_threshold=THRESHOLD,
        min_reference_windows=3,
    )

    baseline_call = baseline_scored[
        "candidate_call"
    ].to_numpy(dtype=bool)

    if baseline_call.sum() != 5 or len(baseline_regions) != 2:
        raise RuntimeError(
            "Baseline detection differs from Step 88"
        )

    candidate_indices = np.flatnonzero(baseline_call)

    region_members = []
    region_coordinates = []

    for region in baseline_regions.itertuples():
        members = [
            int(i) for i in candidate_indices
            if (
                int(profile.iloc[i]["start_bp"]) < region.end_bp
                and int(profile.iloc[i]["end_bp"]) > region.start_bp
            )
        ]

        if not members:
            raise RuntimeError("Candidate has no supporting windows")

        region_members.append(members)
        region_coordinates.append([
            int(region.start_bp),
            int(region.end_bp),
        ])

    rng = np.random.default_rng(SEED)
    replicate_rows = []
    summaries = []

    for retention in RETENTIONS:
        persistence = np.zeros(
            len(candidate_indices),
            dtype=int,
        )
        region_persistence = np.zeros(
            len(region_members),
            dtype=int,
        )

        total_calls = []
        new_calls = []
        merged_counts = []

        for replicate in range(REPLICATES):
            thinned_counts = rng.binomial(
                integer_counts,
                retention,
            )

            simulated = make_profile(thinned_counts)

            scored, merged = module.detect_narrowing(
                simulated,
                score_threshold=THRESHOLD,
                min_reference_windows=3,
            )

            calls = scored[
                "candidate_call"
            ].to_numpy(dtype=bool)

            if not np.all(
                np.isfinite(scored["narrowing_score"])
            ):
                raise RuntimeError(
                    "Undefined narrowing score in simulation"
                )

            for position, index in enumerate(candidate_indices):
                persistence[position] += int(calls[index])

            region_detected = [
                bool(np.any(calls[members]))
                for members in region_members
            ]

            region_persistence += np.array(
                region_detected,
                dtype=int,
            )

            called_count = int(calls.sum())
            newly_called = int(
                np.count_nonzero(calls & ~baseline_call)
            )

            total_calls.append(called_count)
            new_calls.append(newly_called)
            merged_counts.append(len(merged))

            replicate_rows.append({
                "retention": retention,
                "replicate": replicate,
                "called_windows": called_count,
                "newly_called_windows": newly_called,
                "baseline_calls_retained": int(
                    np.count_nonzero(calls & baseline_call)
                ),
                "merged_regions": len(merged),
                "region_a_retained": region_detected[0],
                "region_b_retained": region_detected[1],
            })

        window_persistence = []

        for position, index in enumerate(candidate_indices):
            row = profile.iloc[index]

            window_persistence.append({
                "start_bp": int(row["start_bp"]),
                "end_bp": int(row["end_bp"]),
                "baseline_w80_bins": int(row["w80_bins"]),
                "retention_fraction": float(
                    persistence[position] / REPLICATES
                ),
            })

        summaries.append({
            "retention": retention,
            "replicates": REPLICATES,
            "median_called_windows": float(
                np.median(total_calls)
            ),
            "median_newly_called_windows": float(
                np.median(new_calls)
            ),
            "median_merged_regions": float(
                np.median(merged_counts)
            ),
            "region_persistence": [
                float(value / REPLICATES)
                for value in region_persistence
            ],
            "baseline_window_persistence": window_persistence,
        })

    output_path = root / "results" / f"{OUTPUT}.csv"
    metadata_path = root / "metadata" / f"{OUTPUT}.json"

    if output_path.exists() or metadata_path.exists():
        raise RuntimeError("Output already exists")

    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(replicate_rows[0]),
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(replicate_rows)

    output_bytes = buffer.getvalue().encode("utf-8")
    output_hash = hashlib.sha256(output_bytes).hexdigest()

    metadata = {
        "experiment": "candidate_score_under_joint_contact_thinning",
        "source_sha256": source_hash,
        "profile_sha256": profile_hash,
        "output_sha256": output_hash,
        "threshold": THRESHOLD,
        "seed": SEED,
        "replicates_per_depth": REPLICATES,
        "retentions": list(RETENTIONS),
        "baseline_called_windows": int(baseline_call.sum()),
        "baseline_regions": region_coordinates,
        "summaries": summaries,
        "interpretation": (
            "Conditional technical sensitivity; threshold selected "
            "exploratorily; not biological replication or significance"
        ),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    staging_csv = output_path.with_suffix(".partial")
    staging_meta = metadata_path.with_suffix(".partial")

    staging_csv.write_bytes(output_bytes)
    staging_meta.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    os.replace(staging_csv, output_path)
    os.replace(staging_meta, metadata_path)

    volume.commit()

    return {
        "baseline_regions": region_coordinates,
        "baseline_calls": int(baseline_call.sum()),
        "summary": summaries,
        "result_path": str(output_path),
        "metadata_path": str(metadata_path),
    }


@app.local_entrypoint()
def main():
    result = evaluate_candidate_stability.remote()

    print("\nHCT116 NARROWING CANDIDATE STABILITY")
    print("=" * 65)

    print("Baseline windows:", result["baseline_calls"])
    print("Baseline regions:", result["baseline_regions"])

    for summary in result["summary"]:
        print(
            f"\nCONTACT RETENTION: {summary['retention']:.0%}"
        )
        print(
            "Median called windows:",
            summary["median_called_windows"],
        )
        print(
            "Median newly called windows:",
            summary["median_newly_called_windows"],
        )
        print(
            "Median merged regions:",
            summary["median_merged_regions"],
        )
        print(
            "Region persistence:",
            summary["region_persistence"],
        )
        print("Original candidate-window persistence:")

        for window in summary[
            "baseline_window_persistence"
        ]:
            print(
                f"  {window['start_bp']}-{window['end_bp']}: "
                f"{window['retention_fraction']:.1%}"
            )

    print("\nResults:", result["result_path"])
    print("Metadata:", result["metadata_path"])
    print("\nCandidate stability experiment: PASS")
