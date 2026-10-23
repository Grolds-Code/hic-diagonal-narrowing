"""Read-only comparison of original dense and cloud sparse W80."""

import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-w80-equivalence")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
    .add_local_file(
        "src/hic_diagonal_narrowing/width.py",
        remote_path="/root/width_reference.py",
    )
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile"

RES = 500
WINDOW = 64
MAX_SEP = 64


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def compare_estimators():
    import numpy as np

    spec = importlib.util.spec_from_file_location(
        "width_reference",
        "/root/width_reference.py",
    )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    original_w80 = module.estimate_w80

    source_path = Path(f"/research/matrices/{SOURCE}.npz")
    source_meta_path = Path(f"/research/metadata/{SOURCE}.json")
    profile_path = Path(f"/research/results/{PROFILE}.csv")
    profile_meta_path = Path(f"/research/metadata/{PROFILE}.json")

    source_meta = json.loads(source_meta_path.read_text())
    profile_meta = json.loads(profile_meta_path.read_text())

    source_bytes = source_path.read_bytes()
    profile_bytes = profile_path.read_bytes()

    if hashlib.sha256(source_bytes).hexdigest() != source_meta["sha256"]:
        raise RuntimeError("Source dataset checksum mismatch")

    if hashlib.sha256(profile_bytes).hexdigest() != profile_meta["output_sha256"]:
        raise RuntimeError("Stored W80 profile checksum mismatch")

    with np.load(source_path, allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    rows = list(csv.DictReader(
        io.StringIO(profile_bytes.decode("utf-8"))
    ))

    if len(rows) != 122:
        raise RuntimeError("Unexpected number of W80 windows")

    midpoint_twice = x + y
    separations = (y - x) // RES

    comparisons = []

    for row in rows:
        start_bp = int(row["start_bp"])
        end_bp = int(row["end_bp"])

        if end_bp - start_bp != WINDOW * RES:
            raise RuntimeError("Unexpected window length")

        selected = (
            (midpoint_twice >= 2 * start_bp)
            & (midpoint_twice < 2 * end_bp)
        )

        local_x = x[selected]
        local_y = y[selected]
        local_counts = counts[selected]
        local_dist = separations[selected]

        # The original estimator excludes separation zero.
        non_diagonal = (
            (local_dist >= 1)
            & (local_dist <= MAX_SEP)
        )

        distances = local_dist[non_diagonal]
        values = local_counts[non_diagonal]

        # Sparse calculation using separations 1 through 64.
        distance_signal = np.bincount(
            distances,
            weights=values,
            minlength=MAX_SEP + 1,
        )[1:]

        if distance_signal.sum() <= 0:
            raise RuntimeError("No off-diagonal signal")

        sparse_w80 = int(
            np.searchsorted(
                np.cumsum(distance_signal),
                0.8 * distance_signal.sum(),
            ) + 1
        )

        # Build a padded 192 x 192 reference matrix:
        # 64 left-buffer bins + 64 window bins + 64 right-buffer bins.
        padding_start = start_bp - MAX_SEP * RES

        dense = np.zeros(
            (WINDOW + 2 * MAX_SEP,) * 2,
            dtype=np.float64,
        )

        i = (
            (local_x[non_diagonal] - padding_start) // RES
        ).astype(int)

        j = (
            (local_y[non_diagonal] - padding_start) // RES
        ).astype(int)

        if (
            np.any(i < 0)
            or np.any(j >= dense.shape[0])
            or np.any(i >= j)
        ):
            raise RuntimeError("Invalid padded matrix indexing")

        dense[i, j] = values
        dense[j, i] = values

        dense_w80 = original_w80(
            dense,
            (MAX_SEP, MAX_SEP + WINDOW),
            max_separation_bins=MAX_SEP,
        )

        if dense_w80 != sparse_w80:
            raise RuntimeError(
                f"Dense/sparse disagreement in "
                f"{start_bp}-{end_bp}: "
                f"dense={dense_w80}, sparse={sparse_w80}"
            )

        stored_w80 = int(row["w80_bins"])

        comparisons.append({
            "start_bp": start_bp,
            "end_bp": end_bp,
            "original_dense_w80": dense_w80,
            "corrected_sparse_w80": sparse_w80,
            "previous_modal_w80": stored_w80,
            "diagonal_signal": float(
                local_counts[local_dist == 0].sum()
            ),
        })

    corrected = np.array([
        r["corrected_sparse_w80"]
        for r in comparisons
    ])

    previous = np.array([
        r["previous_modal_w80"]
        for r in comparisons
    ])

    differences = [
        r for r in comparisons
        if (
            r["corrected_sparse_w80"]
            != r["previous_modal_w80"]
        )
    ]

    return {
        "source_checksum": "PASS",
        "profile_checksum": "PASS",
        "windows_tested": len(comparisons),
        "dense_sparse_agreement": "PASS",
        "exact_agreements": len(comparisons),
        "previous_profile_disagreements": len(differences),
        "original_w80_min": int(corrected.min()),
        "original_w80_median": float(np.median(corrected)),
        "original_w80_max": int(corrected.max()),
        "previous_w80_median": float(np.median(previous)),
        "largest_absolute_difference_bins": int(
            np.max(np.abs(corrected - previous))
        ),
        "first_five_differences": differences[:5],
        "lowest_original_w80_windows": sorted(
            comparisons,
            key=lambda r: (
                r["original_dense_w80"],
                r["start_bp"],
            ),
        )[:5],
    }


@app.local_entrypoint()
def main():
    result = compare_estimators.remote()

    print("\nHCT116 W80 NUMERICAL EQUIVALENCE")
    print("=" * 55)

    for key, value in result.items():
        if isinstance(value, list):
            continue
        print(f"{key}: {value}")

    print("\nFIRST FIVE DIFFERENCES FROM STORED PROFILE")
    for item in result["first_five_differences"]:
        print(item)

    print("\nLOWEST W80 WINDOWS UNDER ORIGINAL DEFINITION")
    for item in result["lowest_original_w80_windows"]:
        print(item)

    print("\nOriginal dense vs corrected sparse: PASS")
