"""Depth-matched HCT116 W80 uncertainty under joint binomial thinning.

Conditional technical sensitivity, not biological confidence coverage.
"""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-depth-matched-pilot")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
START = 20_000_000
RES = 500
WINDOW = 64
DMAX = 64
REPLICATES = 100
SEED = 20261023


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=1024,
    timeout=180,
)
def load_verified_contacts():
    import numpy as np

    root = Path("/research")
    source = (root / "matrices" / f"{SOURCE}.npz").read_bytes()
    meta = json.loads(
        (root / "metadata" / f"{SOURCE}.json").read_text()
    )

    digest = hashlib.sha256(source).hexdigest()

    if digest != meta["sha256"]:
        raise RuntimeError("Source checksum mismatch")

    with np.load(io.BytesIO(source), allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    if (
        x.shape != y.shape
        or x.shape != counts.shape
        or counts.size != meta["selected_records"]
        or np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
        or np.any(counts > 2**53 - 1)
    ):
        raise RuntimeError("Invalid stored contact values")

    if (
        np.any(x > y)
        or np.any(x % RES != 0)
        or np.any(y % RES != 0)
        or np.any(y - x > DMAX * RES)
        or np.any(x + y < 2 * START)
        or np.any(x + y >= 42_000_000)
    ):
        raise RuntimeError("Unexpected contact geometry")

    return {
        "source_sha256": digest,
        "x": ((x - START) // RES).tolist(),
        "y": ((y - START) // RES).tolist(),
        "counts": counts.astype(np.int64).tolist(),
    }


@app.local_entrypoint()
def main():
    import numpy as np

    from hic_diagonal_narrowing.poisson_uncertainty import (
        cutoff_exposure_w80_interval_from_bands,
    )

    root = Path("research_artifacts/stage1")

    original_intervals = list(csv.DictReader(
        (root / "hct116_conditional_w80_intervals_v1.csv")
        .open(newline="")
    ))

    thin_path = root / "w80_downsampling_v1.csv"
    thin_meta = json.loads(
        (root / "w80_downsampling_v1.metadata.json").read_text()
    )

    if (
        hashlib.sha256(thin_path.read_bytes()).hexdigest()
        != thin_meta["output_sha256"]
    ):
        raise RuntimeError("Original thinning CSV checksum mismatch")

    if (
        thin_meta["seed"] != SEED
        or thin_meta["replicates_per_depth"] != REPLICATES
        or list(thin_meta["retention_fractions"])
        != [0.5, 0.25, 0.1, 0.05, 0.01]
    ):
        raise RuntimeError("Unexpected original thinning configuration")

    previous = {}

    with thin_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (
                float(row["retention"]),
                int(row["start_bp"]),
                int(row["end_bp"]),
            )
            previous[key] = row

    source = load_verified_contacts.remote()

    if source["source_sha256"] != thin_meta["source_sha256"]:
        raise RuntimeError("Data provenance mismatch")

    x = np.asarray(source["x"], dtype=np.int64)
    y = np.asarray(source["y"], dtype=np.int64)
    counts = np.asarray(source["counts"], dtype=np.int64)

    distances = y - x
    midpoint_twice = x + y

    indices = []
    separation_vectors = []
    baseline_w80 = []
    candidate_indices = {}

    for j, row in enumerate(original_intervals):
        left = int(row["start_bp"])
        right = int(row["end_bp"])

        s = (left - START) // RES
        e = (right - START) // RES

        selected = np.flatnonzero(
            (midpoint_twice >= 2 * s)
            & (midpoint_twice < 2 * e)
            & (distances >= 1)
            & (distances <= DMAX)
        )

        local_d = distances[selected]

        bands = np.bincount(
            local_d,
            weights=counts[selected],
            minlength=DMAX + 1,
        )[1:].astype(np.int64)

        total = int(bands.sum())
        observed = int(
            np.searchsorted(
                np.cumsum(bands),
                0.8 * total,
            ) + 1
        )

        if (
            e - s != WINDOW
            or observed != int(row["w80_bins"])
            or total != int(row["total_contacts"])
        ):
            raise RuntimeError("Baseline validation failed")

        indices.append(selected)
        separation_vectors.append(local_d)
        baseline_w80.append(observed)

        if (left, right) == (20_296_000, 20_328_000):
            candidate_indices["Region A"] = j

        if (left, right) == (20_736_000, 20_768_000):
            candidate_indices["Region B"] = j

    if len(indices) != 122 or len(candidate_indices) != 2:
        raise RuntimeError("Incomplete window setup")

    rng = np.random.default_rng(SEED)
    exposure = np.full(DMAX, WINDOW, dtype=np.int64)

    summaries = {}

    # Include the 50% draws to reproduce the exact RNG sequence
    # used in the previous downsampling experiment.
    for retention in (0.5, 0.25, 0.1):
        simulated_w80 = np.zeros(
            (REPLICATES, len(indices)),
            dtype=np.int64,
        )

        lower = np.zeros_like(simulated_w80)
        upper = np.zeros_like(simulated_w80)

        for replicate in range(REPLICATES):
            # ONE joint thinning draw for all overlapping windows.
            thinned = rng.binomial(counts, retention)

            if retention == 0.5:
                continue

            for j, (positions, ds) in enumerate(
                zip(indices, separation_vectors)
            ):
                bands = np.bincount(
                    ds,
                    weights=thinned[positions],
                    minlength=DMAX + 1,
                )[1:].astype(np.int64)

                if bands.sum() == 0:
                    raise RuntimeError(
                        "Empty thinned window: "
                        f"retention={retention}, "
                        f"replicate={replicate}, window={j}"
                    )

                result = cutoff_exposure_w80_interval_from_bands(
                    bands,
                    exposure,
                    confidence_level=0.95,
                )

                simulated_w80[replicate, j] = int(
                    result["observed_w80"]
                )
                lower[replicate, j] = int(
                    result["cutoff_lower"]
                )
                upper[replicate, j] = int(
                    result["cutoff_upper"]
                )

        if retention == 0.5:
            continue

        # Exact regression against the original thinning experiment.
        checked_quantiles = 0

        for j, row in enumerate(original_intervals):
            key = (
                retention,
                int(row["start_bp"]),
                int(row["end_bp"]),
            )

            reference = previous[key]

            observed_q05 = float(np.percentile(
                simulated_w80[:, j], 5
            ))
            observed_q95 = float(np.percentile(
                simulated_w80[:, j], 95
            ))

            if not (
                np.isclose(
                    observed_q05,
                    float(reference["thinned_w80_q05_bins"]),
                    atol=1e-9,
                )
                and np.isclose(
                    observed_q95,
                    float(reference["thinned_w80_q95_bins"]),
                    atol=1e-9,
                )
            ):
                raise RuntimeError(
                    f"Original thinning regression mismatch: {key}"
                )

            checked_quantiles += 1

        baseline = np.asarray(
            baseline_w80, dtype=np.int64
        )[None, :]

        included = (
            (lower <= baseline) & (baseline <= upper)
        )

        spans = upper - lower

        named = {}

        for name, j in candidate_indices.items():
            named[name] = {
                "baseline_w80_bins": int(baseline_w80[j]),
                "conditional_inclusion_fraction": float(
                    included[:, j].mean()
                ),
                "median_interval_span_bins": float(
                    np.median(spans[:, j])
                ),
                "median_thinned_w80_bins": float(
                    np.median(simulated_w80[:, j])
                ),
            }

        summaries[str(retention)] = {
            "retention": retention,
            "replicates": REPLICATES,
            "windows": len(indices),
            "original_thinning_quantiles_reproduced": (
                checked_quantiles
            ),
            "conditional_baseline_inclusion_fraction": float(
                included.mean()
            ),
            "median_depth_matched_interval_span_bins": float(
                np.median(spans)
            ),
            "p90_depth_matched_interval_span_bins": float(
                np.quantile(spans, 0.9)
            ),
            "candidates": named,
        }

    result = {
        "experiment": "depth_matched_poisson_intervals_under_joint_thinning",
        "source_sha256": source["source_sha256"],
        "seed": SEED,
        "replicates_per_depth": REPLICATES,
        "summaries": summaries,
        "interpretation": (
            "Conditional inclusion of full-depth empirical W80 "
            "under joint binomial thinning. Not biological "
            "confidence-interval calibration."
        ),
    }

    output = root / "hct116_depth_matched_intervals_v1.json"
    output.write_text(json.dumps(result, indent=2) + "\n")

    print("\nHCT116 DEPTH-MATCHED UNCERTAINTY PILOT")
    print("=" * 65)
    print(json.dumps(result, indent=2))
    print("\nSaved:", output)
    print("Original Modal research data unchanged.")
