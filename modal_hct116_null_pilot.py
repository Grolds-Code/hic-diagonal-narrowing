"""Read-only pilot of a midpoint-depth-preserving Hi-C null."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-null-pilot")

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

REPLICATES = 100
SEED = 20261024
THRESHOLD = 0.20


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=2048,
    timeout=600,
)
def run_null_pilot():
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

    if (
        hashlib.sha256(source_bytes).hexdigest()
        != source_meta["sha256"]
    ):
        raise RuntimeError("Source checksum mismatch")

    if (
        hashlib.sha256(profile_bytes).hexdigest()
        != profile_meta["output_sha256"]
    ):
        raise RuntimeError("Profile checksum mismatch")

    with np.load(io.BytesIO(source_bytes), allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    profile = pd.read_csv(
        io.BytesIO(profile_bytes)
    ).sort_values(
        ["chromosome", "start_bp", "end_bp"]
    ).reset_index(drop=True)

    if len(profile) != 122:
        raise RuntimeError("Expected 122 windows")

    if (
        source_meta["resolution_bp"] != 500
        or source_meta["max_separation_bp"] != 32000
    ):
        raise RuntimeError("Unexpected contact-band geometry")

    if (
        np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
        raise RuntimeError("Expected nonnegative integer counts")

    start = int(source_meta["start_bp"])
    end = int(source_meta["end_bp_exclusive"])
    resolution = 500

    # Midpoint indices have 250 bp spacing because they
    # represent the sum of two 500 bp-aligned bin coordinates.
    midpoint_idx = (
        (x + y - 2 * start) // resolution
    ).astype(int)

    separation = (
        (y - x) // resolution
    ).astype(int)

    included = (
        (separation >= 1)
        & (separation <= 64)
    )

    m = midpoint_idx[included]
    d = separation[included]
    c = counts[included].astype(np.int64)

    midpoint_count = 2 * ((end - start) // resolution)

    if (
        np.any(m < 0)
        or np.any(m >= midpoint_count)
        or np.any(m % 2 != d % 2)
    ):
        raise RuntimeError("Invalid midpoint-distance parity")

    observed = np.zeros(
        (midpoint_count, 65),
        dtype=np.int64,
    )

    np.add.at(observed, (m, d), c)

    if observed[:, 0].sum() != 0:
        raise RuntimeError("Main diagonal was not excluded")

    midpoint_totals = observed.sum(axis=1)

    even_midpoints = np.arange(0, midpoint_count, 2)
    odd_midpoints = np.arange(1, midpoint_count, 2)

    even_distances = np.arange(2, 65, 2)
    odd_distances = np.arange(1, 65, 2)

    even_weights = observed[:, even_distances].sum(axis=0)
    odd_weights = observed[:, odd_distances].sum(axis=0)

    if even_weights.sum() <= 0 or odd_weights.sum() <= 0:
        raise RuntimeError("Insufficient background signal")

    p_even = even_weights / even_weights.sum()
    p_odd = odd_weights / odd_weights.sum()

    starts = (
        2 * (
            (profile["start_bp"].to_numpy(dtype=int) - start)
            // resolution
        )
    )
    ends = (
        2 * (
            (profile["end_bp"].to_numpy(dtype=int) - start)
            // resolution
        )
    )

    if (
        np.any(starts < 0)
        or np.any(ends > midpoint_count)
        or np.any(ends - starts != 128)
    ):
        raise RuntimeError("Unexpected W80 window geometry")

    def profile_from_grid(grid):
        cumulative = np.vstack([
            np.zeros((1, 65), dtype=np.int64),
            np.cumsum(grid, axis=0),
        ])

        distance_signals = (
            cumulative[ends] - cumulative[starts]
        )[:, 1:]

        totals = distance_signals.sum(axis=1)

        if np.any(totals <= 0):
            raise RuntimeError("Uninformative simulated window")

        widths = np.array([
            int(
                np.searchsorted(
                    np.cumsum(bands),
                    0.8 * total,
                ) + 1
            )
            for bands, total in zip(
                distance_signals, totals
            )
        ])

        updated = profile.copy()
        updated["w80_bins"] = widths
        updated["total_contacts"] = totals.astype(float)

        return updated

    baseline = profile_from_grid(observed)

    if not np.array_equal(
        baseline["w80_bins"].to_numpy(),
        profile["w80_bins"].to_numpy(),
    ):
        raise RuntimeError("Baseline W80 reconstruction mismatch")

    if not np.allclose(
        baseline["total_contacts"].to_numpy(),
        profile["total_contacts"].to_numpy(),
    ):
        raise RuntimeError("Baseline signal mismatch")

    baseline_scored, baseline_regions = module.detect_narrowing(
        baseline,
        score_threshold=THRESHOLD,
        min_reference_windows=3,
    )

    if (
        baseline_scored["candidate_call"].sum() != 5
        or len(baseline_regions) != 2
    ):
        raise RuntimeError("Baseline detector regression failed")

    observed_max = float(
        baseline_scored["narrowing_score"].max()
    )

    rng = np.random.default_rng(SEED)

    null_max_scores = []
    null_called_windows = []
    null_merged_regions = []

    for replicate in range(REPLICATES):
        simulated = np.zeros_like(observed)

        simulated[np.ix_(
            even_midpoints, even_distances
        )] = rng.multinomial(
            midpoint_totals[even_midpoints],
            p_even,
        )

        simulated[np.ix_(
            odd_midpoints, odd_distances
        )] = rng.multinomial(
            midpoint_totals[odd_midpoints],
            p_odd,
        )

        if not np.array_equal(
            simulated.sum(axis=1),
            midpoint_totals,
        ):
            raise RuntimeError(
                "Midpoint contact-depth preservation failed"
            )

        simulated_profile = profile_from_grid(simulated)

        scored, regions = module.detect_narrowing(
            simulated_profile,
            score_threshold=THRESHOLD,
            min_reference_windows=3,
        )

        scores = scored[
            "narrowing_score"
        ].to_numpy(dtype=float)

        if not np.all(np.isfinite(scores)):
            raise RuntimeError("Nonfinite simulated score")

        null_max_scores.append(float(scores.max()))
        null_called_windows.append(
            int(scored["candidate_call"].sum())
        )
        null_merged_regions.append(len(regions))

    null_max_scores = np.asarray(null_max_scores)

    return {
        "source_checksum": "PASS",
        "baseline_reconstruction": "PASS",
        "midpoint_depth_preservation": "PASS",
        "windows": len(profile),
        "null_replicates": REPLICATES,
        "seed": SEED,
        "observed_max_score": observed_max,
        "observed_called_windows": 5,
        "observed_merged_regions": 2,
        "null_max_score_percentiles": [
            float(q) for q in np.percentile(
                null_max_scores, [5, 50, 95]
            )
        ],
        "null_max_score_range": [
            float(null_max_scores.min()),
            float(null_max_scores.max()),
        ],
        "null_max_exceedances": int(
            np.count_nonzero(
                null_max_scores >= observed_max
            )
        ),
        "median_null_called_windows": float(
            np.median(null_called_windows)
        ),
        "median_null_merged_regions": float(
            np.median(null_merged_regions)
        ),
        "interpretation": (
            "Exploratory homogeneous-distance null; "
            "not a calibrated biological significance test"
        ),
    }


@app.local_entrypoint()
def main():
    result = run_null_pilot.remote()

    print("\nHCT116 HOMOGENEOUS-DISTANCE NULL PILOT")
    print("=" * 65)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nRead-only null-model pilot: PASS")
