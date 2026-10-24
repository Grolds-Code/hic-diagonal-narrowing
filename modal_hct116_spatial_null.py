"""Spatially stratified HCT116 contact-distance null pilot."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-spatial-null")

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

SEED = 20261024
REPLICATES = 100
BLOCK_BP = 100_000
THRESHOLD = 0.20


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=2048,
    timeout=600,
)
def evaluate_spatial_null():
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
        raise RuntimeError("Unexpected profile length")

    if (
        source_meta["resolution_bp"] != 500
        or source_meta["max_separation_bp"] != 32000
    ):
        raise RuntimeError("Unexpected source geometry")

    if (
        np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
        raise RuntimeError(
            "Expected nonnegative integer contact counts"
        )

    start = int(source_meta["start_bp"])
    end = int(source_meta["end_bp_exclusive"])
    resolution = 500

    midpoint_idx = (
        (x + y - 2 * start) // resolution
    ).astype(int)

    separation = (
        (y - x) // resolution
    ).astype(int)

    eligible = (
        (separation >= 1)
        & (separation <= 64)
    )

    midpoints = midpoint_idx[eligible]
    distances = separation[eligible]
    values = counts[eligible].astype(np.int64)

    midpoint_count = 2 * ((end - start) // resolution)
    block_size = 2 * (BLOCK_BP // resolution)

    if midpoint_count % block_size:
        raise RuntimeError("Incomplete spatial blocks")

    n_blocks = midpoint_count // block_size

    if n_blocks != 10:
        raise RuntimeError("Expected ten 100 kb blocks")

    observed = np.zeros(
        (midpoint_count, 65),
        dtype=np.int64,
    )

    if (
        np.any(midpoints < 0)
        or np.any(midpoints >= midpoint_count)
        or np.any(midpoints % 2 != distances % 2)
    ):
        raise RuntimeError("Invalid midpoint geometry")

    np.add.at(
        observed,
        (midpoints, distances),
        values,
    )

    midpoint_totals = observed.sum(axis=1)

    if midpoint_totals.sum() != values.sum():
        raise RuntimeError("Contact count mismatch")

    window_starts = (
        2 * (
            profile["start_bp"].to_numpy(dtype=int) - start
        ) // resolution
    )

    window_ends = (
        2 * (
            profile["end_bp"].to_numpy(dtype=int) - start
        ) // resolution
    )

    if (
        np.any(window_starts < 0)
        or np.any(window_ends > midpoint_count)
        or np.any(window_ends - window_starts != 128)
    ):
        raise RuntimeError("Invalid sliding-window geometry")

    def make_profile(grid):
        cumulative = np.vstack([
            np.zeros((1, 65), dtype=np.int64),
            np.cumsum(grid, axis=0),
        ])

        bands = (
            cumulative[window_ends]
            - cumulative[window_starts]
        )[:, 1:]

        totals = bands.sum(axis=1)

        if np.any(totals <= 0):
            raise RuntimeError("Window without contact signal")

        widths = np.array([
            int(
                np.searchsorted(
                    np.cumsum(signal),
                    0.8 * total,
                ) + 1
            )
            for signal, total in zip(bands, totals)
        ])

        result = profile.copy()
        result["w80_bins"] = widths
        result["total_contacts"] = totals.astype(float)

        return result

    baseline = make_profile(observed)

    if not np.array_equal(
        baseline["w80_bins"].to_numpy(),
        profile["w80_bins"].to_numpy(),
    ):
        raise RuntimeError("Baseline W80 mismatch")

    if not np.allclose(
        baseline["total_contacts"].to_numpy(),
        profile["total_contacts"].to_numpy(),
    ):
        raise RuntimeError("Baseline signal mismatch")

    observed_scored, observed_regions = (
        module.detect_narrowing(
            baseline,
            score_threshold=THRESHOLD,
            min_reference_windows=3,
        )
    )

    if (
        int(observed_scored["candidate_call"].sum()) != 5
        or len(observed_regions) != 2
    ):
        raise RuntimeError("Detector baseline mismatch")

    observed_max = float(
        observed_scored["narrowing_score"].max()
    )

    # Estimate a separate contact-distance distribution
    # for each 100 kb block and midpoint parity.
    strata = []

    for block in range(n_blocks):
        lower = block * block_size
        upper = (block + 1) * block_size

        for parity in (0, 1):
            mids = np.arange(
                lower + parity,
                upper,
                2,
            )

            distances_for_parity = (
                np.arange(2, 65, 2)
                if parity == 0
                else np.arange(1, 65, 2)
            )

            weights = observed[np.ix_(
                mids,
                distances_for_parity,
            )].sum(axis=0)

            if weights.sum() <= 0:
                raise RuntimeError(
                    f"Empty distance background: "
                    f"block={block}, parity={parity}"
                )

            probabilities = weights / weights.sum()

            strata.append((
                mids,
                distances_for_parity,
                probabilities,
            ))

    rng = np.random.default_rng(SEED)

    null_max_scores = []
    null_called_windows = []
    null_region_counts = []
    null_peak_starts = []

    for _ in range(REPLICATES):
        simulated = np.zeros_like(observed)

        for mids, distances_for_parity, probabilities in strata:
            simulated[np.ix_(
                mids,
                distances_for_parity,
            )] = rng.multinomial(
                midpoint_totals[mids],
                probabilities,
            )

        if not np.array_equal(
            simulated.sum(axis=1),
            midpoint_totals,
        ):
            raise RuntimeError(
                "Midpoint contact depth was not preserved"
            )

        simulated_profile = make_profile(simulated)

        scored, regions = module.detect_narrowing(
            simulated_profile,
            score_threshold=THRESHOLD,
            min_reference_windows=3,
        )

        scores = scored[
            "narrowing_score"
        ].to_numpy(dtype=float)

        if not np.all(np.isfinite(scores)):
            raise RuntimeError("Nonfinite null score")

        index = int(np.argmax(scores))

        null_max_scores.append(float(scores[index]))
        null_peak_starts.append(
            int(scored.iloc[index]["start_bp"])
        )
        null_called_windows.append(
            int(scored["candidate_call"].sum())
        )
        null_region_counts.append(len(regions))

    maxima = np.asarray(null_max_scores)

    # Diagnose possible artifacts at artificial
    # 100 kb background-block boundaries.
    internal_seams = np.arange(
        start + BLOCK_BP,
        end,
        BLOCK_BP,
    )

    peak_centers = (
        np.asarray(null_peak_starts) + 16_000
    )

    closest_seam = np.min(
        np.abs(
            peak_centers[:, None]
            - internal_seams[None, :]
        ),
        axis=1,
    )

    return {
        "source_checksum": "PASS",
        "baseline_reconstruction": "PASS",
        "midpoint_depth_preservation": "PASS",
        "background_blocks": n_blocks,
        "block_width_bp": BLOCK_BP,
        "midpoint_parity_strata": len(strata),
        "windows": len(profile),
        "replicates": REPLICATES,
        "seed": SEED,
        "observed_max_score": observed_max,
        "observed_called_windows": 5,
        "observed_merged_regions": 2,
        "null_maximum_range": [
            float(maxima.min()),
            float(maxima.max()),
        ],
        "null_maximum_percentiles": [
            float(v)
            for v in np.percentile(
                maxima, [5, 50, 95]
            )
        ],
        "null_maxima_at_least_observed": int(
            np.count_nonzero(maxima >= observed_max)
        ),
        "median_null_called_windows": float(
            np.median(null_called_windows)
        ),
        "median_null_merged_regions": float(
            np.median(null_region_counts)
        ),
        "fraction_null_peaks_within_32kb_of_seam": float(
            np.mean(closest_seam <= 32_000)
        ),
        "interpretation": (
            "Exploratory block-stratified null. "
            "Not calibrated for biological significance."
        ),
    }


@app.local_entrypoint()
def main():
    result = evaluate_spatial_null.remote()

    print("\nHCT116 SPATIALLY STRATIFIED NULL PILOT")
    print("=" * 65)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nSpatial null-model pilot: PASS")
