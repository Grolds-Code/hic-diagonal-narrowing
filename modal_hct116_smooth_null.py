"""Read-only pilot of smoothly varying HCT116 distance nulls."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-smooth-null")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy==2.5.4",
        "pandas==3.1.0",
        "scipy==1.18.1",
    )
    .add_local_file(
        "src/hic_diagonal_narrowing/detect.py",
        remote_path="/root/detect_reference.py",
    )
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"

SEED = 20261024
REPLICATES = 100
SIGMA_BP = (50_000, 100_000, 200_000)
THRESHOLD = 0.20


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=900,
)
def run_smooth_null():
    import numpy as np
    import pandas as pd
    from scipy.ndimage import gaussian_filter1d

    spec = importlib.util.spec_from_file_location(
        "detect_reference",
        "/root/detect_reference.py",
    )
    detector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(detector)

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
        raise RuntimeError("Contact source checksum mismatch")

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
        raise RuntimeError("Expected 122 W80 windows")

    if (
        source_meta["resolution_bp"] != 500
        or source_meta["max_separation_bp"] != 32000
    ):
        raise RuntimeError("Unexpected matrix geometry")

    if (
        np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
        raise RuntimeError("Invalid contact counts")

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

    m = midpoint_idx[eligible]
    d = separation[eligible]
    c = counts[eligible].astype(np.int64)

    midpoint_count = 2 * ((end - start) // resolution)

    if (
        np.any(m < 0)
        or np.any(m >= midpoint_count)
        or np.any(m % 2 != d % 2)
    ):
        raise RuntimeError("Invalid midpoint-distance geometry")

    observed = np.zeros(
        (midpoint_count, 65),
        dtype=np.int64,
    )

    np.add.at(observed, (m, d), c)

    midpoint_totals = observed.sum(axis=1)

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
        raise RuntimeError("Unexpected W80 window geometry")

    def calculate_profile(grid):
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
            raise RuntimeError("Uninformative simulated window")

        widths = np.array([
            int(
                np.searchsorted(
                    np.cumsum(signal),
                    0.8 * total,
                ) + 1
            )
            for signal, total in zip(bands, totals)
        ])

        updated = profile.copy()
        updated["w80_bins"] = widths
        updated["total_contacts"] = totals.astype(float)

        return updated

    baseline = calculate_profile(observed)

    if not np.array_equal(
        baseline["w80_bins"].to_numpy(),
        profile["w80_bins"].to_numpy(),
    ):
        raise RuntimeError("Baseline W80 reconstruction failed")

    if not np.allclose(
        baseline["total_contacts"].to_numpy(),
        profile["total_contacts"].to_numpy(),
    ):
        raise RuntimeError("Baseline signal mismatch")

    observed_scores, observed_regions = (
        detector.detect_narrowing(
            baseline,
            score_threshold=THRESHOLD,
            min_reference_windows=3,
        )
    )

    if (
        int(observed_scores["candidate_call"].sum()) != 5
        or len(observed_regions) != 2
    ):
        raise RuntimeError("Baseline detector regression failed")

    observed_max = float(
        observed_scores["narrowing_score"].max()
    )

    window_centers = (
        profile["start_bp"].to_numpy(dtype=int)
        + profile["end_bp"].to_numpy(dtype=int)
    ) / 2

    seams = np.arange(
        start + 100_000,
        end,
        100_000,
    )

    nearest_seam_distance = np.min(
        np.abs(
            window_centers[:, None]
            - seams[None, :]
        ),
        axis=1,
    )

    seam_zone = nearest_seam_distance <= 32_000
    geometric_fraction = float(np.mean(seam_zone))

    summaries = []

    for sigma_bp in SIGMA_BP:
        sigma_in_midpoint_steps = sigma_bp / resolution

        # Smooth separately within each midpoint parity.
        # The physical spacing between consecutive positions
        # in either parity subset is 500 bp.
        strata = []

        for parity in (0, 1):
            midpoint_ids = np.arange(
                parity,
                midpoint_count,
                2,
            )

            permitted_distances = (
                np.arange(2, 65, 2)
                if parity == 0
                else np.arange(1, 65, 2)
            )

            local_counts = observed[np.ix_(
                midpoint_ids,
                permitted_distances,
            )].astype(float)

            smoothed = gaussian_filter1d(
                local_counts,
                sigma=sigma_in_midpoint_steps,
                axis=0,
                mode="reflect",
            )

            denominators = smoothed.sum(axis=1)

            if np.any(denominators <= 0):
                raise RuntimeError(
                    f"No smooth background at sigma={sigma_bp}"
                )

            probabilities = (
                smoothed / denominators[:, None]
            )

            if not np.allclose(
                probabilities.sum(axis=1),
                1.0,
            ):
                raise RuntimeError(
                    "Simulated multinomial probabilities invalid"
                )

            strata.append((
                midpoint_ids,
                permitted_distances,
                probabilities,
            ))

        rng = np.random.default_rng(
            np.random.SeedSequence([SEED, sigma_bp])
        )

        maxima = []
        called_counts = []
        merged_counts = []
        peak_near_seam = []

        for _ in range(REPLICATES):
            simulated = np.zeros_like(observed)

            for ids, distances, probabilities in strata:
                simulated[np.ix_(
                    ids,
                    distances,
                )] = rng.multinomial(
                    midpoint_totals[ids],
                    probabilities,
                )

            if not np.array_equal(
                simulated.sum(axis=1),
                midpoint_totals,
            ):
                raise RuntimeError(
                    "Per-midpoint depth preservation failed"
                )

            simulated_profile = calculate_profile(simulated)

            scored, regions = detector.detect_narrowing(
                simulated_profile,
                score_threshold=THRESHOLD,
                min_reference_windows=3,
            )

            scores = scored[
                "narrowing_score"
            ].to_numpy(dtype=float)

            if not np.all(np.isfinite(scores)):
                raise RuntimeError("Nonfinite null-model scores")

            peak_index = int(np.argmax(scores))

            maxima.append(float(scores[peak_index]))
            called_counts.append(
                int(scored["candidate_call"].sum())
            )
            merged_counts.append(len(regions))
            peak_near_seam.append(
                bool(seam_zone[peak_index])
            )

        maxima = np.asarray(maxima)

        summaries.append({
            "sigma_bp": sigma_bp,
            "null_max_range": [
                float(maxima.min()),
                float(maxima.max()),
            ],
            "null_max_percentiles_5_50_95": [
                float(value)
                for value in np.percentile(
                    maxima, [5, 50, 95]
                )
            ],
            "maxima_at_least_observed": int(
                np.count_nonzero(
                    maxima >= observed_max
                )
            ),
            "median_called_windows": float(
                np.median(called_counts)
            ),
            "median_merged_regions": float(
                np.median(merged_counts)
            ),
            "peak_near_100kb_seam_fraction": float(
                np.mean(peak_near_seam)
            ),
        })

    return {
        "verification": "PASS",
        "observed_max_score": observed_max,
        "observed_called_windows": 5,
        "observed_merged_regions": 2,
        "window_count": len(profile),
        "replicates_per_scale": REPLICATES,
        "seed": SEED,
        "geometric_seam_zone_fraction": geometric_fraction,
        "smoothing_results": summaries,
        "significance": "NOT_CALIBRATED",
    }


@app.local_entrypoint()
def main():
    result = run_smooth_null.remote()

    print("\nHCT116 SMOOTH SPATIAL NULL PILOT")
    print("=" * 65)

    print("Verification:", result["verification"])
    print("Observed maximum:", result["observed_max_score"])
    print("Windows:", result["window_count"])
    print(
        "Simulations per scale:",
        result["replicates_per_scale"],
    )
    print(
        "Geometric seam-zone fraction:",
        result["geometric_seam_zone_fraction"],
    )

    for item in result["smoothing_results"]:
        print(
            f"\nGAUSSIAN SIGMA: {item['sigma_bp']:,} bp"
        )

        for key, value in item.items():
            if key != "sigma_bp":
                print(f"{key}: {value}")

    print("\nSignificance:", result["significance"])
    print("Read-only smooth-null pilot: PASS")
