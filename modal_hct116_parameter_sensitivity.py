"""Stage 1 HCT116 parameter sensitivity.

Recomputes midpoint-assigned contact profiles and original
detector references. Exploratory sensitivity, not inference.
"""

import hashlib
import importlib.util
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-stage1-parameter-sensitivity")
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

WINDOWS = (32, 64, 96)
DISTANCES = (32, 48, 64)
QUANTILES = (0.70, 0.80, 0.90)

STEP = 16
THRESHOLD = 0.2

CENTERS = {
    "Region A": 20_312_000,
    "Region B": 20_752_000,
}


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def analyze():
    import numpy as np
    import pandas as pd

    spec = importlib.util.spec_from_file_location(
        "detect_reference",
        "/root/detect_reference.py",
    )
    detector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(detector)

    root = Path("/research")

    source_bytes = (
        root / "matrices" / f"{SOURCE}.npz"
    ).read_bytes()

    profile_bytes = (
        root / "results" / f"{PROFILE}.csv"
    ).read_bytes()

    source_meta = json.loads(
        (root / "metadata" / f"{SOURCE}.json").read_text()
    )
    profile_meta = json.loads(
        (root / "metadata" / f"{PROFILE}.json").read_text()
    )

    source_sha = hashlib.sha256(source_bytes).hexdigest()
    profile_sha = hashlib.sha256(profile_bytes).hexdigest()

    if (
        source_sha != source_meta["sha256"]
        or profile_sha != profile_meta["output_sha256"]
        or profile_meta["source_sha256"] != source_sha
    ):
        raise RuntimeError("Provenance verification failed")

    start = int(source_meta["start_bp"])
    end = int(source_meta["end_bp_exclusive"])
    resolution = int(source_meta["resolution_bp"])

    if (
        start != 20_000_000
        or end != 21_000_000
        or resolution != 500
    ):
        raise RuntimeError("Unexpected source geometry")

    n_bins = (end - start) // resolution

    with np.load(
        io.BytesIO(source_bytes),
        allow_pickle=False,
    ) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        values = data["counts"]

    if (
        x.shape != y.shape
        or x.shape != values.shape
        or len(values) != source_meta["selected_records"]
        or np.any(~np.isfinite(values))
        or np.any(values < 0)
        or np.any(values != np.rint(values))
    ):
        raise RuntimeError("Invalid source records")

    counts = values.astype(np.int64)

    if np.any(counts.astype(float) != values):
        raise RuntimeError("Integer count conversion failed")

    if (
        np.any((x - start) % resolution != 0)
        or np.any((y - start) % resolution != 0)
        or np.any(y < x)
    ):
        raise RuntimeError("Invalid contact coordinates")

    distances = (y - x) // resolution
    midpoints = (x + y - 2 * start) // resolution

    selected = (
        (distances >= 1)
        & (distances <= 64)
    )

    m = midpoints[selected].astype(int)
    d = distances[selected].astype(int)
    c = counts[selected]

    if (
        np.any(m < 0)
        or np.any(m >= 2 * n_bins)
        or np.any(m % 2 != d % 2)
    ):
        raise RuntimeError("Invalid midpoint-distance geometry")

    grid = np.zeros(
        (2 * n_bins, 65),
        dtype=np.int64,
    )
    np.add.at(grid, (m, d), c)

    cumulative = np.vstack([
        np.zeros((1, 65), dtype=np.int64),
        np.cumsum(grid, axis=0),
    ])

    original = pd.read_csv(io.BytesIO(profile_bytes))

    if len(original) != 122:
        raise RuntimeError("Incorrect baseline profile length")

    results = []
    baseline_verified = False

    for window_bins in WINDOWS:
        start_bins = np.arange(
            0,
            n_bins - window_bins + 1,
            STEP,
            dtype=int,
        )

        window_starts = start + resolution * start_bins
        window_ends = (
            window_starts + resolution * window_bins
        )

        if np.any(window_ends > end):
            raise RuntimeError("Window outside analysis region")

        bands_full = (
            cumulative[2 * (start_bins + window_bins)]
            - cumulative[2 * start_bins]
        )[:, 1:]

        # Under complete valid-bin support, midpoint
        # geometry gives exactly window_bins eligible
        # pairs at every separation.
        if np.any(bands_full < 0):
            raise RuntimeError("Negative distance-band count")

        for dmax in DISTANCES:
            if dmax > window_bins:
                continue

            bands = bands_full[:, :dmax]
            totals = bands.sum(axis=1)

            if np.any(totals <= 0):
                raise RuntimeError(
                    "Window lacks informative contacts"
                )

            for q in QUANTILES:
                widths = np.array([
                    int(
                        np.searchsorted(
                            np.cumsum(row),
                            q * total,
                        ) + 1
                    )
                    for row, total in zip(
                        bands, totals
                    )
                ], dtype=int)

                if (
                    np.any(widths < 1)
                    or np.any(widths > dmax)
                ):
                    raise RuntimeError("Invalid Wq estimate")

                # The detector's required w80_bins field
                # stores Wq in non-80% configurations.
                simulated = pd.DataFrame({
                    "chromosome": "chr10",
                    "start_bp": window_starts,
                    "end_bp": window_ends,
                    "w80_bins": widths,
                    "total_contacts": totals.astype(float),
                    "qc_status": "unassessed",
                })

                scored, regions = detector.detect_narrowing(
                    simulated,
                    score_threshold=THRESHOLD,
                    min_reference_windows=3,
                )

                candidates = {}

                for name, center in CENTERS.items():
                    target_start = (
                        center - window_bins * resolution // 2
                    )
                    target_end = (
                        target_start + window_bins * resolution
                    )

                    rows = scored.loc[
                        (scored["start_bp"] == target_start)
                        & (scored["end_bp"] == target_end)
                    ]

                    if len(rows) != 1:
                        raise RuntimeError(
                            f"Candidate center not aligned: {name}"
                        )

                    r = rows.iloc[0]

                    if (
                        not np.isfinite(r["narrowing_score"])
                        or not np.isfinite(
                            r["reference_w80_bins"]
                        )
                    ):
                        raise RuntimeError(
                            "Nonfinite candidate comparison"
                        )

                    candidates[name] = {
                        "start_bp": int(target_start),
                        "end_bp": int(target_end),
                        "width_bins": int(r["w80_bins"]),
                        "reference_width_bins": float(
                            r["reference_w80_bins"]
                        ),
                        "score": float(r["narrowing_score"]),
                        "called": bool(r["candidate_call"]),
                        "reference_windows": int(
                            r["reference_windows"]
                        ),
                        "total_contacts": int(
                            r["total_contacts"]
                        ),
                    }

                if (
                    window_bins == 64
                    and dmax == 64
                    and q == 0.80
                ):
                    if (
                        not np.array_equal(
                            widths,
                            original["w80_bins"].to_numpy(),
                        )
                        or not np.allclose(
                            totals,
                            original["total_contacts"].to_numpy(),
                            rtol=0,
                            atol=1e-9,
                        )
                        or int(scored["candidate_call"].sum()) != 5
                        or len(regions) != 2
                    ):
                        raise RuntimeError(
                            "Original baseline reconstruction failed"
                        )

                    expected = {
                        "Region A": np.log2(17 / 14),
                        "Region B": np.log2(17 / 13),
                    }

                    for name, score in expected.items():
                        if not np.isclose(
                            candidates[name]["score"],
                            score,
                            rtol=0,
                            atol=1e-12,
                        ):
                            raise RuntimeError(
                                f"Baseline score mismatch: {name}"
                            )

                    baseline_verified = True

                results.append({
                    "window_bins": window_bins,
                    "window_bp": window_bins * resolution,
                    "max_separation_bins": dmax,
                    "max_separation_bp": dmax * resolution,
                    "cumulative_fraction": q,
                    "windows_scanned": len(scored),
                    "called_windows": int(
                        scored["candidate_call"].sum()
                    ),
                    "merged_regions": len(regions),
                    "maximum_score": float(
                        scored["narrowing_score"].max()
                    ),
                    "candidates": candidates,
                })

    if not baseline_verified or len(results) != 21:
        raise RuntimeError("Incomplete sensitivity design")

    return {
        "analysis": "hct116_stage1_parameter_sensitivity",
        "source_sha256": source_sha,
        "profile_sha256": profile_sha,
        "baseline_reconstruction": "PASS",
        "configurations": len(results),
        "window_sizes_bins": list(WINDOWS),
        "maximum_separations_bins": list(DISTANCES),
        "cumulative_fractions": list(QUANTILES),
        "scan_step_bins": STEP,
        "score_threshold": THRESHOLD,
        "results": results,
        "limitations": [
            "Candidate centers selected from observed baseline data",
            "Comparisons are descriptive and exploratory",
            "Wq values across different q or Dmax have different definitions",
            "Original QC status remains unassessed",
            "Independent genomic missing-bin validity not established",
            "No biological significance or error-rate calibration",
        ],
    }


@app.local_entrypoint()
def main():
    result = analyze.remote()

    destination = Path(
        "research_artifacts/stage1/"
        "hct116_parameter_sensitivity_v1.json"
    )

    if destination.exists():
        raise RuntimeError("Refusing to overwrite artifact")

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    destination.write_text(
        json.dumps(result, indent=2) + "\n"
    )

    print("\nHCT116 STAGE 1 PARAMETER SENSITIVITY")
    print("=" * 76)
    print(
        "Baseline reconstruction:",
        result["baseline_reconstruction"],
    )
    print(
        "Configurations:",
        result["configurations"],
    )
    print()

    for item in result["results"]:
        a = item["candidates"]["Region A"]
        b = item["candidates"]["Region B"]

        print(
            f"W={item['window_bp']//1000:2d}kb "
            f"D={item['max_separation_bp']//1000:2d}kb "
            f"q={item['cumulative_fraction']:.2f} | "
            f"A={a['score']:+.3f} "
            f"({'Y' if a['called'] else 'N'}) | "
            f"B={b['score']:+.3f} "
            f"({'Y' if b['called'] else 'N'}) | "
            f"called={item['called_windows']}"
        )

    print()
    print("Saved:", destination)
    print("Original HCT116 data unchanged.")
    print("Biological significance: NOT CALIBRATED")
