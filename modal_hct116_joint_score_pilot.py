"""Joint HCT116 candidate-score sensitivity under contact thinning.

Preserves shared random draws across overlapping genomic windows.
Uses the original detector and its original reference selection.

Technical sensitivity only, not a biological significance test.
"""

import hashlib
import importlib.util
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-joint-score-pilot")

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

START = 20_000_000
RES = 500
DMAX = 64
WINDOW = 64

THRESHOLD = 0.2
REPLICATES = 200
SEED = 20261028
RETENTIONS = (0.50, 0.25, 0.10)

PEAKS = {
    "Region A": (20_296_000, 20_328_000),
    "Region B": (20_736_000, 20_768_000),
}


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=900,
)
def evaluate_joint_scores():
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

    source_hash = hashlib.sha256(source_bytes).hexdigest()
    profile_hash = hashlib.sha256(profile_bytes).hexdigest()

    if (
        source_hash != source_meta["sha256"]
        or profile_hash != profile_meta["output_sha256"]
        or profile_meta["source_sha256"] != source_hash
    ):
        raise RuntimeError("Source provenance mismatch")

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
        raise RuntimeError("Invalid original count data")

    counts = values.astype(np.int64)

    if not np.array_equal(counts.astype(float), values):
        raise RuntimeError("Integer conversion mismatch")

    d = (y - x) // RES
    midpoint_twice = x + y

    profile = pd.read_csv(io.BytesIO(profile_bytes))

    if len(profile) != 122:
        raise RuntimeError("Unexpected profile count")

    indices = []
    separations = []
    peak_indices = {}

    for j, row in enumerate(profile.itertuples()):
        left = int(row.start_bp)
        right = int(row.end_bp)

        if right - left != WINDOW * RES:
            raise RuntimeError("Unexpected window width")

        selected = np.flatnonzero(
            (midpoint_twice >= 2 * left)
            & (midpoint_twice < 2 * right)
            & (d >= 1)
            & (d <= DMAX)
        )

        indices.append(selected)
        separations.append(d[selected].astype(int))

        for name, interval in PEAKS.items():
            if (left, right) == interval:
                peak_indices[name] = j

    if len(peak_indices) != 2:
        raise RuntimeError("Candidate windows missing")

    def reconstruct(sample_counts):
        widths = np.zeros(len(indices), dtype=int)
        totals = np.zeros(len(indices), dtype=np.int64)

        for j, (positions, ds) in enumerate(
            zip(indices, separations)
        ):
            bands = np.bincount(
                ds,
                weights=sample_counts[positions],
                minlength=DMAX + 1,
            )[1:]

            total = float(bands.sum())

            if total <= 0:
                raise RuntimeError(
                    "Uninformative thinned window"
                )

            widths[j] = int(
                np.searchsorted(
                    np.cumsum(bands),
                    0.80 * total,
                ) + 1
            )
            totals[j] = int(total)

        simulated = profile.copy()
        simulated["w80_bins"] = widths
        simulated["total_contacts"] = totals

        return simulated

    baseline = reconstruct(counts)

    if (
        not np.array_equal(
            baseline["w80_bins"].to_numpy(),
            profile["w80_bins"].to_numpy(),
        )
        or not np.allclose(
            baseline["total_contacts"].to_numpy(),
            profile["total_contacts"].to_numpy(),
        )
    ):
        raise RuntimeError("Baseline W80 reproduction failed")

    baseline_scored, baseline_regions = (
        detector.detect_narrowing(
            baseline,
            score_threshold=THRESHOLD,
            min_reference_windows=3,
        )
    )

    if (
        int(baseline_scored["candidate_call"].sum()) != 5
        or len(baseline_regions) != 2
    ):
        raise RuntimeError("Detector baseline mismatch")

    baseline_scores = {
        name: float(
            baseline_scored.iloc[j]["narrowing_score"]
        )
        for name, j in peak_indices.items()
    }

    if not (
        np.isclose(
            baseline_scores["Region A"], 0.2801079191927352
        )
        and np.isclose(
            baseline_scores["Region B"],
            np.log2(17 / 13),
        )
    ):
        raise RuntimeError("Unexpected candidate baseline score")

    rng = np.random.default_rng(SEED)
    summaries = {}

    for retention in RETENTIONS:
        scores_a = []
        scores_b = []
        reference_a = []
        reference_b = []
        called_counts = []
        both_called = []

        for replicate in range(REPLICATES):
            # One joint genomic draw for every window.
            thinned = rng.binomial(counts, retention)

            simulated = reconstruct(thinned)

            scored, _ = detector.detect_narrowing(
                simulated,
                score_threshold=THRESHOLD,
                min_reference_windows=3,
            )

            a = scored.iloc[peak_indices["Region A"]]
            b = scored.iloc[peak_indices["Region B"]]

            sa = float(a["narrowing_score"])
            sb = float(b["narrowing_score"])

            if not (
                np.isfinite(sa)
                and np.isfinite(sb)
                and np.isfinite(a["reference_w80_bins"])
                and np.isfinite(b["reference_w80_bins"])
            ):
                raise RuntimeError(
                    "Nonfinite candidate score or reference"
                )

            scores_a.append(sa)
            scores_b.append(sb)
            reference_a.append(float(a["reference_w80_bins"]))
            reference_b.append(float(b["reference_w80_bins"]))

            both_called.append(
                bool(a["candidate_call"])
                and bool(b["candidate_call"])
            )

            called_counts.append(
                int(scored["candidate_call"].sum())
            )

        scores_a = np.asarray(scores_a)
        scores_b = np.asarray(scores_b)
        delta = scores_b - scores_a

        def summarize(values):
            return {
                "q05": float(np.quantile(values, 0.05)),
                "median": float(np.median(values)),
                "q95": float(np.quantile(values, 0.95)),
            }

        call_a = scores_a > THRESHOLD
        call_b = scores_b > THRESHOLD

        correlation = float(np.corrcoef(
            scores_a, scores_b
        )[0, 1])

        if not np.isfinite(correlation):
            raise RuntimeError("Undefined joint score correlation")

        summaries[str(retention)] = {
            "retention": retention,
            "replicates": REPLICATES,
            "region_a_score": summarize(scores_a),
            "region_b_score": summarize(scores_b),
            "score_difference_b_minus_a": summarize(delta),
            "fraction_region_a_called": float(call_a.mean()),
            "fraction_region_b_called": float(call_b.mean()),
            "fraction_both_called": float(
                np.mean(both_called)
            ),
            "fraction_b_score_exceeds_a": float(
                np.mean(delta > 0)
            ),
            "score_correlation": correlation,
            "median_reference_a": float(
                np.median(reference_a)
            ),
            "median_reference_b": float(
                np.median(reference_b)
            ),
            "median_called_windows": float(
                np.median(called_counts)
            ),
        }

    return {
        "analysis": "joint_technical_candidate_score_sensitivity",
        "source_sha256": source_hash,
        "profile_sha256": profile_hash,
        "seed": SEED,
        "replicates_per_retention": REPLICATES,
        "threshold": THRESHOLD,
        "original_called_windows": 5,
        "original_merged_regions": 2,
        "baseline_scores": baseline_scores,
        "summaries": summaries,
        "interpretation": (
            "Joint technical sensitivity conditional on pooled "
            "observed counts. Reference sets are recomputed "
            "by the original detector in every thinning draw. "
            "This is not a p-value, post-selection confidence "
            "interval, or biological replication experiment."
        ),
    }


@app.local_entrypoint()
def main():
    result = evaluate_joint_scores.remote()

    output = Path(
        "research_artifacts/stage1/"
        "hct116_joint_score_thinning_v1.json"
    )

    if output.exists():
        raise RuntimeError(
            "Output already exists; refusing to overwrite"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2) + "\n"
    )

    print("\nHCT116 JOINT CANDIDATE-SCORE SENSITIVITY")
    print("=" * 65)
    print(json.dumps(result, indent=2))
    print("\nSaved:", output)
    print("Original research data unchanged.")
