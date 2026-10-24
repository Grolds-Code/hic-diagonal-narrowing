"""Read-only mask-aware W80 sensitivity analysis."""

import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-mask-sensitivity")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = modal.Image.debian_slim(
    python_version="3.12"
).pip_install(
    "numpy==2.5.4",
    "pandas==3.1.0",
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"

CASES = [
    {
        "name": "Region A",
        "peak": (20296000, 20328000),
        "controls": [
            (20104000, 20136000),
            (20512000, 20544000),
            (20560000, 20592000),
            (20336000, 20368000),
            (20432000, 20464000),
        ],
    },
    {
        "name": "Region B",
        "peak": (20736000, 20768000),
        "controls": [
            (20576000, 20608000),
            (20624000, 20656000),
            (20032000, 20064000),
            (20664000, 20696000),
            (20208000, 20240000),
        ],
    },
]


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def evaluate_mask_sensitivity():
    import numpy as np
    import pandas as pd

    root = Path("/research")

    source_path = root / "matrices" / f"{SOURCE}.npz"
    profile_path = root / "results" / f"{PROFILE}.csv"

    source_bytes = source_path.read_bytes()
    profile_bytes = profile_path.read_bytes()

    source_meta = json.loads(
        (root / "metadata" / f"{SOURCE}.json").read_text()
    )
    profile_meta = json.loads(
        (root / "metadata" / f"{PROFILE}.json").read_text()
    )

    if (
        hashlib.sha256(source_bytes).hexdigest()
        != source_meta["sha256"]
    ):
        raise RuntimeError("Contact source checksum mismatch")

    if (
        hashlib.sha256(profile_bytes).hexdigest()
        != profile_meta["output_sha256"]
    ):
        raise RuntimeError("W80 profile checksum mismatch")

    with np.load(io.BytesIO(source_bytes), allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    profile = pd.read_csv(io.BytesIO(profile_bytes))

    if len(profile) != 122:
        raise RuntimeError("Expected 122 profile windows")

    start = int(source_meta["start_bp"])
    end = int(source_meta["end_bp_exclusive"])
    resolution = int(source_meta["resolution_bp"])

    if (
        resolution != 500
        or source_meta["max_separation_bp"] != 32000
    ):
        raise RuntimeError("Unexpected contact geometry")

    if (
        np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
        raise RuntimeError("Invalid contact counts")

    if np.any(x > y):
        raise RuntimeError("Invalid upper-triangle coordinates")

    n_bins = (end - start) // resolution
    separations = (y - x) // resolution

    selected = (
        (separations >= 1)
        & (separations <= 64)
        & (counts > 0)
    )

    contact_lookup = {}

    for a, b, value in zip(
        x[selected],
        y[selected],
        counts[selected],
    ):
        key = (int(a), int(b))

        if key in contact_lookup:
            raise RuntimeError("Duplicate sparse contact record")

        contact_lookup[key] = int(value)

    # Anchor support is measured in the same 0.5-32 kb
    # band as the original extracted sparse contacts.
    anchor_signal = np.zeros(
        n_bins,
        dtype=np.int64,
    )

    for endpoint in (x[selected], y[selected]):
        inside = (
            (endpoint >= start)
            & (endpoint < end)
        )

        ids = (
            (endpoint[inside] - start) // resolution
        ).astype(int)

        np.add.at(
            anchor_signal,
            ids,
            counts[selected][inside].astype(np.int64),
        )

    profile_lookup = {
        (int(row.start_bp), int(row.end_bp)): row
        for row in profile.itertuples()
    }

    def evaluate_window(interval, masked):
        left, right = interval

        if interval not in profile_lookup:
            raise RuntimeError("Window absent from profile")

        s = (left - start) // resolution
        e = (right - start) // resolution

        if e - s != 64:
            raise RuntimeError("Expected 64-bin window")

        # At each distance, select exactly the genomic
        # pairs whose midpoint belongs to this window.
        distance_means = []
        valid_pair_counts = []
        observed_contact_total = 0

        for separation in range(1, 65):
            possible_x = np.arange(
                s - separation,
                e + 1,
            )

            midpoint_inside = (
                (2 * possible_x + separation >= 2 * s)
                & (2 * possible_x + separation < 2 * e)
            )

            possible_x = possible_x[midpoint_inside]
            possible_y = possible_x + separation

            if len(possible_x) != 64:
                raise RuntimeError(
                    "Unexpected midpoint pair count"
                )

            if (
                np.any(possible_x < 0)
                or np.any(possible_y >= n_bins)
            ):
                raise RuntimeError(
                    "Contact endpoints outside extracted region"
                )

            valid = (
                ~masked[possible_x]
                & ~masked[possible_y]
            )

            valid_x = possible_x[valid]
            valid_y = possible_y[valid]

            n_valid = len(valid_x)

            if n_valid == 0:
                raise RuntimeError(
                    "No valid pairs at a separation"
                )

            distance_signal = sum(
                contact_lookup.get(
                    (
                        int(start + resolution * a),
                        int(start + resolution * b),
                    ),
                    0,
                )
                for a, b in zip(valid_x, valid_y)
            )

            observed_contact_total += distance_signal

            distance_means.append(
                distance_signal / n_valid
            )
            valid_pair_counts.append(n_valid)

        means = np.asarray(
            distance_means,
            dtype=float,
        )

        if means.sum() <= 0:
            raise RuntimeError("No usable contact signal")

        w80 = int(
            np.searchsorted(
                np.cumsum(means),
                0.8 * means.sum(),
            ) + 1
        )

        if not np.isfinite(w80):
            raise RuntimeError("Invalid masked W80")

        # Independently verify the original unmasked
        # estimator before interpreting mask effects.
        if not np.any(masked):
            original = profile_lookup[interval]

            if (
                w80 != int(original.w80_bins)
                or not np.isclose(
                    observed_contact_total,
                    original.total_contacts,
                )
            ):
                raise RuntimeError(
                    f"Baseline reconstruction failed: {interval}"
                )

        return {
            "w80_bins": w80,
            "valid_pairs": int(
                sum(valid_pair_counts)
            ),
            "observed_contact_total": int(
                observed_contact_total
            ),
            "minimum_valid_pairs_per_distance": int(
                min(valid_pair_counts)
            ),
            "masked_anchor_bins": int(
                masked[s:e].sum()
            ),
            "masked_bin_positions": [
                int(start + resolution * index)
                for index in range(s, e)
                if masked[index]
            ],
        }

    conditions = [
        ("unmasked", None),
        ("zero_signal_bins", 0),
        ("anchor_signal_at_most_3", 3),
    ]

    results = []

    for label, threshold in conditions:
        if threshold is None:
            masked = np.zeros(
                n_bins,
                dtype=bool,
            )
        else:
            masked = anchor_signal <= threshold

        case_results = []

        for case in CASES:
            peak = evaluate_window(
                case["peak"],
                masked,
            )

            controls = [
                evaluate_window(interval, masked)
                for interval in case["controls"]
            ]

            reference_w80 = float(
                np.median([
                    item["w80_bins"]
                    for item in controls
                ])
            )

            score = float(
                np.log2(
                    reference_w80 / peak["w80_bins"]
                )
            )

            case_results.append({
                "region": case["name"],
                "peak": list(case["peak"]),
                "peak_w80_bins": peak["w80_bins"],
                "control_w80_bins": [
                    item["w80_bins"]
                    for item in controls
                ],
                "control_median_w80_bins": reference_w80,
                "matched_narrowing_score": score,
                "peak_valid_pairs": peak["valid_pairs"],
                "peak_removed_pair_fraction": float(
                    1 - peak["valid_pairs"] / 4096
                ),
                "peak_masked_anchor_bins": (
                    peak["masked_anchor_bins"]
                ),
                "peak_masked_bin_positions": (
                    peak["masked_bin_positions"]
                ),
                "minimum_valid_pairs_per_distance": (
                    peak[
                        "minimum_valid_pairs_per_distance"
                    ]
                ),
                "control_masked_anchor_bins": [
                    item["masked_anchor_bins"]
                    for item in controls
                ],
            })

        results.append({
            "condition": label,
            "genome_bins_masked": int(masked.sum()),
            "comparisons": case_results,
        })

    return {
        "source_checksum": "PASS",
        "profile_checksum": "PASS",
        "baseline_reconstruction": "PASS",
        "distance_normalization": (
            "MEAN_OVER_REMAINING_VALID_PAIRS"
        ),
        "results": results,
        "independent_missing_bin_annotation": "NOT_AVAILABLE",
    }


@app.local_entrypoint()
def main():
    result = evaluate_mask_sensitivity.remote()

    print("\nHCT116 MASK-AWARE W80 SENSITIVITY")
    print("=" * 65)

    print("Source checksum:", result["source_checksum"])
    print("Profile checksum:", result["profile_checksum"])
    print(
        "Baseline reconstruction:",
        result["baseline_reconstruction"],
    )
    print(
        "Distance normalization:",
        result["distance_normalization"],
    )

    for scenario in result["results"]:
        print(
            f"\nMASK CONDITION: {scenario['condition']}"
        )
        print(
            "Genome bins masked:",
            scenario["genome_bins_masked"],
        )

        for case in scenario["comparisons"]:
            print(f"\n  {case['region']}")

            for key, value in case.items():
                if key not in ("region", "peak"):
                    print(f"    {key}: {value}")

    print(
        "\nIndependent missing-bin annotation:",
        result["independent_missing_bin_annotation"],
    )
    print("Read-only mask sensitivity: PASS")
