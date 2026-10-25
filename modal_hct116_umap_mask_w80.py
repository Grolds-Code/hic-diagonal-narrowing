"""Independent Umap-derived mask sensitivity for HCT116 W80."""

import hashlib
import io
import json
import math
import subprocess

from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-umap-mask-w80")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("ca-certificates", "curl")
    .pip_install("numpy==2.5.4", "pandas==3.1.0")
    .run_commands(
        "curl -fL --retry 3 --max-time 120 "
        "-o /usr/local/bin/bigWigToBedGraph "
        "https://hgdownload.soe.ucsc.edu/admin/exe/"
        "linux.x86_64/bigWigToBedGraph",
        "chmod +x /usr/local/bin/bigWigToBedGraph",
    )
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"
OUTPUT = SOURCE + "_umap_mask_w80_v1"

START = 20_000_000
END = 21_000_000
RESOLUTION = 500

ROOT = (
    "https://hgdownload.soe.ucsc.edu/gbdb/"
    "hg38/hoffmanMappability/"
)

TRACKS = {
    "k50": ROOT + "k50.Umap.MultiTrackMappability.bw",
    "k100": ROOT + "k100.Umap.MultiTrackMappability.bw",
}

THRESHOLDS = (0.50, 0.80, 0.95)

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
    cpu=2,
    memory=4096,
    timeout=900,
)
def evaluate_umap_masked_w80():
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

    source_hash = hashlib.sha256(source_bytes).hexdigest()
    profile_hash = hashlib.sha256(profile_bytes).hexdigest()

    if source_hash != source_meta["sha256"]:
        raise RuntimeError("Source checksum mismatch")

    if profile_hash != profile_meta["output_sha256"]:
        raise RuntimeError("Profile checksum mismatch")

    if (
        int(source_meta["start_bp"]) != START
        or int(source_meta["end_bp_exclusive"]) != END
        or int(source_meta["resolution_bp"]) != RESOLUTION
        or int(source_meta["max_separation_bp"]) != 32000
    ):
        raise RuntimeError("Unexpected genomic geometry")

    with np.load(io.BytesIO(source_bytes), allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    if (
        np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
        or np.any(x > y)
    ):
        raise RuntimeError("Invalid sparse contact records")

    profile = pd.read_csv(io.BytesIO(profile_bytes))

    if len(profile) != 122:
        raise RuntimeError("Expected 122 W80 windows")

    lookup = {}

    for a, b, value in zip(x, y, counts):
        a, b = int(a), int(b)

        if 0 < b - a <= 32000 and value > 0:
            if (a, b) in lookup:
                raise RuntimeError("Duplicate sparse contact")

            lookup[(a, b)] = int(value)

    profile_lookup = {
        (int(r.start_bp), int(r.end_bp)): r
        for r in profile.itertuples()
    }

    n_bins = (END - START) // RESOLUTION

    all_intervals = sorted(set(
        interval
        for case in CASES
        for interval in [case["peak"]] + case["controls"]
    ))

    # Cache the 64-by-64 midpoint-assigned pair geometry
    # and the original count of every possible pair.
    windows = {}

    for interval in all_intervals:
        left, right = interval

        if interval not in profile_lookup:
            raise RuntimeError("Window absent from W80 profile")

        s = (left - START) // RESOLUTION
        e = (right - START) // RESOLUTION

        if e - s != 64:
            raise RuntimeError("Unexpected window length")

        rows = []
        columns = []
        signals = []

        for distance in range(1, 65):
            possible_x = np.arange(s - distance, e + 1)

            selected = (
                (2 * possible_x + distance >= 2 * s)
                & (2 * possible_x + distance < 2 * e)
            )

            possible_x = possible_x[selected]
            possible_y = possible_x + distance

            if (
                len(possible_x) != 64
                or np.any(possible_x < 0)
                or np.any(possible_y >= n_bins)
            ):
                raise RuntimeError("Invalid midpoint pair geometry")

            observed = np.fromiter(
                (
                    lookup.get(
                        (
                            START + RESOLUTION * int(a),
                            START + RESOLUTION * int(b),
                        ),
                        0,
                    )
                    for a, b in zip(possible_x, possible_y)
                ),
                dtype=np.int64,
                count=64,
            )

            rows.append(possible_x)
            columns.append(possible_y)
            signals.append(observed)

        windows[interval] = {
            "rows": np.stack(rows),
            "columns": np.stack(columns),
            "signals": np.stack(signals),
        }

    def compute_width(interval, invalid_bins):
        window = windows[interval]

        a = window["rows"]
        b = window["columns"]
        values = window["signals"]

        valid = (
            ~invalid_bins[a]
            & ~invalid_bins[b]
        )

        denominators = valid.sum(axis=1)

        if np.any(denominators == 0):
            raise RuntimeError(
                f"Unscorable distance: {interval}"
            )

        numerators = (
            values * valid
        ).sum(axis=1)

        means = numerators / denominators
        total = float(means.sum())

        if total <= 0:
            raise RuntimeError("Uninformative masked window")

        width = int(
            np.searchsorted(
                np.cumsum(means),
                0.8 * total,
            ) + 1
        )

        s = (interval[0] - START) // RESOLUTION
        e = (interval[1] - START) // RESOLUTION

        return {
            "w80_bins": width,
            "valid_pairs": int(valid.sum()),
            "removed_pair_fraction": float(
                1 - valid.sum() / 4096
            ),
            "minimum_valid_pairs_per_distance": int(
                denominators.min()
            ),
            "masked_bins_within_window": int(
                invalid_bins[s:e].sum()
            ),
        }

    no_mask = np.zeros(n_bins, dtype=bool)

    for interval in all_intervals:
        result = compute_width(interval, no_mask)
        original = profile_lookup[interval]

        if result["w80_bins"] != int(original.w80_bins):
            raise RuntimeError(
                f"Baseline W80 mismatch: {interval}"
            )

        if not np.isclose(
            windows[interval]["signals"].sum(),
            original.total_contacts,
        ):
            raise RuntimeError(
                f"Baseline signal mismatch: {interval}"
            )

    def compare_regions(invalid_bins):
        comparisons = []

        for case in CASES:
            peak = compute_width(
                case["peak"], invalid_bins
            )

            controls = [
                compute_width(interval, invalid_bins)
                for interval in case["controls"]
            ]

            reference = float(
                np.median([
                    item["w80_bins"]
                    for item in controls
                ])
            )

            comparisons.append({
                "region": case["name"],
                "peak_interval": list(case["peak"]),
                "peak_w80_bins": peak["w80_bins"],
                "control_w80_bins": [
                    item["w80_bins"]
                    for item in controls
                ],
                "reference_median_w80_bins": reference,
                "matched_narrowing_score": float(
                    np.log2(
                        reference / peak["w80_bins"]
                    )
                ),
                "peak_masked_bins": (
                    peak["masked_bins_within_window"]
                ),
                "peak_removed_pair_fraction": (
                    peak["removed_pair_fraction"]
                ),
                "peak_minimum_valid_pairs": (
                    peak["minimum_valid_pairs_per_distance"]
                ),
                "control_masked_bins": [
                    item["masked_bins_within_window"]
                    for item in controls
                ],
            })

        return comparisons

    baseline_comparison = compare_regions(no_mask)

    track_results = {}

    for name, url in TRACKS.items():
        process = subprocess.run(
            [
                "/usr/local/bin/bigWigToBedGraph",
                "-chrom=chr10",
                f"-start={START}",
                f"-end={END}",
                url,
                "stdout",
            ],
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )

        if process.returncode != 0:
            raise RuntimeError(
                f"{name} extraction failed: "
                f"{process.stderr[-500:]}"
            )

        if not process.stdout.strip():
            raise RuntimeError(f"{name}: empty Umap output")

        # Published Umap BigWigs omit zero-valued bases.
        # Reconstruct those bases as zero rather than NaN.
        base_values = np.zeros(
            END - START,
            dtype=np.float32,
        )

        coverage = np.zeros(
            END - START,
            dtype=bool,
        )

        previous_right = START
        interval_count = 0

        for line in process.stdout.splitlines():
            fields = line.split()

            if len(fields) != 4 or fields[0] != "chr10":
                raise RuntimeError("Invalid BedGraph record")

            left = int(fields[1])
            right = int(fields[2])
            value = float(fields[3])

            if (
                left >= right
                or left < previous_right
                or not math.isfinite(value)
                or not 0 <= value <= 1
            ):
                raise RuntimeError("Invalid Umap interval")

            previous_right = right

            left = max(START, left)
            right = min(END, right)

            if right <= left:
                continue

            first = left - START
            last = right - START

            base_values[first:last] = value
            coverage[first:last] = True
            interval_count += 1

        bin_means = (
            base_values.reshape(n_bins, RESOLUTION)
            .mean(axis=1)
        )

        bin_coverage = (
            coverage.reshape(n_bins, RESOLUTION)
            .mean(axis=1)
        )

        def bin_index(position):
            return (position - START) // RESOLUTION

        # Verify against the previously observed extreme bins.
        if bin_means[bin_index(20754000)] > 1e-6:
            raise RuntimeError(
                f"{name}: low-mappability sentinel mismatch"
            )

        if bin_means[bin_index(20738000)] < 0.99:
            raise RuntimeError(
                f"{name}: control sentinel mismatch"
            )

        scenarios = []

        for threshold in THRESHOLDS:
            invalid = bin_means < threshold

            scenarios.append({
                "validity_threshold": threshold,
                "masked_bins_across_1mb": int(
                    invalid.sum()
                ),
                "comparisons": compare_regions(invalid),
            })

        track_results[name] = {
            "url": url,
            "bedgraph_sha256": hashlib.sha256(
                process.stdout.encode()
            ).hexdigest(),
            "source_intervals": interval_count,
            "reported_base_fraction": float(
                coverage.mean()
            ),
            "zero_filled_mean_mappability": float(
                base_values.mean()
            ),
            "bin_coverage_consistency": float(
                bin_coverage.mean()
            ),
            "scenarios": scenarios,
        }

    result = {
        "study": "HCT116 independent Umap-mask W80 sensitivity",
        "assembly": "GRCh38",
        "chromosome": "chr10",
        "region": [START, END],
        "resolution_bp": RESOLUTION,
        "source_sha256": source_hash,
        "profile_sha256": profile_hash,
        "mask_rule": "mean_500bp_Umap_below_threshold",
        "omitted_umap_values": "interpreted_as_zero",
        "distance_normalization": (
            "mean_contact_count_over_remaining_valid_pairs"
        ),
        "baseline": baseline_comparison,
        "tracks": track_results,
        "limitations": (
            "Read lengths are hypothetical; thresholds are "
            "sensitivity settings; missing-bin identity and "
            "biological significance are not established."
        ),
    }

    output_path = (
        root / "results" / f"{OUTPUT}.json"
    )

    metadata_path = (
        root / "metadata" / f"{OUTPUT}.json"
    )

    if output_path.exists() or metadata_path.exists():
        raise RuntimeError("Versioned output already exists")

    serialized = (
        json.dumps(result, indent=2) + "\n"
    ).encode()

    digest = hashlib.sha256(serialized).hexdigest()

    metadata = {
        "output_sha256": digest,
        "source_sha256": source_hash,
        "profile_sha256": profile_hash,
        "output_file": str(output_path),
    }

    output_path.write_bytes(serialized)
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    volume.commit()

    return {
        "baseline_validation": "PASS",
        "baseline": baseline_comparison,
        "tracks": {
            name: {
                "reported_base_fraction": (
                    data["reported_base_fraction"]
                ),
                "scenarios": data["scenarios"],
            }
            for name, data in track_results.items()
        },
        "output_path": str(output_path),
        "metadata_path": str(metadata_path),
    }


@app.local_entrypoint()
def main():
    result = evaluate_umap_masked_w80.remote()

    print("\nHCT116 INDEPENDENT UMAP MASK W80 ANALYSIS")
    print("=" * 65)
    print("Baseline validation:", result["baseline_validation"])

    for item in result["baseline"]:
        print(
            f"BASELINE {item['region']}: "
            f"W80={item['peak_w80_bins']}, "
            f"reference={item['reference_median_w80_bins']}, "
            f"score={item['matched_narrowing_score']:.4f}"
        )

    for name, track in result["tracks"].items():
        print(f"\nUMAP TRACK: {name}")
        print(
            "Reported base fraction:",
            track["reported_base_fraction"],
        )

        for scenario in track["scenarios"]:
            print(
                f"\n  Threshold: "
                f"{scenario['validity_threshold']}"
            )

            print(
                "  Masked bins in 1 Mb:",
                scenario["masked_bins_across_1mb"],
            )

            for item in scenario["comparisons"]:
                print(
                    f"  {item['region']}: "
                    f"W80={item['peak_w80_bins']}, "
                    f"ref={item['reference_median_w80_bins']}, "
                    f"score={item['matched_narrowing_score']:.4f}, "
                    f"masked={item['peak_masked_bins']}, "
                    f"removed_pairs="
                    f"{item['peak_removed_pair_fraction']:.1%}"
                )

    print("\nResults:", result["output_path"])
    print("Metadata:", result["metadata_path"])
    print("\nUmap-mask sensitivity: PASS")
    print("Statistical significance: NOT ASSESSED")
