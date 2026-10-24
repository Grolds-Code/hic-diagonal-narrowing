"""Read-only HCT116 per-bin contact-support audit."""

import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-bin-support-audit")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4", "pandas==3.1.0")
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
def audit_bin_support():
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

    profile = pd.read_csv(io.BytesIO(profile_bytes))

    if len(profile) != 122:
        raise RuntimeError("Unexpected W80 profile length")

    start = int(source_meta["start_bp"])
    end = int(source_meta["end_bp_exclusive"])
    resolution = int(source_meta["resolution_bp"])

    if (
        resolution != 500
        or source_meta["max_separation_bp"] != 32000
    ):
        raise RuntimeError("Unexpected data geometry")

    if (
        np.any(~np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
        raise RuntimeError("Invalid contact counts")

    if np.any(x > y):
        raise RuntimeError("Expected upper-triangle contacts")

    n_bins = (end - start) // resolution
    c = counts.astype(np.int64)
    distances = (y - x) // resolution
    midpoints_twice = x + y

    eligible = (
        (distances >= 1)
        & (distances <= 64)
        & (c > 0)
    )

    # Each off-diagonal record contributes to both
    # endpoint genomic bins.
    anchor_signal = np.zeros(n_bins, dtype=np.int64)
    anchor_partners = np.zeros(n_bins, dtype=np.int64)

    for endpoints in (x, y):
        inside = (
            eligible
            & (endpoints >= start)
            & (endpoints < end)
        )

        bin_ids = (
            (endpoints[inside] - start) // resolution
        ).astype(int)

        anchor_signal += np.bincount(
            bin_ids,
            weights=c[inside],
            minlength=n_bins,
        ).astype(np.int64)

        anchor_partners += np.bincount(
            bin_ids,
            minlength=n_bins,
        ).astype(np.int64)

    profile_lookup = {
        (int(row.start_bp), int(row.end_bp)): row
        for row in profile.itertuples()
    }

    def window_measurements(interval):
        left, right = interval

        if (
            left < start
            or right > end
            or (right - left) != 32000
        ):
            raise RuntimeError("Invalid audit window")

        if interval not in profile_lookup:
            raise RuntimeError("Window missing from profile")

        row = profile_lookup[interval]

        selected = (
            eligible
            & (midpoints_twice >= 2 * left)
            & (midpoints_twice < 2 * right)
        )

        d = distances[selected]
        v = c[selected]

        band_signal = np.bincount(
            d,
            weights=v,
            minlength=65,
        )

        band_pairs = np.bincount(
            d,
            minlength=65,
        )

        total = float(band_signal[1:65].sum())

        if total <= 0:
            raise RuntimeError("Empty candidate/control window")

        w80 = int(
            np.searchsorted(
                np.cumsum(band_signal[1:65]),
                0.8 * total,
            ) + 1
        )

        pair_occupancy = float(
            band_pairs[1:65].sum() / (64 * 64)
        )

        if (
            w80 != int(row.w80_bins)
            or not np.isclose(total, row.total_contacts)
            or not np.isclose(
                pair_occupancy,
                row.pair_occupancy,
                atol=1e-10,
            )
        ):
            raise RuntimeError(
                f"W80, signal or occupancy mismatch: {interval}"
            )

        first = (left - start) // resolution
        last = (right - start) // resolution

        signals = anchor_signal[first:last]
        partners = anchor_partners[first:last]

        if len(signals) != 64:
            raise RuntimeError("Incorrect anchor-bin count")

        bands = {}

        for name, lo, hi in (
            ("0.5_to_5kb", 1, 10),
            ("5.5_to_16kb", 11, 32),
            ("16.5_to_32kb", 33, 64),
        ):
            bands[name] = {
                "signal_fraction": float(
                    band_signal[lo:hi + 1].sum() / total
                ),
                "pair_occupancy": float(
                    band_pairs[lo:hi + 1].sum()
                    / (64 * (hi - lo + 1))
                ),
            }

        return {
            "interval": [left, right],
            "w80_bins": w80,
            "total_contacts": total,
            "pair_occupancy": pair_occupancy,
            "anchor_signal": signals,
            "anchor_partners": partners,
            "bands": bands,
        }

    results = []

    for case in CASES:
        peak = window_measurements(case["peak"])

        controls = [
            window_measurements(interval)
            for interval in case["controls"]
        ]

        all_reference_signals = np.concatenate([
            item["anchor_signal"]
            for item in controls
        ])

        reference_p10 = float(
            np.percentile(all_reference_signals, 10)
        )

        candidate_signals = peak["anchor_signal"]
        candidate_partners = peak["anchor_partners"]

        lowest_ids = np.argsort(
            candidate_signals,
            kind="stable",
        )[:8]

        weakest_bins = [
            {
                "start_bp": int(
                    case["peak"][0] + index * resolution
                ),
                "signal": int(candidate_signals[index]),
                "observed_partners": int(
                    candidate_partners[index]
                ),
            }
            for index in lowest_ids
        ]

        band_comparison = {}

        for band in peak["bands"]:
            band_comparison[band] = {
                "candidate_signal_fraction": (
                    peak["bands"][band]["signal_fraction"]
                ),
                "control_median_signal_fraction": float(
                    np.median([
                        item["bands"][band]["signal_fraction"]
                        for item in controls
                    ])
                ),
                "candidate_pair_occupancy": (
                    peak["bands"][band]["pair_occupancy"]
                ),
                "control_median_pair_occupancy": float(
                    np.median([
                        item["bands"][band]["pair_occupancy"]
                        for item in controls
                    ])
                ),
            }

        results.append({
            "region": case["name"],
            "peak": list(case["peak"]),
            "peak_w80_bins": peak["w80_bins"],
            "reference_median_w80_bins": float(
                np.median([
                    item["w80_bins"]
                    for item in controls
                ])
            ),
            "peak_zero_signal_bins": int(
                np.count_nonzero(candidate_signals == 0)
            ),
            "peak_zero_partner_bins": int(
                np.count_nonzero(candidate_partners == 0)
            ),
            "peak_min_anchor_signal": int(
                candidate_signals.min()
            ),
            "peak_median_anchor_signal": float(
                np.median(candidate_signals)
            ),
            "reference_median_anchor_signal": float(
                np.median(all_reference_signals)
            ),
            "reference_anchor_signal_p10": reference_p10,
            "peak_bins_below_reference_p10": int(
                np.count_nonzero(
                    candidate_signals < reference_p10
                )
            ),
            "peak_median_observed_partners": float(
                np.median(candidate_partners)
            ),
            "weakest_eight_bins": weakest_bins,
            "distance_band_comparison": band_comparison,
        })

    return {
        "source_checksum": "PASS",
        "profile_checksum": "PASS",
        "peak_profile_reconstruction": "PASS",
        "reference_windows_per_region": 5,
        "resolution_bp": resolution,
        "results": results,
        "missing_bin_status": (
            "NOT_DETERMINED_FROM_CONTACT_SUPPORT_ALONE"
        ),
    }


@app.local_entrypoint()
def main():
    result = audit_bin_support.remote()

    print("\nHCT116 GENOMIC BIN SUPPORT AUDIT")
    print("=" * 65)

    for key in (
        "source_checksum",
        "profile_checksum",
        "peak_profile_reconstruction",
        "reference_windows_per_region",
        "resolution_bp",
    ):
        print(f"{key}: {result[key]}")

    for item in result["results"]:
        print(f"\n{item['region']}: {item['peak']}")

        for key, value in item.items():
            if key not in (
                "region",
                "peak",
                "weakest_eight_bins",
                "distance_band_comparison",
            ):
                print(f"{key}: {value}")

        print("\nDistance-band comparisons:")

        for name, comparison in item[
            "distance_band_comparison"
        ].items():
            print(f"  {name}: {comparison}")

        print("\nEight least-supported genomic bins:")

        for bin_info in item["weakest_eight_bins"]:
            print(f"  {bin_info}")

    print(
        "\nMissing-bin assessment:",
        result["missing_bin_status"],
    )
    print("Read-only support audit: PASS")
