"""Read-only, unmasked HCT116 pair-heterogeneity diagnostic.

Fano ratios across genomic pairs are descriptive, not a
formal test of Poisson sampling or biological overdispersion.
"""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-heterogeneity-pilot")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"

START = 20_000_000
END = 21_000_000
RES = 500
WIDTH = 64
STEP = 16
DMAX = 64


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def analyze():
    import numpy as np

    root = Path("/research")

    source_bytes = (
        root / "matrices" / f"{SOURCE}.npz"
    ).read_bytes()

    source_meta = json.loads(
        (root / "metadata" / f"{SOURCE}.json").read_text()
    )

    profile_bytes = (
        root / "results" / f"{PROFILE}.csv"
    ).read_bytes()

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

    if (
        int(source_meta["start_bp"]) != START
        or int(source_meta["end_bp_exclusive"]) != END
        or int(source_meta["resolution_bp"]) != RES
        or int(source_meta["max_separation_bp"]) != DMAX * RES
    ):
        raise RuntimeError("Unexpected source geometry")

    with np.load(
        io.BytesIO(source_bytes),
        allow_pickle=False,
    ) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    if not (
        x.shape == y.shape == counts.shape
        and len(counts) == source_meta["selected_records"]
    ):
        raise RuntimeError("Source array mismatch")

    if (
        not np.all(np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
        raise RuntimeError("Counts are not valid integers")

    lookup = {}

    for a, b, count in zip(x, y, counts):
        a = int(a)
        b = int(b)
        count = int(count)

        if (
            a % RES or b % RES or a > b
            or b - a > DMAX * RES
            or not 2 * START <= a + b < 2 * END
        ):
            raise RuntimeError("Invalid contact geometry")

        if a == b:
            continue

        key = ((a - START) // RES, (b - START) // RES)

        if key in lookup:
            raise RuntimeError("Duplicate off-diagonal pair")

        lookup[key] = count

    profile_rows = list(
        csv.DictReader(io.StringIO(profile_bytes.decode()))
    )

    if len(profile_rows) != 122:
        raise RuntimeError("Unexpected profile length")

    profile_lookup = {
        (int(r["start_bp"]), int(r["end_bp"])): r
        for r in profile_rows
    }

    fano_values = []
    zero_bands = 0
    total_possible_pairs = 0
    total_positive_pairs = 0
    confirmed = 0
    regions = {}

    peaks = {
        "Region A": (20_296_000, 20_328_000),
        "Region B": (20_736_000, 20_768_000),
    }

    for s in range(0, (END - START) // RES - WIDTH + 1, STEP):
        left = START + s * RES
        right = left + WIDTH * RES
        reference = profile_lookup[(left, right)]

        band_totals = []
        band_fanos = []
        positive_pairs = 0

        for d in range(1, DMAX + 1):
            # Exactly WIDTH pairs have midpoints in [s, s+WIDTH).
            first_i = s - d // 2

            values = np.fromiter(
                (
                    lookup.get((i, i + d), 0)
                    for i in range(first_i, first_i + WIDTH)
                ),
                dtype=np.int64,
                count=WIDTH,
            )

            band_total = int(values.sum())
            band_totals.append(band_total)

            positives = int(np.count_nonzero(values))
            positive_pairs += positives

            mean = float(values.mean())

            if mean > 0:
                variance = float(values.var(ddof=1))
                fano = variance / mean
                band_fanos.append(fano)
                fano_values.append(fano)
            else:
                zero_bands += 1

        means = np.asarray(band_totals, dtype=float) / WIDTH
        total_signal = float(sum(band_totals))

        if total_signal <= 0:
            raise RuntimeError("Uninformative W80 window")

        w80 = int(
            np.searchsorted(
                np.cumsum(means),
                0.80 * means.sum(),
            ) + 1
        )

        if (
            w80 != int(reference["w80_bins"])
            or not np.isclose(
                total_signal,
                float(reference["total_contacts"]),
                rtol=0,
                atol=1e-8,
            )
        ):
            raise RuntimeError(
                f"Existing W80 profile mismatch: {left}-{right}"
            )

        confirmed += 1
        total_possible_pairs += WIDTH * DMAX
        total_positive_pairs += positive_pairs

        for name, peak in peaks.items():
            if (left, right) == peak:
                regions[name] = {
                    "interval_bp": [left, right],
                    "w80_bins": w80,
                    "positive_pair_fraction": (
                        positive_pairs / (WIDTH * DMAX)
                    ),
                    "bands_with_positive_signal": len(band_fanos),
                    "median_pair_fano": (
                        float(np.median(band_fanos))
                        if band_fanos else None
                    ),
                    "p90_pair_fano": (
                        float(np.quantile(band_fanos, 0.9))
                        if band_fanos else None
                    ),
                }

    if confirmed != 122 or len(regions) != 2:
        raise RuntimeError("Incomplete window validation")

    f = np.asarray(fano_values, dtype=float)

    return {
        "analysis": "descriptive_unmasked_pair_heterogeneity",
        "source_sha256": source_meta["sha256"],
        "profile_sha256": profile_meta["output_sha256"],
        "source_checksum_verified": True,
        "profile_checksum_verified": True,
        "w80_windows_reproduced": confirmed,
        "distance_bands_total": confirmed * DMAX,
        "distance_bands_with_positive_signal": int(f.size),
        "distance_bands_with_zero_signal": zero_bands,
        "candidate_pair_count": total_possible_pairs,
        "positive_pair_count": total_positive_pairs,
        "positive_pair_fraction": (
            total_positive_pairs / total_possible_pairs
        ),
        "median_pair_fano": float(np.median(f)),
        "p90_pair_fano": float(np.quantile(f, 0.9)),
        "fraction_fano_above_2": float(np.mean(f > 2)),
        "fraction_fano_above_5": float(np.mean(f > 5)),
        "regions": regions,
        "limitations": [
            "Unmasked pairs include assumed implicit zeros",
            "Overlapping windows reuse genomic contact pairs",
            "Genomic pair heterogeneity can inflate Fano ratios",
            "A single pooled matrix cannot estimate biological replicate variance",
            "This is not a formal Poisson goodness-of-fit test",
        ],
    }


@app.local_entrypoint()
def main():
    result = analyze.remote()

    destination = Path(
        "research_artifacts/stage1/"
        "hct116_pair_heterogeneity_pilot.json"
    )

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(result, indent=2) + "\n"
    )

    print("\nHCT116 UNMASKED PAIR-HETEROGENEITY PILOT")
    print("=" * 65)
    print(json.dumps(result, indent=2))
    print("\nSaved:", destination)
    print("Source data unchanged.")
