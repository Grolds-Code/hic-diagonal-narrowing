"""Read-only HCT116 W80 intervals under a conditional Poisson model.

These intervals do not establish biological 95% coverage.
"""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-conditional-intervals")

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
WINDOW = 64
STEP = 16
DMAX = 64


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def extract_verified_bands():
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

    source_sha = hashlib.sha256(source_bytes).hexdigest()
    profile_sha = hashlib.sha256(profile_bytes).hexdigest()

    if (
        source_sha != source_meta["sha256"]
        or profile_sha != profile_meta["output_sha256"]
        or profile_meta["source_sha256"] != source_sha
    ):
        raise RuntimeError("Research provenance checksum mismatch")

    if (
        int(source_meta["start_bp"]) != START
        or int(source_meta["end_bp_exclusive"]) != END
        or int(source_meta["resolution_bp"]) != RES
        or int(source_meta["max_separation_bp"]) != DMAX * RES
    ):
        raise RuntimeError("Unexpected genomic geometry")

    with np.load(
        io.BytesIO(source_bytes), allow_pickle=False
    ) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    if not (x.shape == y.shape == counts.shape):
        raise RuntimeError("Invalid source dimensions")

    if len(counts) != source_meta["selected_records"]:
        raise RuntimeError("Unexpected source record count")

    if (
        not np.all(np.isfinite(counts))
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
        raise RuntimeError("Invalid integer contact counts")

    lookup = {}

    for a, b, value in zip(x, y, counts):
        a, b, value = int(a), int(b), int(value)

        if (
            a > b
            or a % RES
            or b % RES
            or b - a > DMAX * RES
            or not 2 * START <= a + b < 2 * END
        ):
            raise RuntimeError("Invalid contact geometry")

        if a == b:
            continue

        key = ((a - START) // RES, (b - START) // RES)

        if key in lookup:
            raise RuntimeError("Duplicate off-diagonal contact")

        lookup[key] = value

    profile_rows = list(
        csv.DictReader(io.StringIO(profile_bytes.decode()))
    )

    reference = {
        (int(row["start_bp"]), int(row["end_bp"])): row
        for row in profile_rows
    }

    if len(profile_rows) != 122 or len(reference) != 122:
        raise RuntimeError("Unexpected corrected W80 profile")

    windows = []

    for s in range(
        0, (END - START) // RES - WINDOW + 1, STEP
    ):
        left = START + s * RES
        right = left + WINDOW * RES

        bands = []

        for d in range(1, DMAX + 1):
            first_i = s - d // 2

            values = [
                lookup.get((i, i + d), 0)
                for i in range(first_i, first_i + WINDOW)
            ]

            if any(v < 0 for v in values):
                raise RuntimeError("Invalid reconstructed contacts")

            if not all(
                2 * s <= 2 * i + d < 2 * (s + WINDOW)
                for i in range(first_i, first_i + WINDOW)
            ):
                raise RuntimeError("Midpoint geometry mismatch")

            bands.append(int(sum(values)))

        if sum(bands) <= 0:
            raise RuntimeError("Empty-signal genomic window")

        cumulative = np.cumsum(bands, dtype=np.int64)
        w80 = int(
            np.searchsorted(
                cumulative, 0.80 * cumulative[-1]
            ) + 1
        )

        original = reference[(left, right)]

        if (
            w80 != int(original["w80_bins"])
            or sum(bands) != int(
                float(original["total_contacts"])
            )
        ):
            raise RuntimeError(
                f"Corrected profile mismatch: {left}-{right}"
            )

        windows.append({
            "start_bp": left,
            "end_bp": right,
            "w80_bins": w80,
            "band_counts": bands,
        })

    if len(windows) != 122:
        raise RuntimeError("Incomplete window extraction")

    return {
        "source_sha256": source_sha,
        "profile_sha256": profile_sha,
        "windows": windows,
    }


@app.local_entrypoint()
def main():
    import numpy as np

    from hic_diagonal_narrowing.poisson_uncertainty import (
        cutoff_exposure_w80_interval_from_bands,
    )

    extracted = extract_verified_bands.remote()

    output = Path(
        "research_artifacts/stage1/"
        "hct116_conditional_w80_intervals_v1.csv"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    for item in extracted["windows"]:
        bands = np.asarray(
            item["band_counts"], dtype=np.int64
        )
        exposures = np.full(DMAX, WINDOW, dtype=np.int64)

        result = cutoff_exposure_w80_interval_from_bands(
            bands,
            exposures,
            confidence_level=0.95,
        )

        observed = int(result["observed_w80"])
        lower = int(result["cutoff_lower"])
        upper = int(result["cutoff_upper"])

        if observed != item["w80_bins"]:
            raise RuntimeError("W80 implementation mismatch")

        if not 1 <= lower <= observed <= upper <= DMAX:
            raise RuntimeError("Invalid uncertainty interval")

        rows.append({
            "start_bp": item["start_bp"],
            "end_bp": item["end_bp"],
            "w80_bins": observed,
            "lower_95_bins": lower,
            "upper_95_bins": upper,
            "interval_span_bins": upper - lower,
            "total_contacts": int(bands.sum()),
            "minimum_eligible_pairs": WINDOW,
            "maximum_eligible_pairs": WINDOW,
            "model": "independent_poisson_band_counts",
            "mask": "none",
        })

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    spans = np.array(
        [row["interval_span_bins"] for row in rows]
    )

    named = {
        "Region A": (20_296_000, 20_328_000),
        "Region B": (20_736_000, 20_768_000),
    }

    candidates = {}

    for name, interval in named.items():
        matches = [
            r for r in rows
            if (r["start_bp"], r["end_bp"]) == interval
        ]

        if len(matches) != 1:
            raise RuntimeError(f"Missing {name}")

        candidates[name] = matches[0]

    summary = {
        "verified_windows": len(rows),
        "source_sha256": extracted["source_sha256"],
        "profile_sha256": extracted["profile_sha256"],
        "median_interval_span_bins": float(np.median(spans)),
        "p90_interval_span_bins": float(
            np.quantile(spans, 0.9)
        ),
        "max_interval_span_bins": int(spans.max()),
        "full_range_interval_count": int(
            sum(
                r["lower_95_bins"] == 1
                and r["upper_95_bins"] == DMAX
                for r in rows
            )
        ),
        "candidates": candidates,
        "assumptions": [
            "Independent Poisson distance-band counts",
            "Fixed equal exposure of 64 pairs per distance",
            "Unreported eligible contacts treated as zero",
            "No genomic missing-bin mask applied",
        ],
        "interpretation": (
            "Model-conditional sensitivity analysis; "
            "not validated biological confidence coverage"
        ),
    }

    print("\nREAL HCT116 MODEL-CONDITIONAL W80 INTERVALS")
    print("=" * 65)
    print(json.dumps(summary, indent=2))
    print("\nSaved:", output)
    print("Source data unchanged.")
