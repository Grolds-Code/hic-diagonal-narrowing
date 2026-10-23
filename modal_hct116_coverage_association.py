"""Read-only HCT116 W80 and coverage association diagnostics."""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal


app = modal.App("hic-diagonal-narrowing-coverage-association")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"

RESOLUTION = 500
MAX_SEPARATION_BINS = 64
WINDOW_BINS = 64


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def analyze_coverage_association():
    import numpy as np

    source_path = Path(f"/research/matrices/{SOURCE}.npz")
    source_meta_path = Path(f"/research/metadata/{SOURCE}.json")

    profile_path = Path(f"/research/results/{PROFILE}.csv")
    profile_meta_path = Path(f"/research/metadata/{PROFILE}.json")

    source_meta = json.loads(source_meta_path.read_text())
    profile_meta = json.loads(profile_meta_path.read_text())

    source_bytes = source_path.read_bytes()
    profile_bytes = profile_path.read_bytes()

    if hashlib.sha256(source_bytes).hexdigest() != source_meta["sha256"]:
        raise RuntimeError("Source checksum mismatch")

    if hashlib.sha256(profile_bytes).hexdigest() != profile_meta["output_sha256"]:
        raise RuntimeError("Profile checksum mismatch")

    with np.load(source_path, allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    rows = list(csv.DictReader(
        io.StringIO(profile_bytes.decode("utf-8"))
    ))

    if len(rows) != 122:
        raise RuntimeError("Expected 122 windows")

    midpoints_twice = x + y
    separations = (y - x) // RESOLUTION

    widths = []
    signals = []
    occupancies = []

    for row in rows:
        start = int(row["start_bp"])
        end = int(row["end_bp"])

        if end - start != WINDOW_BINS * RESOLUTION:
            raise RuntimeError("Unexpected window size")

        selected = (
            (midpoints_twice >= 2 * start)
            & (midpoints_twice < 2 * end)
            & (separations >= 1)
            & (separations <= MAX_SEPARATION_BINS)
        )

        local_counts = counts[selected]
        local_separations = separations[selected]

        distance_signal = np.bincount(
            local_separations,
            weights=local_counts,
            minlength=MAX_SEPARATION_BINS + 1,
        )[1:]

        total = float(distance_signal.sum())

        if total <= 0:
            raise RuntimeError("Window has no off-diagonal signal")

        recomputed_w80 = int(
            np.searchsorted(
                np.cumsum(distance_signal),
                0.8 * total,
                side="left",
            ) + 1
        )

        reported_w80 = int(row["w80_bins"])

        if recomputed_w80 != reported_w80:
            raise RuntimeError(
                f"W80 mismatch at {start}: "
                f"{recomputed_w80} != {reported_w80}"
            )

        if not np.isclose(
            total,
            float(row["total_contacts"])
        ):
            raise RuntimeError("Window contact signal mismatch")

        possible_pairs = (
            WINDOW_BINS * MAX_SEPARATION_BINS
        )

        occupancy = (
            np.count_nonzero(local_counts > 0)
            / possible_pairs
        )

        if not np.isclose(
            occupancy,
            float(row["pair_occupancy"])
        ):
            raise RuntimeError("Window occupancy mismatch")

        widths.append(reported_w80)
        signals.append(total)
        occupancies.append(occupancy)

    widths = np.asarray(widths, dtype=float)
    signals = np.asarray(signals, dtype=float)
    occupancies = np.asarray(occupancies, dtype=float)

    def average_ranks(values):
        order = np.argsort(values, kind="stable")
        sorted_values = values[order]

        starts = np.r_[
            0,
            np.flatnonzero(
                sorted_values[1:] != sorted_values[:-1]
            ) + 1,
        ]

        ends = np.r_[starts[1:], len(values)]
        ranks = np.empty(len(values), dtype=float)

        for first, last in zip(starts, ends):
            ranks[order[first:last]] = (
                first + last + 1
            ) / 2

        return ranks

    def spearman(a, b):
        ranks_a = average_ranks(a)
        ranks_b = average_ranks(b)

        if (
            np.std(ranks_a) == 0
            or np.std(ranks_b) == 0
        ):
            return None

        return float(
            np.corrcoef(ranks_a, ranks_b)[0, 1]
        )

    quartile_size = len(rows) // 4
    sorted_indices = np.argsort(signals)

    low_indices = sorted_indices[:quartile_size]
    high_indices = sorted_indices[-quartile_size:]

    def summarize_group(indices):
        return {
            "windows": len(indices),
            "median_contact_signal": float(
                np.median(signals[indices])
            ),
            "median_w80_bins": float(
                np.median(widths[indices])
            ),
            "median_pair_occupancy": float(
                np.median(occupancies[indices])
            ),
        }

    return {
        "source_checksum": "PASS",
        "profile_checksum": "PASS",
        "w80_recomputation": "PASS",
        "windows_verified": len(rows),
        "w80_range_bins": [
            int(widths.min()),
            int(widths.max()),
        ],
        "contact_signal_range": [
            float(signals.min()),
            float(signals.max()),
        ],
        "occupancy_range": [
            float(occupancies.min()),
            float(occupancies.max()),
        ],
        "spearman_w80_vs_signal": spearman(
            widths, signals
        ),
        "spearman_w80_vs_occupancy": spearman(
            widths, occupancies
        ),
        "spearman_signal_vs_occupancy": spearman(
            signals, occupancies
        ),
        "lowest_signal_quartile": summarize_group(
            low_indices
        ),
        "highest_signal_quartile": summarize_group(
            high_indices
        ),
        "missing_bin_status": source_meta[
            "missing_bin_status"
        ],
        "statistical_significance": "not_assessed",
    }


@app.local_entrypoint()
def main():
    result = analyze_coverage_association.remote()

    print("\nHCT116 W80–COVERAGE ASSOCIATION")
    print("=" * 55)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nRead-only association diagnostics: PASS")
