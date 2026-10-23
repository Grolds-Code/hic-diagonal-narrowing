"""Read-only, window-level HCT116 W80 sensitivity diagnostics."""

import csv
import hashlib
import io
import json
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-window-reliability")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"
PROFILE = SOURCE + "_w80_profile_v2_offdiag"
DOWNSAMPLE = SOURCE + "_w80_downsampling_v1"


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=300,
)
def analyze_window_reliability():
    import numpy as np

    root = Path("/research")

    def verified_bytes(data_path, meta_path, hash_field):
        content = data_path.read_bytes()
        metadata = json.loads(meta_path.read_text())
        digest = hashlib.sha256(content).hexdigest()

        if digest != metadata[hash_field]:
            raise RuntimeError(
                f"Checksum mismatch: {data_path.name}"
            )

        return content, metadata, digest

    source_bytes, source_meta, source_hash = verified_bytes(
        root / "matrices" / f"{SOURCE}.npz",
        root / "metadata" / f"{SOURCE}.json",
        "sha256",
    )

    profile_bytes, profile_meta, profile_hash = verified_bytes(
        root / "results" / f"{PROFILE}.csv",
        root / "metadata" / f"{PROFILE}.json",
        "output_sha256",
    )

    down_bytes, down_meta, _ = verified_bytes(
        root / "results" / f"{DOWNSAMPLE}.csv",
        root / "metadata" / f"{DOWNSAMPLE}.json",
        "output_sha256",
    )

    if (
        down_meta["source_sha256"] != source_hash
        or down_meta["reference_sha256"] != profile_hash
    ):
        raise RuntimeError("Experiment provenance mismatch")

    with np.load(io.BytesIO(source_bytes), allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    profile = list(csv.DictReader(
        io.StringIO(profile_bytes.decode())
    ))

    downsampled = list(csv.DictReader(
        io.StringIO(down_bytes.decode())
    ))

    if len(profile) != 122 or len(downsampled) != 610:
        raise RuntimeError("Unexpected experiment dimensions")

    midpoint_twice = x + y
    separation = (y - x) // 500

    features = {}

    for row in profile:
        start = int(row["start_bp"])
        end = int(row["end_bp"])
        key = (start, end)

        if key in features:
            raise RuntimeError("Duplicate profile window")

        eligible = (
            (midpoint_twice >= 2 * start)
            & (midpoint_twice < 2 * end)
            & (separation >= 1)
            & (separation <= 64)
        )

        near = eligible & (separation <= 4)

        signal = float(counts[eligible].sum())

        if signal <= 0 or not np.isclose(
            signal, float(row["total_contacts"])
        ):
            raise RuntimeError("Baseline signal mismatch")

        features[key] = {
            "signal": signal,
            "occupancy": float(row["pair_occupancy"]),
            "w80": int(row["w80_bins"]),
            "near_fraction": float(
                counts[near].sum() / signal
            ),
        }

    indexed = {}

    for row in downsampled:
        key = (
            int(row["start_bp"]),
            int(row["end_bp"]),
            float(row["retention"]),
        )

        if key in indexed:
            raise RuntimeError("Duplicate downsampling row")

        window_key = key[:2]

        if window_key not in features:
            raise RuntimeError("Unknown genomic window")

        baseline = features[window_key]

        if (
            int(row["baseline_w80_bins"]) != baseline["w80"]
            or not np.isclose(
                float(row["baseline_contact_signal"]),
                baseline["signal"],
            )
        ):
            raise RuntimeError("Downsampling baseline mismatch")

        if (
            int(row["replicates"]) != 100
            or int(row["valid_replicates"]) != 100
        ):
            raise RuntimeError("Unexpected replicate coverage")

        indexed[key] = row

    def average_ranks(values):
        values = np.asarray(values, dtype=float)
        unique, inverse, sizes = np.unique(
            values,
            return_inverse=True,
            return_counts=True,
        )
        cumulative = np.cumsum(sizes)
        ranks = (2 * cumulative - sizes + 1) / 2
        return ranks[inverse]

    def spearman(a, b):
        ra = average_ranks(a)
        rb = average_ranks(b)

        if np.std(ra) == 0 or np.std(rb) == 0:
            return None

        return float(np.corrcoef(ra, rb)[0, 1])

    reports = []

    for retention in (0.25, 0.10):
        windows = []

        for key, baseline in features.items():
            row = indexed.get((*key, retention))

            if row is None:
                raise RuntimeError("Missing downsampling result")

            windows.append({
                "start_bp": key[0],
                "end_bp": key[1],
                **baseline,
                "narrowing_risk": float(
                    row["fraction_narrower_by_two_or_more_bins"]
                ),
                "exact_fraction": float(
                    row["fraction_exact"]
                ),
                "q05": float(row["thinned_w80_q05_bins"]),
                "q95": float(row["thinned_w80_q95_bins"]),
            })

        risk = np.array([
            row["narrowing_risk"] for row in windows
        ])

        associations = {}

        for feature in (
            "signal", "occupancy", "w80", "near_fraction"
        ):
            associations[feature] = spearman(
                [row[feature] for row in windows],
                risk,
            )

        by_signal = sorted(
            windows,
            key=lambda row: row["signal"],
        )

        low = by_signal[:30]
        high = by_signal[-30:]

        def group_summary(group):
            return {
                "median_signal": float(np.median(
                    [r["signal"] for r in group]
                )),
                "mean_narrowing_risk": float(np.mean(
                    [r["narrowing_risk"] for r in group]
                )),
                "mean_exact_fraction": float(np.mean(
                    [r["exact_fraction"] for r in group]
                )),
            }

        highest_risk = sorted(
            windows,
            key=lambda row: (
                -row["narrowing_risk"],
                row["start_bp"],
            ),
        )[:5]

        reports.append({
            "retention": retention,
            "windows": len(windows),
            "risk_range": [
                float(risk.min()),
                float(risk.max()),
            ],
            "risk_correlations": associations,
            "lowest_signal_30": group_summary(low),
            "highest_signal_30": group_summary(high),
            "highest_risk_windows": highest_risk,
        })

    return {
        "checksums": "PASS",
        "provenance": "PASS",
        "windows_verified": len(features),
        "downsampling_rows_verified": len(indexed),
        "reports": reports,
    }


@app.local_entrypoint()
def main():
    result = analyze_window_reliability.remote()

    print("\nHCT116 WINDOW-LEVEL W80 RELIABILITY")
    print("=" * 65)
    print("Checksums:", result["checksums"])
    print("Provenance:", result["provenance"])
    print("Windows:", result["windows_verified"])
    print(
        "Downsampling rows:",
        result["downsampling_rows_verified"],
    )

    for report in result["reports"]:
        print(
            f"\nCONTACT RETENTION: {report['retention']:.0%}"
        )
        print("Narrowing-risk range:", report["risk_range"])
        print("Spearman associations with narrowing risk:")

        for feature, correlation in report[
            "risk_correlations"
        ].items():
            print(f"  {feature}: {correlation:.4f}")

        print(
            "Lowest-signal 30:",
            report["lowest_signal_30"],
        )
        print(
            "Highest-signal 30:",
            report["highest_signal_30"],
        )

        print("Five most vulnerable windows:")

        for row in report["highest_risk_windows"]:
            print(
                f"  chr10:{row['start_bp']}-{row['end_bp']}"
                f" | W80={row['w80']} bins"
                f" | signal={row['signal']:.0f}"
                f" | occupancy={row['occupancy']:.3f}"
                f" | near_fraction={row['near_fraction']:.3f}"
                f" | narrowing_risk={row['narrowing_risk']:.1%}"
            )

    print("\nWindow-level sensitivity diagnostics: PASS")
