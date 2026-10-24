"""Read-only coverage-matched reference feasibility study."""

import csv
import hashlib
import io
import json
import statistics
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-matched-references")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = modal.Image.debian_slim(
    python_version="3.12"
)

PROFILE = (
    "hct116_chr10_20000000_21000000_500bp_"
    "midpoint_32kb_w80_profile_v2_offdiag"
)

CANDIDATE_FOOTPRINTS = [
    (20_296_000, 20_336_000),
    (20_720_000, 20_792_000),
]

PEAK_WINDOWS = [
    (20_296_000, 20_328_000),
    (20_736_000, 20_768_000),
]

SIGNAL_TOLERANCE = 0.30
OCCUPANCY_TOLERANCE = 0.075
MAX_CONTROLS = 5


def overlaps(a, b):
    return a[0] < b[1] and b[0] < a[1]


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=512,
    timeout=120,
)
def inspect_matched_references():
    profile_path = Path(
        f"/research/results/{PROFILE}.csv"
    )

    metadata_path = Path(
        f"/research/metadata/{PROFILE}.json"
    )

    content = profile_path.read_bytes()
    metadata = json.loads(metadata_path.read_text())

    if (
        hashlib.sha256(content).hexdigest()
        != metadata["output_sha256"]
    ):
        raise RuntimeError("Profile checksum mismatch")

    raw_rows = list(csv.DictReader(
        io.StringIO(content.decode("utf-8"))
    ))

    if len(raw_rows) != 122:
        raise RuntimeError("Unexpected profile length")

    windows = []

    for row in raw_rows:
        window = {
            "start": int(row["start_bp"]),
            "end": int(row["end_bp"]),
            "w80": int(row["w80_bins"]),
            "signal": float(row["total_contacts"]),
            "occupancy": float(row["pair_occupancy"]),
        }

        if window["signal"] <= 0:
            raise RuntimeError("Nonpositive window signal")

        windows.append(window)

    results = []

    for candidate_index, peak in enumerate(
        PEAK_WINDOWS, start=1
    ):
        targets = [
            w for w in windows
            if (w["start"], w["end"]) == peak
        ]

        if len(targets) != 1:
            raise RuntimeError(
                "Candidate peak not uniquely identified"
            )

        target = targets[0]

        eligible = [
            w for w in windows
            if all(
                not overlaps(
                    (w["start"], w["end"]),
                    footprint,
                )
                for footprint in CANDIDATE_FOOTPRINTS
            )
        ]

        matches = [
            w for w in eligible
            if (
                abs(w["signal"] / target["signal"] - 1)
                <= SIGNAL_TOLERANCE
                and abs(
                    w["occupancy"] - target["occupancy"]
                ) <= OCCUPANCY_TOLERANCE
            )
        ]

        # Rank by relative signal and occupancy differences.
        matches.sort(
            key=lambda w: (
                abs(
                    w["signal"] / target["signal"] - 1
                ) / SIGNAL_TOLERANCE
                + abs(
                    w["occupancy"] - target["occupancy"]
                ) / OCCUPANCY_TOLERANCE,
                w["start"],
            )
        )

        # Select mutually non-overlapping reference windows.
        chosen = []

        for window in matches:
            if all(
                not overlaps(
                    (window["start"], window["end"]),
                    (other["start"], other["end"]),
                )
                for other in chosen
            ):
                chosen.append(window)

            if len(chosen) == MAX_CONTROLS:
                break

        if chosen:
            reference_w80 = statistics.median(
                w["w80"] for w in chosen
            )
        else:
            reference_w80 = None

        results.append({
            "candidate": candidate_index,
            "peak": target,
            "eligible_pool": len(eligible),
            "coverage_matched_pool": len(matches),
            "nonoverlapping_controls": chosen,
            "matched_reference_w80": reference_w80,
        })

    return {
        "checksum": "PASS",
        "windows_verified": len(windows),
        "matching": {
            "signal_relative_tolerance": SIGNAL_TOLERANCE,
            "occupancy_absolute_tolerance": OCCUPANCY_TOLERANCE,
            "maximum_controls": MAX_CONTROLS,
        },
        "results": results,
    }


@app.local_entrypoint()
def main():
    result = inspect_matched_references.remote()

    print("\nHCT116 COVERAGE-MATCHED REFERENCE PROBE")
    print("=" * 65)

    print("Checksum:", result["checksum"])
    print("Windows verified:", result["windows_verified"])
    print("Matching criteria:", result["matching"])

    for item in result["results"]:
        target = item["peak"]

        print(
            f"\nREGION {item['candidate']}: "
            f"chr10:{target['start']}-{target['end']}"
        )

        print("Target W80:", target["w80"], "bins")
        print("Target signal:", target["signal"])
        print("Target occupancy:", target["occupancy"])
        print("Eligible pool:", item["eligible_pool"])
        print(
            "Coverage-matched pool:",
            item["coverage_matched_pool"],
        )
        print(
            "Nonoverlapping controls:",
            len(item["nonoverlapping_controls"]),
        )
        print(
            "Matched reference W80:",
            item["matched_reference_w80"],
        )

        for control in item["nonoverlapping_controls"]:
            print(
                f"  chr10:{control['start']}-"
                f"{control['end']}"
                f" | W80={control['w80']}"
                f" | signal={control['signal']:.0f}"
                f" | occupancy={control['occupancy']:.3f}"
            )

    print("\nReference feasibility analysis: PASS")
    print("Statistical significance: NOT ASSESSED")
