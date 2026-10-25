"""Read-only HCT116 bin visibility across 0.5 kb to 1 Mb."""

import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-long-range-visibility")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install(
        "libcurl4-openssl-dev",
        "zlib1g-dev",
    )
    .pip_install(
        "hic-straw==1.3.1",
        "numpy==2.5.4",
    )
)

SOURCE = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"

LOW_SUPPORT_BINS = [
    20752000,
    20753000,
    20754000,
    20755000,
    20755500,
    20756000,
    20756500,
    20757000,
]

CHILD_CODE = r'''
import json
import math
import sys

import hicstraw

position = int(sys.argv[1])

url = (
    "https://huggingface.co/buckets/FatimaInstitute/"
    "3d-genomics-prediction/resolve/HCT116/"
    "ENCFF573OPJ.hic?download=true"
)

bands = [
    ("0.5-32kb", 500, 32000),
    ("32.5-100kb", 32500, 100000),
    ("100.5-250kb", 100500, 250000),
    ("250.5-1000kb", 250500, 1000000),
]

result = {
    "position_bp": position,
    "status": "QUERY_COMPLETED",
    "bands": {
        name: {"pairs": 0, "signal": 0.0}
        for name, _, _ in bands
    },
}

try:
    hic = hicstraw.HiCFile(url)

    matrix = hic.getMatrixZoomData(
        "chr10",
        "chr10",
        "observed",
        "NONE",
        "BP",
        500,
    )

    # Query the upstream and downstream sides separately.
    # This avoids constructing a large square Hi-C matrix.
    upstream = matrix.getRecords(
        position - 1000000,
        position - 1,
        position,
        position + 499,
    )

    downstream = matrix.getRecords(
        position,
        position + 499,
        position + 500,
        position + 1000000 + 499,
    )

    observed_pairs = set()

    for record in list(upstream) + list(downstream):
        x = int(record.binX)
        y = int(record.binY)
        value = float(record.counts)

        if position not in (x, y):
            raise ValueError(
                f"Unexpected pair outside anchor: {x}, {y}"
            )

        other = y if x == position else x
        separation = abs(other - position)

        if (
            separation < 500
            or separation > 1000000
            or separation % 500 != 0
        ):
            raise ValueError("Invalid contact separation")

        if not math.isfinite(value) or value < 0:
            raise ValueError("Nonfinite or negative contact")

        pair = (min(x, y), max(x, y))

        if pair in observed_pairs:
            raise ValueError("Duplicate pair in query")

        observed_pairs.add(pair)

        for name, lower, upper in bands:
            if lower <= separation <= upper:
                result["bands"][name]["pairs"] += 1
                result["bands"][name]["signal"] += value
                break

    result["total_pairs"] = len(observed_pairs)
    result["total_signal"] = sum(
        item["signal"]
        for item in result["bands"].values()
    )

except Exception as exc:
    result = {
        "position_bp": position,
        "status": "QUERY_ERROR",
        "error_type": type(exc).__name__,
        "error": str(exc)[:500],
    }

print(json.dumps(result, sort_keys=True))
'''


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=1500,
)
def audit_long_range_visibility():
    import numpy as np

    root = Path("/research")
    source_path = root / "matrices" / f"{SOURCE}.npz"

    source_bytes = source_path.read_bytes()
    metadata = json.loads(
        (root / "metadata" / f"{SOURCE}.json").read_text()
    )

    if (
        hashlib.sha256(source_bytes).hexdigest()
        != metadata["sha256"]
    ):
        raise RuntimeError("Source checksum mismatch")

    with np.load(io.BytesIO(source_bytes), allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    # Select three observed high-support comparison bins
    # from separated portions of Region B's peak window.
    peak_start = 20736000
    peak_positions = np.arange(
        peak_start,
        peak_start + 32000,
        500,
    )

    local_signal = []
    local_pairs = []

    def expected_local_support(position):
        selected = (
            (y > x)
            & (y - x <= 32000)
            & ((x == position) | (y == position))
        )

        return {
            "pairs": int(np.count_nonzero(selected)),
            "signal": float(counts[selected].sum()),
        }

    for position in peak_positions:
        expected = expected_local_support(int(position))
        local_signal.append(expected["signal"])
        local_pairs.append(expected["pairs"])

    local_signal = np.asarray(local_signal)

    controls = []

    for lower, upper in ((0, 16), (16, 32), (48, 64)):
        selected_index = (
            lower
            + int(np.argmax(local_signal[lower:upper]))
        )

        controls.append(int(peak_positions[selected_index]))

    if any(
        expected_local_support(position)["signal"] <= 0
        for position in controls
    ):
        raise RuntimeError("Invalid comparison-bin selection")

    positions = LOW_SUPPORT_BINS + controls
    results = []

    for position in positions:
        try:
            completed = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    CHILD_CODE,
                    str(position),
                ],
                capture_output=True,
                text=True,
                timeout=100,
                check=False,
            )

            if completed.returncode != 0:
                outcome = {
                    "position_bp": position,
                    "status": "CHILD_PROCESS_FAILED",
                    "returncode": completed.returncode,
                    "stderr_tail": completed.stderr[-500:],
                }
            else:
                lines = [
                    line.strip()
                    for line in completed.stdout.splitlines()
                    if line.strip()
                ]

                try:
                    outcome = json.loads(lines[-1])
                except (IndexError, json.JSONDecodeError):
                    outcome = {
                        "position_bp": position,
                        "status": "UNPARSEABLE_OUTPUT",
                        "stdout_tail": completed.stdout[-500:],
                    }

        except subprocess.TimeoutExpired:
            outcome = {
                "position_bp": position,
                "status": "TIMEOUT",
            }

        expected = expected_local_support(position)
        outcome["source_32kb_expected"] = expected

        if outcome["status"] == "QUERY_COMPLETED":
            observed_local = outcome["bands"]["0.5-32kb"]

            outcome["source_32kb_agreement"] = bool(
                observed_local["pairs"] == expected["pairs"]
                and np.isclose(
                    observed_local["signal"],
                    expected["signal"],
                )
            )
        else:
            outcome["source_32kb_agreement"] = None

        outcome["group"] = (
            "low_support"
            if position in LOW_SUPPORT_BINS
            else "comparison"
        )

        results.append(outcome)

    all_success = all(
        result["status"] == "QUERY_COMPLETED"
        and result["source_32kb_agreement"] is True
        for result in results
    )

    return {
        "source_checksum": "PASS",
        "resolution_bp": 500,
        "normalization": "NONE",
        "maximum_separation_bp": 1000000,
        "low_support_positions": LOW_SUPPORT_BINS,
        "comparison_positions": controls,
        "results": results,
        "complete_validated_audit": all_success,
        "independent_bin_validity": "NOT_ASSESSED",
    }


@app.local_entrypoint()
def main():
    result = audit_long_range_visibility.remote()

    print("\nHCT116 LONG-RANGE BIN VISIBILITY AUDIT")
    print("=" * 65)

    print("Source checksum:", result["source_checksum"])
    print("Resolution:", result["resolution_bp"])
    print("Normalization:", result["normalization"])

    print(
        "Maximum separation:",
        result["maximum_separation_bp"],
    )

    print(
        "Comparison positions:",
        result["comparison_positions"],
    )

    for item in result["results"]:
        print(
            f"\n{item['group'].upper()} BIN: "
            f"chr10:{item['position_bp']}"
        )

        print("Status:", item["status"])
        print(
            "Original 0.5-32kb support:",
            item["source_32kb_expected"],
        )
        print(
            "Original-data agreement:",
            item["source_32kb_agreement"],
        )

        if item["status"] == "QUERY_COMPLETED":
            for band, values in item["bands"].items():
                print(f"  {band}: {values}")

            print(
                "Total distinct pairs:",
                item["total_pairs"],
            )
        else:
            print(
                "Error details:",
                {
                    key: value
                    for key, value in item.items()
                    if key in (
                        "returncode",
                        "stderr_tail",
                        "error",
                        "error_type",
                    )
                },
            )

    print(
        "\nComplete validated audit:",
        result["complete_validated_audit"],
    )

    print(
        "Independent bin validity:",
        result["independent_bin_validity"],
    )

    print("\nRead-only long-range audit finished")
