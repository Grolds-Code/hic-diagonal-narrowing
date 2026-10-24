"""Read-only, subprocess-isolated HCT116 normalization probe."""

import json
import subprocess
import sys

import modal

app = modal.App("hic-diagonal-narrowing-normalization-probe")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install(
        "libcurl4-openssl-dev",
        "zlib1g-dev",
    )
    .pip_install("hic-straw==1.3.1")
)

NORMALIZATIONS = (
    "NONE",
    "VC",
    "VC_SQRT",
    "KR",
    "SCALE",
)

CHILD_CODE = r'''
import hashlib
import json
import math
import sys

import hicstraw

normalization = sys.argv[1]

url = (
    "https://huggingface.co/buckets/FatimaInstitute/"
    "3d-genomics-prediction/resolve/HCT116/"
    "ENCFF573OPJ.hic?download=true"
)

start = 20_744_000
end = 20_764_000

problem_bins = [
    20_752_000,
    20_753_000,
    20_754_000,
    20_755_500,
    20_756_000,
    20_756_500,
]

try:
    hic = hicstraw.HiCFile(url)

    matrix = hic.getMatrixZoomData(
        "chr10",
        "chr10",
        "observed",
        normalization,
        "BP",
        500,
    )

    records = matrix.getRecords(
        start,
        end - 1,
        start,
        end - 1,
    )

    pairs = []
    values = []
    problem_bin_hits = {
        str(position): 0
        for position in problem_bins
    }

    for record in records:
        x = int(record.binX)
        y = int(record.binY)
        value = float(record.counts)

        pairs.append((min(x, y), max(x, y)))
        values.append(value)

        for position in problem_bins:
            if x == position or y == position:
                problem_bin_hits[str(position)] += 1

    canonical_pairs = sorted(set(pairs))

    digest = hashlib.sha256(
        "\n".join(
            f"{x},{y}"
            for x, y in canonical_pairs
        ).encode("utf-8")
    ).hexdigest()

    finite_values = [
        value for value in values
        if math.isfinite(value)
    ]

    result = {
        "normalization": normalization,
        "status": "QUERY_COMPLETED",
        "records": len(records),
        "distinct_pairs": len(canonical_pairs),
        "coordinate_digest": digest,
        "finite_records": len(finite_values),
        "nonfinite_records": (
            len(values) - len(finite_values)
        ),
        "positive_finite_records": sum(
            value > 0 for value in finite_values
        ),
        "problem_bin_record_counts": problem_bin_hits,
        "vector_values_directly_accessed": False,
    }

except Exception as exc:
    result = {
        "normalization": normalization,
        "status": "QUERY_ERROR",
        "error_type": type(exc).__name__,
        "error_message": str(exc)[:500],
    }

print(json.dumps(result, sort_keys=True))
'''


@app.function(
    image=image,
    cpu=1,
    memory=2048,
    timeout=600,
)
def inspect_normalizations():
    outcomes = []

    for normalization in NORMALIZATIONS:
        try:
            completed = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    CHILD_CODE,
                    normalization,
                ],
                capture_output=True,
                text=True,
                timeout=90,
                check=False,
            )

            if completed.returncode != 0:
                result = {
                    "normalization": normalization,
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

                if not lines:
                    result = {
                        "normalization": normalization,
                        "status": "NO_CHILD_OUTPUT",
                    }
                else:
                    try:
                        result = json.loads(lines[-1])
                    except json.JSONDecodeError:
                        result = {
                            "normalization": normalization,
                            "status": "UNPARSEABLE_OUTPUT",
                            "stdout_tail": (
                                completed.stdout[-500:]
                            ),
                        }

        except subprocess.TimeoutExpired:
            result = {
                "normalization": normalization,
                "status": "TIMEOUT",
            }

        outcomes.append(result)

    return {
        "genome": "hg38",
        "chromosome": "chr10",
        "resolution_bp": 500,
        "query_start_bp": 20_744_000,
        "query_end_bp": 20_764_000,
        "results": outcomes,
        "normalization_vector_status": (
            "NOT_DIRECTLY_ACCESSED"
        ),
        "missing_bin_status": (
            "NOT_DETERMINED"
        ),
    }


@app.local_entrypoint()
def main():
    result = inspect_normalizations.remote()

    print("\nHCT116 500 BP NORMALIZATION QUERY PROBE")
    print("=" * 65)

    print("Genome:", result["genome"])
    print("Chromosome:", result["chromosome"])
    print("Resolution:", result["resolution_bp"])

    print(
        "Interval:",
        result["query_start_bp"],
        result["query_end_bp"],
    )

    for outcome in result["results"]:
        print(
            f"\nNORMALIZATION: "
            f"{outcome['normalization']}"
        )

        for key, value in outcome.items():
            if key != "normalization":
                print(f"  {key}: {value}")

    print(
        "\nDirect vector access:",
        result["normalization_vector_status"],
    )

    print(
        "Missing-bin status:",
        result["missing_bin_status"],
    )

    print("\nRead-only normalization probe completed")
