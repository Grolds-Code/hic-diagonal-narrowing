"""Isolated small-region HCT116 sparse extraction on Modal."""

import json
import subprocess
import sys

import modal

app = modal.App("hic-diagonal-narrowing-sparse-probe")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("libcurl4-openssl-dev", "zlib1g-dev")
    .pip_install("hic-straw==1.3.1")
)

WORKER_CODE = r'''
import json
import math
import hicstraw

url = (
    "https://huggingface.co/buckets/"
    "FatimaInstitute/3d-genomics-prediction/"
    "resolve/HCT116/ENCFF573OPJ.hic?download=true"
)

chrom = "chr10"
start = 20_000_000
end = 20_100_000
resolution = 500

hic = hicstraw.HiCFile(url)

mzd = hic.getMatrixZoomData(
    chrom, chrom, "observed", "NONE", "BP", resolution
)

# Query the last included base, rather than the exclusive end.
records = mzd.getRecords(
    start, end - 1, start, end - 1
)

outside = 0
misaligned = 0
invalid_counts = 0
upper = 0
lower = 0
diagonal = 0
total_signal = 0.0
coordinates = []

for rec in records:
    x = int(rec.binX)
    y = int(rec.binY)
    count = float(rec.counts)

    if not (start <= x < end and start <= y < end):
        outside += 1

    if x % resolution or y % resolution:
        misaligned += 1

    if not math.isfinite(count) or count < 0:
        invalid_counts += 1
    else:
        total_signal += count

    if x < y:
        upper += 1
    elif x > y:
        lower += 1
    else:
        diagonal += 1

    coordinates.append((x, y))

result = {
    "chromosome": chrom,
    "start_bp": start,
    "end_bp_exclusive": end,
    "resolution_bp": resolution,
    "expected_bins": (end - start) // resolution,
    "sparse_records": len(records),
    "upper_records": upper,
    "lower_records": lower,
    "diagonal_records": diagonal,
    "outside_interval": outside,
    "misaligned_records": misaligned,
    "invalid_counts": invalid_counts,
    "sum_record_counts": total_signal,
    "coordinate_range": [
        min(min(x, y) for x, y in coordinates)
        if coordinates else None,
        max(max(x, y) for x, y in coordinates)
        if coordinates else None,
    ],
}

print("PROBE_JSON=" + json.dumps(result), flush=True)

if not records or outside or misaligned or invalid_counts:
    raise RuntimeError("Sparse-record validation failed")
'''


@app.function(
    image=image,
    cpu=1,
    memory=2048,
    timeout=360,
)
def probe_sparse_records():
    """Isolate native hic-straw execution in a child process."""

    proc = subprocess.run(
        [sys.executable, "-c", WORKER_CODE],
        capture_output=True,
        text=True,
        timeout=240,
        check=False,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Extraction worker exited: {proc.returncode}\n"
            f"stdout: {proc.stdout[-3000:]}\n"
            f"stderr: {proc.stderr[-3000:]}"
        )

    lines = [
        line for line in proc.stdout.splitlines()
        if line.startswith("PROBE_JSON=")
    ]

    if not lines:
        raise RuntimeError(
            "No structured extraction result returned"
        )

    return json.loads(lines[-1].split("=", 1)[1])


@app.local_entrypoint()
def main():
    result = probe_sparse_records.remote()

    print("\nHCT116 SPARSE EXTRACTION PROBE")
    print("=" * 45)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nRemote sparse extraction: PASS")
