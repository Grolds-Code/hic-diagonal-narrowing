"""Validate boundary-aware Hi-C extraction on Modal."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-tile-probe")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("libcurl4-openssl-dev", "zlib1g-dev")
    .pip_install("hic-straw==1.3.1", "numpy==2.5.4")
)

BASENAME = "hct116_chr10_20000000_20100000_500bp"

WORKER_CODE = r'''
import json
import math
import hicstraw

url = (
    "https://huggingface.co/buckets/"
    "FatimaInstitute/3d-genomics-prediction/"
    "resolve/HCT116/ENCFF573OPJ.hic?download=true"
)

start = 20_000_000
tile_end = 20_100_000
buffer_end = 20_132_000
resolution = 500
max_separation = 32_000

hic = hicstraw.HiCFile(url)
mzd = hic.getMatrixZoomData(
    "chr10", "chr10", "observed", "NONE", "BP",
    resolution
)

records = mzd.getRecords(
    start, tile_end - 1,
    start, buffer_end - 1
)

within_count = 0
within_signal = 0.0
crossing_count = 0
near_crossing_count = 0
near_crossing_signal = 0.0
seen = set()

for rec in records:
    x = int(rec.binX)
    y = int(rec.binY)
    value = float(rec.counts)

    if not (start <= x < tile_end):
        raise RuntimeError(f"Unexpected x coordinate: {x}")

    if not (start <= y < buffer_end):
        raise RuntimeError(f"Unexpected y coordinate: {y}")

    if x > y:
        raise RuntimeError("Unexpected lower triangle")

    if x % resolution or y % resolution:
        raise RuntimeError("Unaligned genomic coordinate")

    if not math.isfinite(value) or value < 0:
        raise RuntimeError("Invalid contact count")

    if (x, y) in seen:
        raise RuntimeError("Duplicate contact record")

    seen.add((x, y))

    if y < tile_end:
        within_count += 1
        within_signal += value
    else:
        crossing_count += 1

        if y - x <= max_separation:
            near_crossing_count += 1
            near_crossing_signal += value

print("RESULT_JSON=" + json.dumps({
    "query_records": len(records),
    "within_tile_records": within_count,
    "within_tile_signal": within_signal,
    "cross_boundary_records": crossing_count,
    "cross_boundary_near_32kb_records": near_crossing_count,
    "cross_boundary_near_32kb_signal": near_crossing_signal,
}), flush=True)
'''


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=1,
    memory=2048,
    timeout=420,
)
def verify_tile_boundary():
    metadata_path = Path(
        f"/research/metadata/{BASENAME}.json"
    )
    matrix_path = Path(
        f"/research/matrices/{BASENAME}.npz"
    )

    metadata = json.loads(metadata_path.read_text())

    digest = hashlib.sha256(
        matrix_path.read_bytes()
    ).hexdigest()

    if digest != metadata["sha256"]:
        raise RuntimeError("Stored matrix checksum mismatch")

    proc = subprocess.run(
        [sys.executable, "-c", WORKER_CODE],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Worker exit code: {proc.returncode}\n"
            f"stdout: {proc.stdout[-2000:]}\n"
            f"stderr: {proc.stderr[-2000:]}"
        )

    lines = [
        line for line in proc.stdout.splitlines()
        if line.startswith("RESULT_JSON=")
    ]

    if len(lines) != 1:
        raise RuntimeError("Missing extraction results")

    result = json.loads(lines[0].split("=", 1)[1])

    if (
        result["within_tile_records"]
        != metadata["sparse_records"]
    ):
        raise RuntimeError(
            "Tiled extraction disagrees with original "
            "100 kb sparse-record count"
        )

    if (
        abs(
            result["within_tile_signal"]
            - metadata["unique_contact_signal"]
        ) > 1e-6
    ):
        raise RuntimeError(
            "Tiled extraction disagrees with original "
            "100 kb contact signal"
        )

    result["original_100kb_agreement"] = "PASS"

    return result


@app.local_entrypoint()
def main():
    result = verify_tile_boundary.remote()

    print("\nHCT116 BOUNDARY-AWARE TILE PROBE")
    print("=" * 50)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nBoundary-aware extraction test: PASS")
