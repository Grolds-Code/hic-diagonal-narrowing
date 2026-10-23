"""Reconstruct real HCT116 Hi-C contacts on Modal.

Dataset: ENCFF573OPJ.hic
Assembly: GRCh38
Coordinates: zero-based, half-open
Resolution: 500 bp
"""

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import modal


app = modal.App("hic-diagonal-narrowing-matrix")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("libcurl4-openssl-dev", "zlib1g-dev")
    .pip_install("hic-straw==1.3.1", "numpy==2.5.4")
)

BASENAME = "hct116_chr10_20000000_20100000_500bp"

URL = (
    "https://huggingface.co/buckets/"
    "FatimaInstitute/3d-genomics-prediction/"
    "resolve/HCT116/ENCFF573OPJ.hic?download=true"
)


WORKER_CODE = r'''
import json
import math

import hicstraw
import numpy as np

url = (
    "https://huggingface.co/buckets/"
    "FatimaInstitute/3d-genomics-prediction/"
    "resolve/HCT116/ENCFF573OPJ.hic?download=true"
)

start = 20_000_000
end = 20_100_000
resolution = 500
n_bins = (end - start) // resolution

hic = hicstraw.HiCFile(url)

mzd = hic.getMatrixZoomData(
    "chr10", "chr10",
    "observed", "NONE", "BP", resolution
)

records = mzd.getRecords(
    start, end - 1, start, end - 1
)

matrix = np.zeros(
    (n_bins, n_bins),
    dtype=np.float64
)

seen = set()
total_signal = 0.0
diagonal_records = 0

for record in records:
    x = int(record.binX)
    y = int(record.binY)
    count = float(record.counts)

    if not (start <= x < end and start <= y < end):
        raise ValueError("Record outside requested interval")

    if x % resolution or y % resolution:
        raise ValueError("Misaligned genomic coordinate")

    if not math.isfinite(count) or count < 0:
        raise ValueError("Invalid contact count")

    i = (x - start) // resolution
    j = (y - start) // resolution

    if i > j:
        raise ValueError("Unexpected lower-triangle record")

    if (i, j) in seen:
        raise ValueError("Duplicate contact record")

    seen.add((i, j))

    matrix[i, j] = count
    matrix[j, i] = count

    total_signal += count

    if i == j:
        diagonal_records += 1

if len(records) == 0:
    raise RuntimeError("No contact records retrieved")

assert matrix.shape == (200, 200)
assert np.array_equal(matrix, matrix.T)
assert np.all(np.isfinite(matrix))
assert np.all(matrix >= 0)

unique_signal = float(
    np.triu(matrix).sum()
)

if not np.isclose(unique_signal, total_signal):
    raise RuntimeError("Contact signal mismatch")

np.savez_compressed(
    "/tmp/hct116_matrix.npz",
    contacts=matrix,
)

result = {
    "shape": list(matrix.shape),
    "sparse_records": len(records),
    "diagonal_records": diagonal_records,
    "nonzero_cells": int(np.count_nonzero(matrix)),
    "nonzero_fraction": float(
        np.count_nonzero(matrix) / matrix.size
    ),
    "unique_contact_signal": unique_signal,
    "symmetric": True,
    "coordinate_validation": "PASS",
}

print("WORKER_JSON=" + json.dumps(result), flush=True)
'''


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def reconstruct_hct116_matrix():
    """Extract, validate, and persist a contact matrix."""

    import numpy as np

    proc = subprocess.run(
        [sys.executable, "-c", WORKER_CODE],
        capture_output=True,
        text=True,
        timeout=420,
        check=False,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Extraction failed: {proc.returncode}\n"
            f"stdout: {proc.stdout[-3000:]}\n"
            f"stderr: {proc.stderr[-3000:]}"
        )

    lines = [
        line for line in proc.stdout.splitlines()
        if line.startswith("WORKER_JSON=")
    ]

    if len(lines) != 1:
        raise RuntimeError("Missing worker metadata")

    result = json.loads(lines[0].split("=", 1)[1])

    temporary = Path("/tmp/hct116_matrix.npz")

    with np.load(temporary, allow_pickle=False) as data:
        matrix = data["contacts"]

    if matrix.shape != (200, 200):
        raise RuntimeError("Unexpected matrix dimensions")

    if not np.array_equal(matrix, matrix.T):
        raise RuntimeError("Matrix symmetry check failed")

    if not np.all(np.isfinite(matrix)):
        raise RuntimeError("Non-finite matrix values")

    if np.any(matrix < 0):
        raise RuntimeError("Negative contact counts")

    matrix_dir = Path("/research/matrices")
    metadata_dir = Path("/research/metadata")

    matrix_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    destination = matrix_dir / f"{BASENAME}.npz"
    metadata_path = metadata_dir / f"{BASENAME}.json"

    if destination.exists() or metadata_path.exists():
        raise RuntimeError(
            "Output already exists; refusing to overwrite"
        )

    digest = hashlib.sha256(
        temporary.read_bytes()
    ).hexdigest()

    metadata = {
        "dataset": "ENCFF573OPJ.hic",
        "source_url": URL,
        "cell_line": "HCT116",
        "genome_assembly": "GRCh38",
        "chromosome": "chr10",
        "start_bp": 20_000_000,
        "end_bp_exclusive": 20_100_000,
        "resolution_bp": 500,
        "normalization": "NONE",
        "coordinate_convention": "zero_based_half_open",
        "matrix_dtype": "float64",
        "missing_bin_status": "not_assessed",
        "sha256": digest,
        **result,
    }

    staging = destination.with_suffix(".partial")

    shutil.copyfile(temporary, staging)
    os.replace(staging, destination)

    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    volume.commit()

    return {
        "matrix_path": str(destination),
        "metadata_path": str(metadata_path),
        "matrix_size_bytes": destination.stat().st_size,
        **result,
    }


@app.local_entrypoint()
def main():
    result = reconstruct_hct116_matrix.remote()

    print("\nHCT116 MATRIX RECONSTRUCTION")
    print("=" * 45)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nReal Hi-C matrix reconstruction: PASS")
