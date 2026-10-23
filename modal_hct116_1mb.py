"""Extract a midpoint-owned 1 Mb HCT116 Hi-C contact band on Modal."""

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import modal

app = modal.App("hic-diagonal-narrowing-1mb")

volume = modal.Volume.from_name(
    "hic-diagonal-narrowing-data"
)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("libcurl4-openssl-dev", "zlib1g-dev")
    .pip_install("hic-straw==1.3.1", "numpy==2.5.4")
)

BASENAME = "hct116_chr10_20000000_21000000_500bp_midpoint_32kb"

WORKER_CODE = r'''
import json
import math

import hicstraw
import numpy as np

URL = (
    "https://huggingface.co/buckets/"
    "FatimaInstitute/3d-genomics-prediction/"
    "resolve/HCT116/ENCFF573OPJ.hic?download=true"
)

START = 20_000_000
END = 21_000_000
TILE_BP = 100_000
RES = 500
MAX_SEP = 32_000

hic = hicstraw.HiCFile(URL)

mzd = hic.getMatrixZoomData(
    "chr10", "chr10",
    "observed", "NONE", "BP", RES
)

selected = {}
tile_stats = []

for tile_index, tile_start in enumerate(
    range(START, END, TILE_BP)
):
    tile_end = tile_start + TILE_BP

    x0 = tile_start - MAX_SEP
    x1 = tile_end
    y0 = tile_start
    y1 = tile_end + MAX_SEP

    records = mzd.getRecords(
        x0, x1 - 1,
        y0, y1 - 1
    )

    raw_seen = set()
    tile_selected = {}
    tile_signal = 0.0

    for rec in records:
        x = int(rec.binX)
        y = int(rec.binY)
        count = float(rec.counts)

        if not (x0 <= x < x1 and y0 <= y < y1):
            raise RuntimeError(
                f"Out-of-range record in tile {tile_index}"
            )

        if x > y or x % RES or y % RES:
            raise RuntimeError(
                f"Invalid coordinates in tile {tile_index}"
            )

        if not math.isfinite(count) or count < 0:
            raise RuntimeError("Invalid contact count")

        key = (x, y)

        if key in raw_seen:
            raise RuntimeError("Duplicate raw record")

        raw_seen.add(key)

        midpoint_inside = (
            2 * tile_start <= x + y < 2 * tile_end
        )

        if (
            midpoint_inside
            and 0 <= y - x <= MAX_SEP
        ):
            if key in tile_selected:
                raise RuntimeError(
                    "Duplicate contact within tile"
                )

            tile_selected[key] = count
            tile_signal += count

    overlap = set(selected) & set(tile_selected)

    if overlap:
        raise RuntimeError(
            f"Duplicate ownership across tiles: {len(overlap)}"
        )

    selected.update(tile_selected)

    tile_stats.append({
        "tile_index": tile_index,
        "start_bp": tile_start,
        "end_bp_exclusive": tile_end,
        "raw_query_records": len(records),
        "selected_records": len(tile_selected),
        "selected_signal": tile_signal,
    })

# Regression check against previously validated adjacent tiles.
first_two_count = sum(
    x["selected_records"] for x in tile_stats[:2]
)
first_two_signal = sum(
    x["selected_signal"] for x in tile_stats[:2]
)

if first_two_count != 8166:
    raise RuntimeError(
        f"Two-tile count regression: {first_two_count}"
    )

if not math.isclose(
    first_two_signal, 21638.0, abs_tol=1e-6
):
    raise RuntimeError(
        f"Two-tile signal regression: {first_two_signal}"
    )

if len(tile_stats) != 10 or not selected:
    raise RuntimeError("Incomplete 1 Mb extraction")

if sum(
    x["selected_records"] for x in tile_stats
) != len(selected):
    raise RuntimeError("Tile counts do not reconcile")

items = sorted(selected.items())

bin_x_bp = np.array(
    [key[0] for key, _ in items],
    dtype=np.int64
)

bin_y_bp = np.array(
    [key[1] for key, _ in items],
    dtype=np.int64
)

counts = np.array(
    [value for _, value in items],
    dtype=np.float64
)

if not np.all(
    (2 * START <= bin_x_bp + bin_y_bp)
    & (bin_x_bp + bin_y_bp < 2 * END)
):
    raise RuntimeError("Invalid midpoint ownership")

if not np.all(
    (bin_y_bp >= bin_x_bp)
    & (bin_y_bp - bin_x_bp <= MAX_SEP)
):
    raise RuntimeError("Invalid separation")

np.savez_compressed(
    "/tmp/hct116_1mb_band.npz",
    bin_x_bp=bin_x_bp,
    bin_y_bp=bin_y_bp,
    counts=counts,
)

seam_crossings = {}

for seam in range(
    START + TILE_BP, END, TILE_BP
):
    crossing = (
        (bin_x_bp < seam)
        & (bin_y_bp >= seam)
    )

    seam_crossings[str(seam)] = int(
        np.count_nonzero(crossing)
    )

result = {
    "chromosome": "chr10",
    "start_bp": START,
    "end_bp_exclusive": END,
    "resolution_bp": RES,
    "max_separation_bp": MAX_SEP,
    "tiles_processed": len(tile_stats),
    "selected_records": len(selected),
    "selected_contact_signal": float(counts.sum()),
    "first_two_tile_regression": "PASS",
    "duplicate_ownership": 0,
    "outside_left_records": int(
        np.count_nonzero(bin_x_bp < START)
    ),
    "outside_right_records": int(
        np.count_nonzero(bin_y_bp >= END)
    ),
    "seam_crossings": seam_crossings,
    "tile_stats": tile_stats,
}

print(
    "RESULT_JSON=" + json.dumps(result),
    flush=True
)
'''


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=1200,
)
def extract_one_megabase():
    import numpy as np

    matrix_dir = Path("/research/matrices")
    metadata_dir = Path("/research/metadata")

    destination = matrix_dir / f"{BASENAME}.npz"
    metadata_path = metadata_dir / f"{BASENAME}.json"

    if destination.exists() or metadata_path.exists():
        raise RuntimeError(
            "Existing output detected; refusing to overwrite"
        )

    proc = subprocess.run(
        [sys.executable, "-c", WORKER_CODE],
        capture_output=True,
        text=True,
        timeout=900,
        check=False,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Extraction worker exit: {proc.returncode}\n"
            f"stdout: {proc.stdout[-4000:]}\n"
            f"stderr: {proc.stderr[-4000:]}"
        )

    lines = [
        line for line in proc.stdout.splitlines()
        if line.startswith("RESULT_JSON=")
    ]

    if len(lines) != 1:
        raise RuntimeError("Missing worker result")

    result = json.loads(lines[0].split("=", 1)[1])

    temporary = Path("/tmp/hct116_1mb_band.npz")

    with np.load(temporary, allow_pickle=False) as data:
        x = data["bin_x_bp"]
        y = data["bin_y_bp"]
        counts = data["counts"]

    n = result["selected_records"]

    if any(arr.shape != (n,) for arr in (x, y, counts)):
        raise RuntimeError("Unexpected sparse array dimensions")

    if not np.all(np.isfinite(counts)) or np.any(counts < 0):
        raise RuntimeError("Invalid stored contact counts")

    if len(set(zip(x.tolist(), y.tolist()))) != n:
        raise RuntimeError("Duplicate pairs in stored arrays")

    if not np.all(
        (x + y >= 40_000_000)
        & (x + y < 42_000_000)
        & (y >= x)
        & (y - x <= 32_000)
    ):
        raise RuntimeError("Sparse array coordinate validation failed")

    if not np.isclose(
        counts.sum(),
        result["selected_contact_signal"],
    ):
        raise RuntimeError("Contact total mismatch")

    digest = hashlib.sha256(
        temporary.read_bytes()
    ).hexdigest()

    metadata = {
        "dataset": "ENCFF573OPJ.hic",
        "source_url": (
            "https://huggingface.co/buckets/"
            "FatimaInstitute/3d-genomics-prediction/"
            "resolve/HCT116/ENCFF573OPJ.hic?download=true"
        ),
        "cell_line": "HCT116",
        "assembly": "GRCh38",
        "normalization": "NONE",
        "coordinate_convention": "zero_based_half_open",
        "representation": "sparse_upper_triangle_midpoint_band",
        "selection": (
            "midpoint in target region; "
            "0 <= separation <= 32000 bp"
        ),
        "counts_dtype": "float64",
        "missing_bin_status": "not_assessed",
        "sha256": digest,
        **result,
    }

    matrix_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    staging = destination.with_suffix(".partial")

    shutil.copyfile(temporary, staging)
    os.replace(staging, destination)

    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )

    volume.commit()

    return {
        "data_path": str(destination),
        "metadata_path": str(metadata_path),
        "size_bytes": destination.stat().st_size,
        "sha256": digest,
        **result,
    }


@app.local_entrypoint()
def main():
    result = extract_one_megabase.remote()

    print("\nHCT116 COMPLETE 1 MB EXTRACTION")
    print("=" * 55)

    for key, value in result.items():
        if key == "tile_stats":
            print("\nPER-TILE RESULTS")
            for tile in value:
                print(
                    f"Tile {tile['tile_index']:02d}: "
                    f"{tile['selected_records']} records, "
                    f"signal={tile['selected_signal']:.1f}"
                )
        elif key == "seam_crossings":
            print("\nBOUNDARY-CROSSING RECORDS")
            for seam, count in value.items():
                print(f"{seam}: {count}")
        else:
            print(f"{key}: {value}")

    print("\n1 Mb cloud extraction: PASS")
