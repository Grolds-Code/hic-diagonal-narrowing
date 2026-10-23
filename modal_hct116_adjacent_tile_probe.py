"""Validate adjacent midpoint-owned Hi-C tiles on Modal."""

import json
import subprocess
import sys

import modal

app = modal.App("hic-diagonal-narrowing-adjacent-tiles")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("libcurl4-openssl-dev", "zlib1g-dev")
    .pip_install("hic-straw==1.3.1")
)

WORKER_CODE = r'''
import json
import math
import hicstraw

URL = (
    "https://huggingface.co/buckets/"
    "FatimaInstitute/3d-genomics-prediction/"
    "resolve/HCT116/ENCFF573OPJ.hic?download=true"
)

START = 20_000_000
SEAM = 20_100_000
END = 20_200_000
MAX_SEP = 32_000
RES = 500

hic = hicstraw.HiCFile(URL)
mzd = hic.getMatrixZoomData(
    "chr10", "chr10",
    "observed", "NONE", "BP", RES
)


def extract(x0, x1, y0, y1):
    records = mzd.getRecords(
        x0, x1 - 1, y0, y1 - 1
    )

    result = {}

    for rec in records:
        x = int(rec.binX)
        y = int(rec.binY)
        count = float(rec.counts)

        if not (
            x0 <= x < x1
            and y0 <= y < y1
        ):
            raise RuntimeError("Out-of-range record")

        if x > y:
            raise RuntimeError("Unexpected triangle")

        if x % RES or y % RES:
            raise RuntimeError("Unaligned coordinate")

        if not math.isfinite(count) or count < 0:
            raise RuntimeError("Invalid count")

        key = (x, y)

        if key in result:
            raise RuntimeError("Duplicate sparse record")

        result[key] = count

    return result


def select(pairs, start, end):
    return {
        (x, y): count
        for (x, y), count in pairs.items()
        if (
            2 * start <= x + y < 2 * end
            and 0 <= y - x <= MAX_SEP
        )
    }


def buffered_tile(start, end):
    records = extract(
        start - MAX_SEP,
        end,
        start,
        end + MAX_SEP
    )
    return records, select(records, start, end)


raw_a, tile_a = buffered_tile(START, SEAM)
raw_b, tile_b = buffered_tile(SEAM, END)

# A midpoint belongs to precisely one half-open tile.
overlap = set(tile_a) & set(tile_b)

if overlap:
    raise RuntimeError(
        f"Duplicate ownership: {len(overlap)}"
    )

combined = dict(tile_a)
combined.update(tile_b)

# Continuous reference covering both buffered tiles.
reference_raw = extract(
    START - MAX_SEP,
    END + MAX_SEP,
    START - MAX_SEP,
    END + MAX_SEP
)

reference = select(reference_raw, START, END)

missing = set(reference) - set(combined)
extra = set(combined) - set(reference)

value_mismatches = sum(
    combined[key] != reference[key]
    for key in set(combined) & set(reference)
)

if missing or extra or value_mismatches:
    raise RuntimeError(
        "Adjacent tiles disagree with reference: "
        f"missing={len(missing)}, "
        f"extra={len(extra)}, "
        f"value_mismatches={value_mismatches}"
    )

# Contacts straddling the shared genomic boundary.
crossing = {
    key
    for key in reference
    if key[0] < SEAM <= key[1]
}

crossing_a = crossing & set(tile_a)
crossing_b = crossing & set(tile_b)

if len(crossing_a) + len(crossing_b) != len(crossing):
    raise RuntimeError("Boundary-contact ownership failure")

# Exact boundary midpoints belong to Tile B.
exact_midpoints = {
    key for key in reference
    if key[0] + key[1] == 2 * SEAM
}

if not exact_midpoints.issubset(set(tile_b)):
    raise RuntimeError("Half-open midpoint rule failed")

if not combined:
    raise RuntimeError("No eligible contacts retrieved")

result = {
    "tile_a_interval": [START, SEAM],
    "tile_b_interval": [SEAM, END],
    "max_separation_bp": MAX_SEP,
    "tile_a_raw_records": len(raw_a),
    "tile_b_raw_records": len(raw_b),
    "tile_a_selected_records": len(tile_a),
    "tile_b_selected_records": len(tile_b),
    "combined_unique_records": len(combined),
    "reference_selected_records": len(reference),
    "combined_contact_signal": sum(combined.values()),
    "reference_contact_signal": sum(reference.values()),
    "shared_boundary_crossing_records": len(crossing),
    "crossing_assigned_tile_a": len(crossing_a),
    "crossing_assigned_tile_b": len(crossing_b),
    "exact_boundary_midpoint_records": len(exact_midpoints),
    "duplicate_ownership": len(overlap),
    "missing_records": len(missing),
    "extra_records": len(extra),
    "value_mismatches": value_mismatches,
    "reference_agreement": "PASS",
}

print("RESULT_JSON=" + json.dumps(result), flush=True)
'''


@app.function(
    image=image,
    cpu=1,
    memory=2048,
    timeout=600,
)
def validate_adjacent_tiles():
    proc = subprocess.run(
        [sys.executable, "-c", WORKER_CODE],
        capture_output=True,
        text=True,
        timeout=480,
        check=False,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Worker failed: {proc.returncode}\n"
            f"stdout: {proc.stdout[-3000:]}\n"
            f"stderr: {proc.stderr[-3000:]}"
        )

    lines = [
        line for line in proc.stdout.splitlines()
        if line.startswith("RESULT_JSON=")
    ]

    if len(lines) != 1:
        raise RuntimeError("Missing validation result")

    return json.loads(lines[0].split("=", 1)[1])


@app.local_entrypoint()
def main():
    result = validate_adjacent_tiles.remote()

    print("\nHCT116 ADJACENT-TILE VALIDATION")
    print("=" * 55)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nAdjacent-tile validation: PASS")
