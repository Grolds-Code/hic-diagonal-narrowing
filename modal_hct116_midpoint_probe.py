"""Validate two-sided, midpoint-based Hi-C tiling."""

import json
import subprocess
import sys

import modal

app = modal.App("hic-diagonal-narrowing-midpoint-probe")

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
END = 20_100_000
MAX_SEP = 32_000
RESOLUTION = 500

LEFT = START - MAX_SEP
RIGHT = END + MAX_SEP

hic = hicstraw.HiCFile(URL)

mzd = hic.getMatrixZoomData(
    "chr10", "chr10",
    "observed", "NONE", "BP", RESOLUTION
)


def extract(x_start, x_end, y_start, y_end):
    records = mzd.getRecords(
        x_start, x_end - 1,
        y_start, y_end - 1
    )

    pairs = {}

    for rec in records:
        x = int(rec.binX)
        y = int(rec.binY)
        value = float(rec.counts)

        if not (
            x_start <= x < x_end
            and y_start <= y < y_end
        ):
            raise RuntimeError(
                f"Record outside requested rectangle: {x}, {y}"
            )

        if x > y:
            raise RuntimeError(
                "Unexpected lower-triangle record"
            )

        if x % RESOLUTION or y % RESOLUTION:
            raise RuntimeError("Misaligned coordinates")

        if not math.isfinite(value) or value < 0:
            raise RuntimeError("Invalid contact count")

        key = (x, y)

        if key in pairs:
            raise RuntimeError(
                f"Duplicate contact pair: {key}"
            )

        pairs[key] = value

    return pairs


def select_by_midpoint(pairs):
    selected = {}

    for (x, y), value in pairs.items():
        separation = y - x

        # Integer arithmetic avoids midpoint rounding.
        midpoint_inside = (
            2 * START <= x + y < 2 * END
        )

        if (
            midpoint_inside
            and separation <= MAX_SEP
        ):
            selected[(x, y)] = value

    return selected


# Buffered rectangular extraction:
# first endpoint may lie in the left buffer;
# second endpoint may lie in the right buffer.
rectangular = extract(
    LEFT, END,
    START, RIGHT
)

# Independent reference:
# retrieve a complete square spanning both buffers.
full_square = extract(
    LEFT, RIGHT,
    LEFT, RIGHT
)

tile_contacts = select_by_midpoint(rectangular)
reference_contacts = select_by_midpoint(full_square)

if tile_contacts != reference_contacts:
    missing = set(reference_contacts) - set(tile_contacts)
    extra = set(tile_contacts) - set(reference_contacts)

    mismatched = [
        key for key in
        set(tile_contacts) & set(reference_contacts)
        if tile_contacts[key] != reference_contacts[key]
    ]

    raise RuntimeError(
        "Midpoint tiling disagrees with square reference: "
        f"missing={len(missing)}, "
        f"extra={len(extra)}, "
        f"mismatched={len(mismatched)}"
    )

left_crossing = sum(
    x < START
    for x, y in tile_contacts
)

right_crossing = sum(
    y >= END
    for x, y in tile_contacts
)

internal = sum(
    x >= START and y < END
    for x, y in tile_contacts
)

if (
    left_crossing + right_crossing + internal
    != len(tile_contacts)
):
    raise RuntimeError(
        "Contact classification does not reconcile"
    )

result = {
    "core_interval": [START, END],
    "buffered_interval": [LEFT, RIGHT],
    "max_separation_bp": MAX_SEP,
    "rectangular_query_records": len(rectangular),
    "square_reference_records": len(full_square),
    "midpoint_selected_records": len(tile_contacts),
    "internal_records": internal,
    "left_crossing_records": left_crossing,
    "right_crossing_records": right_crossing,
    "midpoint_selected_signal": sum(
        tile_contacts.values()
    ),
    "reference_agreement": "PASS",
    "duplicate_check": "PASS",
}

print(
    "RESULT_JSON=" + json.dumps(result),
    flush=True
)
'''


@app.function(
    image=image,
    cpu=1,
    memory=2048,
    timeout=600,
)
def validate_midpoint_tiling():
    proc = subprocess.run(
        [sys.executable, "-c", WORKER_CODE],
        capture_output=True,
        text=True,
        timeout=420,
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
        raise RuntimeError(
            "Missing structured validation result"
        )

    return json.loads(lines[0].split("=", 1)[1])


@app.local_entrypoint()
def main():
    result = validate_midpoint_tiling.remote()

    print("\nHCT116 TWO-SIDED MIDPOINT TILING")
    print("=" * 55)

    for key, value in result.items():
        print(f"{key}: {value}")

    print("\nMidpoint tiling validation: PASS")
