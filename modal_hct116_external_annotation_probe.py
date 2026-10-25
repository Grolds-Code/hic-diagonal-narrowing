"""Read-only GRCh38 blacklist and sequence-mappability audit."""

import gzip
import hashlib
import math
from urllib.request import Request, urlopen

import modal

app = modal.App("hic-diagonal-narrowing-external-annotations")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("libcurl4-openssl-dev", "zlib1g-dev")
    .pip_install("numpy==2.5.4", "pyBigWig==0.3.26")
)

BLACKLIST_URL = (
    "https://raw.githubusercontent.com/Boyle-Lab/"
    "Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
)

TRACK_ROOT = (
    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/"
    "hoffmanMappability/"
)

TRACKS = {
    "k50": TRACK_ROOT + "k50.Umap.MultiTrackMappability.bw",
    "k100": TRACK_ROOT + "k100.Umap.MultiTrackMappability.bw",
}

PEAKS = {
    "Region A": (20296000, 20328000),
    "Region B": (20736000, 20768000),
}

LOW_SUPPORT = [
    20752000, 20753000, 20754000, 20755000,
    20755500, 20756000, 20756500, 20757000,
]

COMPARISONS = [
    20738000, 20744000, 20764000,
]


def overlap(a, b):
    return max(a[0], b[0]) < min(a[1], b[1])


def summarize_values(values):
    if values is None:
        return {
            "status": "NO_TRACK_DATA",
        }

    valid = [
        float(value)
        for value in values
        if value is not None and math.isfinite(float(value))
    ]

    count = len(values)

    return {
        "bases": count,
        "reported_bases": len(valid),
        "unreported_fraction": round(
            1 - len(valid) / count, 4
        ) if count else None,
        "mean_reported": round(
            sum(valid) / len(valid), 4
        ) if valid else None,
        "reported_zero_fraction": round(
            sum(v == 0 for v in valid) / len(valid), 4
        ) if valid else None,
    }


@app.function(
    image=image,
    cpu=1,
    memory=1536,
    timeout=600,
)
def audit_external_annotations():
    import pyBigWig

    output = {
        "assembly": "GRCh38/hg38",
        "coordinate_convention": "0-based half-open",
        "blacklist": {},
        "mappability": {},
    }

    # The blacklist is a small compressed BED file.
    try:
        request = Request(
            BLACKLIST_URL,
            headers={"User-Agent": "HiC-Research-QC/1.0"},
        )

        with urlopen(request, timeout=45) as response:
            raw = response.read()

        if raw[:2] != b"\x1f\x8b":
            raise RuntimeError("Expected gzip-compressed BED")

        lines = gzip.decompress(raw).decode().splitlines()
        intervals = []

        for line in lines:
            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3 or fields[0] != "chr10":
                continue

            left = int(fields[1])
            right = int(fields[2])

            if right <= left:
                raise RuntimeError("Invalid blacklist interval")

            intervals.append((left, right))

        output["blacklist"] = {
            "status": "SUCCESS",
            "url": BLACKLIST_URL,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "chr10_blacklist_intervals": len(intervals),
            "peak_overlaps": {
                name: [
                    list(interval)
                    for interval in intervals
                    if overlap(peak, interval)
                ]
                for name, peak in PEAKS.items()
            },
            "low_support_bin_overlaps": {
                str(pos): any(
                    overlap((pos, pos + 500), interval)
                    for interval in intervals
                )
                for pos in LOW_SUPPORT
            },
        }

    except Exception as exc:
        output["blacklist"] = {
            "status": "FAILED",
            "error": f"{type(exc).__name__}: {exc}",
        }

    # Remote BigWig queries use HTTP byte ranges.
    # No full genome-wide track is downloaded.
    remote_enabled = bool(pyBigWig.remote)

    for name, url in TRACKS.items():
        if not remote_enabled:
            output["mappability"][name] = {
                "status": "REMOTE_ACCESS_UNAVAILABLE",
                "url": url,
            }
            continue

        bw = None

        try:
            bw = pyBigWig.open(url)

            if bw is None:
                raise RuntimeError("Could not open remote BigWig")

            if "chr10" not in bw.chroms():
                raise RuntimeError("chr10 absent from BigWig")

            peak_values = {
                peak_name: bw.values(
                    "chr10",
                    left,
                    right,
                    numpy=False,
                )
                for peak_name, (left, right) in PEAKS.items()
            }

            region_b_start = PEAKS["Region B"][0]
            region_b_values = peak_values["Region B"]

            examined_bins = {}

            for pos in LOW_SUPPORT + COMPARISONS:
                offset = pos - region_b_start

                if not 0 <= offset <= len(region_b_values) - 500:
                    raise RuntimeError("Bin outside Region B peak")

                examined_bins[str(pos)] = summarize_values(
                    region_b_values[offset:offset + 500]
                )

            output["mappability"][name] = {
                "status": "SUCCESS",
                "url": url,
                "peak_summaries": {
                    peak_name: summarize_values(values)
                    for peak_name, values in peak_values.items()
                },
                "examined_500bp_bins": examined_bins,
            }

        except Exception as exc:
            output["mappability"][name] = {
                "status": "FAILED",
                "url": url,
                "error": f"{type(exc).__name__}: {exc}",
            }

        finally:
            if bw is not None:
                bw.close()

    output["pybigwig_remote_enabled"] = remote_enabled

    return output


@app.local_entrypoint()
def main():
    result = audit_external_annotations.remote()

    print("\nHCT116 INDEPENDENT GENOMIC ANNOTATION AUDIT")
    print("=" * 65)

    print("Assembly:", result["assembly"])
    print("Coordinates:", result["coordinate_convention"])
    print(
        "Remote BigWig support:",
        result["pybigwig_remote_enabled"],
    )

    print("\nENCODE BLACKLIST")
    for key, value in result["blacklist"].items():
        print(f"{key}: {value}")

    for track, record in result["mappability"].items():
        print(f"\nUMAP TRACK: {track}")
        print("Status:", record["status"])
        print("Source:", record["url"])

        if record["status"] != "SUCCESS":
            print("Error:", record.get("error"))
            continue

        print("\nPeak summaries:")
        for peak, values in record["peak_summaries"].items():
            print(f"  {peak}: {values}")

        print("\nLow-support and comparison bins:")
        for pos, values in record[
            "examined_500bp_bins"
        ].items():
            group = (
                "LOW"
                if int(pos) in LOW_SUPPORT
                else "CONTROL"
            )
            print(f"  {group} {pos}: {values}")

    print("\nExternal annotation audit finished")
    print("Hi-C-specific bin validity: NOT DETERMINED")
