"""Read-only remote UCSC Umap interval extraction."""

import math
import subprocess

import modal

app = modal.App("hic-diagonal-narrowing-umap-range")

BIN_URL = (
    "https://hgdownload.soe.ucsc.edu/admin/exe/"
    "linux.x86_64/bigWigToBedGraph"
)

ROOT = (
    "https://hgdownload.soe.ucsc.edu/gbdb/"
    "hg38/hoffmanMappability/"
)

TRACKS = {
    "k50": ROOT + "k50.Umap.MultiTrackMappability.bw",
    "k100": ROOT + "k100.Umap.MultiTrackMappability.bw",
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

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("ca-certificates", "curl")
    .run_commands(
        f"curl -fL --retry 3 --max-time 120 "
        f"-o /usr/local/bin/bigWigToBedGraph {BIN_URL}",
        "chmod +x /usr/local/bin/bigWigToBedGraph",
    )
)


def summarize_interval(segments, start, end):
    covered = 0
    weighted_sum = 0.0
    zero_bases = 0

    for left, right, value in segments:
        a = max(left, start)
        b = min(right, end)

        if b <= a:
            continue

        length = b - a
        covered += length
        weighted_sum += value * length

        if value == 0:
            zero_bases += length

    size = end - start

    if covered > size:
        raise RuntimeError("Overlapping BigWig intervals")

    return {
        "interval_bp": [start, end],
        "covered_bases": covered,
        "unreported_bases": size - covered,
        "covered_fraction": round(covered / size, 4),
        "mean_over_reported": (
            round(weighted_sum / covered, 4)
            if covered else None
        ),
        "explicit_zero_bases": zero_bases,
    }


@app.function(
    image=image,
    cpu=1,
    memory=1024,
    timeout=600,
)
def inspect_umap_tracks():
    results = {}

    for track, url in TRACKS.items():
        track_results = {
            "url": url,
            "regions": {},
        }

        for name, (start, end) in PEAKS.items():
            command = [
                "/usr/local/bin/bigWigToBedGraph",
                "-chrom=chr10",
                f"-start={start}",
                f"-end={end}",
                url,
                "stdout",
            ]

            try:
                completed = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=120,
                    check=False,
                )

                if completed.returncode != 0:
                    track_results["regions"][name] = {
                        "status": "QUERY_FAILED",
                        "returncode": completed.returncode,
                        "error": completed.stderr[-1000:],
                    }
                    continue

                segments = []

                for line in completed.stdout.splitlines():
                    fields = line.split()

                    if len(fields) != 4:
                        raise RuntimeError(
                            "Unexpected BedGraph output"
                        )

                    chrom, left, right, raw_value = fields

                    if chrom != "chr10":
                        raise RuntimeError(
                            "Unexpected chromosome"
                        )

                    left = int(left)
                    right = int(right)
                    value = float(raw_value)

                    if (
                        left >= right
                        or not math.isfinite(value)
                        or value < 0
                        or value > 1
                    ):
                        raise RuntimeError(
                            "Invalid Umap interval or value"
                        )

                    segments.append((left, right, value))

                segments.sort()

                for previous, current in zip(
                    segments, segments[1:]
                ):
                    if previous[1] > current[0]:
                        raise RuntimeError(
                            "Overlapping source intervals"
                        )

                region_result = {
                    "status": "SUCCESS",
                    "source_intervals": len(segments),
                    "overall": summarize_interval(
                        segments, start, end
                    ),
                }

                if name == "Region B":
                    region_result["bins"] = {}

                    for position in LOW_SUPPORT + COMPARISONS:
                        region_result["bins"][
                            str(position)
                        ] = summarize_interval(
                            segments,
                            position,
                            position + 500,
                        )

                track_results["regions"][name] = region_result

            except Exception as exc:
                track_results["regions"][name] = {
                    "status": "QUERY_ERROR",
                    "error": (
                        f"{type(exc).__name__}: {exc}"
                    ),
                }

        results[track] = track_results

    return {
        "assembly": "GRCh38/hg38",
        "coordinate_system": "0-based half-open",
        "tracks": results,
        "bin_validity": "NOT_DETERMINED",
    }


@app.local_entrypoint()
def main():
    result = inspect_umap_tracks.remote()

    print("\nHCT116 UCSC UMAP MAPPABILITY AUDIT")
    print("=" * 65)

    print("Assembly:", result["assembly"])
    print("Coordinates:", result["coordinate_system"])

    for track, data in result["tracks"].items():
        print(f"\nUMAP TRACK: {track}")
        print("Source:", data["url"])

        for region, record in data["regions"].items():
            print(f"\n  {region}: {record['status']}")

            if record["status"] != "SUCCESS":
                print("  Details:", record.get("error"))
                continue

            print(
                "  Source intervals:",
                record["source_intervals"],
            )
            print(
                "  Overall:",
                record["overall"],
            )

            for position, summary in (
                record.get("bins", {}).items()
            ):
                label = (
                    "LOW"
                    if int(position) in LOW_SUPPORT
                    else "CONTROL"
                )

                print(
                    f"  {label} {position}: {summary}"
                )

    print(
        "\nHi-C bin validity:",
        result["bin_validity"],
    )
    print("Remote Umap interval audit finished")
