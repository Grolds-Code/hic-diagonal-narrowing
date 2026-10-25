"""Read-only HCT116 mapping-filter comparison at Region B."""

import json
import subprocess
import sys
from urllib.request import Request, urlopen

import modal

app = modal.App("hic-diagonal-narrowing-mapq-probe")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("libcurl4-openssl-dev", "zlib1g-dev")
    .pip_install("hic-straw==1.3.1")
)

CURRENT_URL = (
    "https://huggingface.co/buckets/FatimaInstitute/"
    "3d-genomics-prediction/resolve/HCT116/"
    "ENCFF573OPJ.hic?download=true"
)

ACCESSIONS = (
    "ENCFF573OPJ",
    "ENCFF726YRX",
)

PEAK_START = 20_736_000
PEAK_END = 20_768_000

LOW_SUPPORT = (
    20752000, 20753000, 20754000, 20755000,
    20755500, 20756000, 20756500, 20757000,
)

CHILD_CODE = r'''
import json
import math
import sys

import hicstraw

accession = sys.argv[1]
url = sys.argv[2]

start = 20_736_000
end = 20_768_000
max_sep = 32_000

low_bins = [
    20752000, 20753000, 20754000, 20755000,
    20755500, 20756000, 20756500, 20757000,
]

result = {
    "accession": accession,
    "status": "STARTED",
}

try:
    hic = hicstraw.HiCFile(url)

    result["genome"] = hic.getGenomeID()
    result["resolution_500bp_available"] = (
        500 in hic.getResolutions()
    )

    if not result["resolution_500bp_available"]:
        result["status"] = "RESOLUTION_UNAVAILABLE"
    else:
        matrix = hic.getMatrixZoomData(
            "chr10", "chr10",
            "observed", "NONE", "BP", 500,
        )

        records = matrix.getRecords(
            start - max_sep,
            end + max_sep - 1,
            start - max_sep,
            end + max_sep - 1,
        )

        distances = [0] * 65
        support = {str(p): 0 for p in low_bins}
        signal = {str(p): 0 for p in low_bins}

        seen = set()
        retained_pairs = 0

        for record in records:
            a = int(record.binX)
            b = int(record.binY)
            value = float(record.counts)

            if a > b:
                a, b = b, a

            separation = b - a
            midpoint_twice = a + b

            if not (
                500 <= separation <= max_sep
                and 2 * start <= midpoint_twice < 2 * end
            ):
                continue

            if (
                not math.isfinite(value)
                or value < 0
                or not value.is_integer()
            ):
                raise RuntimeError(
                    "Invalid unnormalized contact count"
                )

            key = (a, b)

            if key in seen:
                raise RuntimeError(
                    "Duplicate upper-triangle contact"
                )

            seen.add(key)
            distance_bin = separation // 500
            distances[distance_bin] += int(value)
            retained_pairs += 1

            for position in (a, b):
                label = str(position)

                if label in support:
                    support[label] += 1
                    signal[label] += int(value)

        total = sum(distances[1:])

        if total <= 0:
            raise RuntimeError("No off-diagonal signal")

        cumulative = 0
        w80 = None

        for distance_bin in range(1, 65):
            cumulative += distances[distance_bin]

            if cumulative >= 0.8 * total:
                w80 = distance_bin
                break

        result.update({
            "status": "SUCCESS",
            "retained_sparse_pairs": retained_pairs,
            "total_contact_signal": total,
            "w80_bins": w80,
            "w80_bp": w80 * 500,
            "low_bin_observed_pairs": support,
            "low_bin_signal": signal,
            "signal_by_band": {
                "0.5-5kb": sum(distances[1:11]),
                "5.5-16kb": sum(distances[11:33]),
                "16.5-32kb": sum(distances[33:65]),
            },
        })

except Exception as exc:
    result.update({
        "status": "QUERY_ERROR",
        "error": f"{type(exc).__name__}: {exc}",
    })

print(json.dumps(result, sort_keys=True))
'''


def fetch_metadata(accession):
    """Previously verified ENCODE metadata, Steps 114-115.

    Cached to avoid repeated ENCODE API requests.
    Contact file access is tested independently.
    """
    parents = [
        "/files/ENCFF583BZD/",
        "/files/ENCFF624NZV/",
    ]

    cached = {
        "ENCFF573OPJ": {
            "accession": "ENCFF573OPJ",
            "assembly": "GRCh38",
            "output_type": "mapping quality thresholded contact matrix",
            "biological_replicates": [1, 2],
            "derived_from": parents,
            "aliases": [
                "encode-processing-pipeline:inter_30"
            ],
            "cloud_metadata": {},
        },
        "ENCFF726YRX": {
            "accession": "ENCFF726YRX",
            "assembly": "GRCh38",
            "output_type": "contact matrix",
            "biological_replicates": [1, 2],
            "derived_from": parents,
            "aliases": [
                "encode-processing-pipeline:inter_1"
            ],
            "cloud_metadata": {},
        },
    }

    if accession not in cached:
        raise ValueError(f"Unverified accession: {accession}")

    return cached[accession]


def byte_range_probe(url):
    request = Request(
        url,
        headers={
            "Range": "bytes=0-4095",
            "User-Agent": "HiC-Diagonal-Research/1.0",
        },
    )

    with urlopen(request, timeout=60) as response:
        status = response.status
        content_range = response.headers.get(
            "Content-Range"
        )
        payload = response.read(4097)

    # Refuse endpoints that ignore HTTP Range requests.
    if status != 206 or len(payload) != 4096:
        raise RuntimeError(
            f"HTTP byte range unsupported: "
            f"status={status}, bytes={len(payload)}"
        )

    if not content_range or not content_range.startswith(
        "bytes 0-4095/"
    ):
        raise RuntimeError(
            "Unexpected Content-Range response"
        )

    return {
        "status": "PASS",
        "http_status": status,
        "returned_bytes": len(payload),
        "content_range": content_range,
    }


@app.function(
    image=image,
    cpu=1,
    memory=2048,
    timeout=600,
)
def compare_matrices():
    results = []

    for accession in ACCESSIONS:
        item = {
            "accession": accession,
            "status": "STARTED",
        }

        try:
            metadata = fetch_metadata(accession)

            cloud = metadata.get("cloud_metadata") or {}

            if not isinstance(cloud, dict):
                cloud = {}

            cloud_url = cloud.get("url")

            if accession == "ENCFF573OPJ":
                candidate_urls = [
                    CURRENT_URL,
                    cloud_url,
                ]
            else:
                candidate_urls = [cloud_url]

            # NCBI GEO lists this exact alternative matrix
            # under experiment GSE237721. Probe HTTP Range
            # support before attempting any matrix reads.
            if accession == "ENCFF726YRX":
                candidate_urls.append(
                    "https://ftp.ncbi.nlm.nih.gov/geo/"
                    "series/GSE237nnn/GSE237721/suppl/"
                    "GSE237721_ENCFF726YRX_"
                    "contact_matrix_GRCh38.hic"
                )

            # Fallback to ENCODE's documented direct
            # file download endpoint if necessary.
            candidate_urls.append(
                f"https://www.encodeproject.org/files/"
                f"{accession}/@@download/{accession}.hic"
            )

            selected_url = None
            range_errors = []

            for url in candidate_urls:
                if not url:
                    continue

                try:
                    range_info = byte_range_probe(url)
                    selected_url = url
                    item["range_probe"] = range_info
                    break

                except Exception as exc:
                    range_errors.append(
                        f"{type(exc).__name__}: {exc}"
                    )

            item["metadata"] = {
                "output_type": metadata.get("output_type"),
                "assembly": metadata.get("assembly"),
                "biological_replicates": metadata.get(
                    "biological_replicates"
                ),
                "derived_from": metadata.get("derived_from"),
                "aliases": metadata.get("aliases"),
            }

            if selected_url is None:
                item["status"] = "NO_WORKING_RANGE_URL"
                item["range_errors"] = range_errors
                results.append(item)
                continue

            item["selected_url"] = selected_url

            completed = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    CHILD_CODE,
                    accession,
                    selected_url,
                ],
                capture_output=True,
                text=True,
                timeout=180,
                check=False,
            )

            if completed.returncode != 0:
                item["status"] = "CHILD_FAILED"
                item["returncode"] = completed.returncode
                item["stderr"] = completed.stderr[-500:]
            else:
                lines = [
                    line for line in
                    completed.stdout.splitlines()
                    if line.strip()
                ]

                item["query"] = json.loads(lines[-1])
                item["status"] = item["query"]["status"]

        except Exception as exc:
            item["status"] = "ERROR"
            item["error"] = (
                f"{type(exc).__name__}: {exc}"
            )

        results.append(item)

    original = next(
        item for item in results
        if item["accession"] == "ENCFF573OPJ"
    )

    baseline_pass = (
        original["status"] == "SUCCESS"
        and original["query"]["w80_bins"] == 13
        and original["query"]["total_contact_signal"] == 4695
    )

    return {
        "experiment": "ENCSR477GZK",
        "peak": [PEAK_START, PEAK_END],
        "resolution_bp": 500,
        "max_separation_bp": 32000,
        "normalization": "NONE",
        "original_baseline_validation": baseline_pass,
        "results": results,
        "interpretation": (
            "Technical comparison of mapping filters; "
            "not independent biological replication"
        ),
    }


@app.local_entrypoint()
def main():
    result = compare_matrices.remote()

    print("\nHCT116 MAPPING-FILTER COMPARISON PROBE")
    print("=" * 65)

    print("Experiment:", result["experiment"])
    print("Peak:", result["peak"])
    print(
        "Original baseline validation:",
        result["original_baseline_validation"],
    )

    for item in result["results"]:
        print(
            f"\nFILE: {item['accession']}"
        )
        print("Status:", item["status"])
        print("Metadata:", item.get("metadata"))
        print("Range probe:", item.get("range_probe"))

        if "query" in item:
            print("Contact query:", item["query"])

        if "error" in item:
            print("Error:", item["error"])

        if "range_errors" in item:
            print("Range errors:", item["range_errors"])

        if "stderr" in item:
            print("Child stderr:", item["stderr"])

    print("\nRead-only mapping-filter probe finished")
