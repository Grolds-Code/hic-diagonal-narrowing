"""Read-only provenance and sequencing metadata for ENCFF573OPJ."""

import json
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import modal

app = modal.App("hic-diagonal-narrowing-source-metadata")

image = modal.Image.debian_slim(python_version="3.12")

BASE = "https://www.encodeproject.org"

ENDPOINTS = {
    "hic_file": (
        "/files/ENCFF573OPJ/?format=json&frame=object"
    ),
    "experiment": (
        "/experiments/ENCSR477GZK/"
        "?format=json&frame=object"
    ),
    "fastqs": (
        "/search/?type=File"
        "&dataset=%2Fexperiments%2FENCSR477GZK%2F"
        "&file_format=fastq"
        "&format=json&frame=object&limit=all"
    ),
}


def fetch_json(path):
    request = Request(
        BASE + path,
        headers={
            "Accept": "application/json",
            "User-Agent": "HiC-Diagonal-Research/1.0",
        },
    )

    try:
        with urlopen(request, timeout=45) as response:
            content = response.read(20_000_001)

            if len(content) > 20_000_000:
                raise RuntimeError(
                    "Metadata response exceeded 20 MB limit"
                )

            return {
                "status": "SUCCESS",
                "data": json.loads(content),
            }

    except (HTTPError, URLError, ValueError, RuntimeError) as exc:
        return {
            "status": "FAILED",
            "error": f"{type(exc).__name__}: {exc}",
        }


@app.function(image=image, cpu=1, memory=512, timeout=240)
def inspect_source():
    responses = {
        name: fetch_json(path)
        for name, path in ENDPOINTS.items()
    }

    result = {
        "accession": "ENCFF573OPJ",
        "expected_experiment": "ENCSR477GZK",
        "api_status": {
            name: response["status"]
            for name, response in responses.items()
        },
        "errors": {
            name: response.get("error")
            for name, response in responses.items()
            if response["status"] != "SUCCESS"
        },
    }

    file_response = responses["hic_file"]

    if file_response["status"] == "SUCCESS":
        data = file_response["data"]

        result["hic_file"] = {
            field: data.get(field)
            for field in (
                "accession",
                "dataset",
                "file_format",
                "output_type",
                "assembly",
                "status",
                "file_size",
                "md5sum",
                "derived_from",
                "biological_replicates",
                "technical_replicates",
                "analysis_step_version",
            )
        }

    experiment_response = responses["experiment"]

    if experiment_response["status"] == "SUCCESS":
        data = experiment_response["data"]

        result["experiment"] = {
            field: data.get(field)
            for field in (
                "accession",
                "description",
                "assay_title",
                "biosample_summary",
                "replication_type",
                "bio_replicate_count",
                "tech_replicate_count",
                "default_analysis",
            )
        }

    fastq_response = responses["fastqs"]

    if fastq_response["status"] == "SUCCESS":
        data = fastq_response["data"]
        files = data.get("@graph", [])

        fastqs = []

        for item in files:
            if item.get("file_format") != "fastq":
                continue

            fastqs.append({
                field: item.get(field)
                for field in (
                    "accession",
                    "status",
                    "read_length",
                    "run_type",
                    "paired_end",
                    "platform",
                    "replicate",
                    "biological_replicates",
                    "technical_replicates",
                )
            })

        result["fastqs"] = fastqs

        result["read_lengths_reported"] = sorted({
            int(item["read_length"])
            for item in fastqs
            if isinstance(item["read_length"], int)
        })

        result["fastq_records"] = len(fastqs)

        result["search_total"] = data.get("total")

    result["metadata_only"] = True
    result["read_length_verified"] = bool(
        result.get("read_lengths_reported")
    )

    return result


@app.local_entrypoint()
def main():
    result = inspect_source.remote()

    print("\nHCT116 SOURCE SEQUENCING METADATA AUDIT")
    print("=" * 65)

    print("File:", result["accession"])
    print("Expected experiment:", result["expected_experiment"])
    print("Metadata only:", result["metadata_only"])
    print("API statuses:", result["api_status"])

    for name, error in result["errors"].items():
        print(f"API error ({name}): {error}")

    print("\nHIC FILE PROVENANCE")
    for key, value in result.get("hic_file", {}).items():
        print(f"{key}: {value}")

    print("\nEXPERIMENT METADATA")
    for key, value in result.get("experiment", {}).items():
        print(f"{key}: {value}")

    print("\nRAW FASTQ METADATA")
    print("FASTQ records:", result.get("fastq_records"))
    print(
        "Reported read lengths:",
        result.get("read_lengths_reported"),
    )

    for item in result.get("fastqs", []):
        print(json.dumps(item, sort_keys=True))

    print(
        "\nRead length verified:",
        result["read_length_verified"],
    )
    print("Source metadata audit completed")
