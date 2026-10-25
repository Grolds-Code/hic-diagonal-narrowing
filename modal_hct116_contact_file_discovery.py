"""Discover replicate-level contact files in ENCODE metadata."""

import json
from collections import Counter
from urllib.request import Request, urlopen

import modal

app = modal.App("hic-diagonal-narrowing-contact-discovery")
image = modal.Image.debian_slim(python_version="3.12")

BASE = "https://www.encodeproject.org"

ENDPOINTS = {
    "extra_hic_1": "/files/ENCFF743HCX/?format=json&frame=object",
    "extra_hic_2": "/files/ENCFF876NPA/?format=json&frame=object",
    "experiment_files": (
        "/search/?type=File"
        "&dataset=%2Fexperiments%2FENCSR477GZK%2F"
        "&format=json&frame=object&limit=all"
    ),
}

FIELDS = (
    "accession",
    "dataset",
    "file_format",
    "output_type",
    "status",
    "assembly",
    "file_size",
    "biological_replicates",
    "technical_replicates",
    "derived_from",
    "analysis_step_version",
    "description",
    "aliases",
)


def fetch(path):
    request = Request(
        BASE + path,
        headers={
            "Accept": "application/json",
            "User-Agent": "HiC-Diagonal-Research/1.0",
        },
    )

    try:
        with urlopen(request, timeout=60) as response:
            raw = response.read(15_000_001)

        if len(raw) > 15_000_000:
            raise RuntimeError("Metadata response exceeds limit")

        return {"status": "SUCCESS", "data": json.loads(raw)}

    except Exception as exc:
        return {
            "status": "FAILED",
            "error": f"{type(exc).__name__}: {exc}",
        }


def summarize(record):
    return {
        key: record.get(key)
        for key in FIELDS
    }


@app.function(image=image, cpu=1, memory=512, timeout=300)
def inspect_files():
    responses = {
        name: fetch(endpoint)
        for name, endpoint in ENDPOINTS.items()
    }

    result = {
        "experiment": "ENCSR477GZK",
        "metadata_only": True,
        "extra_hic_files": [],
        "file_format_counts": {},
        "contact_candidates": [],
        "single_replicate_candidates": [],
        "errors": {},
    }

    for name in ("extra_hic_1", "extra_hic_2"):
        response = responses[name]

        if response["status"] == "SUCCESS":
            result["extra_hic_files"].append(
                summarize(response["data"])
            )
        else:
            result["errors"][name] = response["error"]

    search = responses["experiment_files"]

    if search["status"] != "SUCCESS":
        result["errors"]["experiment_files"] = search["error"]
        return result

    data = search["data"]
    records = data.get("@graph", [])

    result["search_total"] = data.get("total")
    result["returned_records"] = len(records)

    counts = Counter(
        str(item.get("file_format"))
        for item in records
    )

    result["file_format_counts"] = dict(
        sorted(counts.items())
    )

    for item in records:
        fmt = str(item.get("file_format", "")).lower()
        output = str(item.get("output_type", "")).lower()

        if fmt not in ("hic", "pairs") and output != "pairs":
            continue

        record = summarize(item)
        result["contact_candidates"].append(record)

        replicates = record.get("biological_replicates") or []

        if (
            record.get("status") == "released"
            and len(replicates) == 1
        ):
            result["single_replicate_candidates"].append(
                record
            )

    result["contact_candidates"].sort(
        key=lambda item: item.get("accession") or ""
    )

    return result


@app.local_entrypoint()
def main():
    result = inspect_files.remote()

    print("\nHCT116 CONTACT FILE DISCOVERY")
    print("=" * 65)

    print("Experiment:", result["experiment"])
    print("Metadata only:", result["metadata_only"])

    print("\nADDITIONAL HIC FILE PROVENANCE")
    for item in result["extra_hic_files"]:
        print(json.dumps(item, indent=2))

    print("\nEXPERIMENT FILE INVENTORY")
    print("Search total:", result.get("search_total"))
    print("Records returned:", result.get("returned_records"))
    print("File formats:", result["file_format_counts"])

    print("\nCONTACT FILE CANDIDATES")
    for item in result["contact_candidates"]:
        print(json.dumps(item, sort_keys=True))

    print("\nSINGLE-REPLICATE CONTACT CANDIDATES")
    print("Count:", len(result["single_replicate_candidates"]))

    for item in result["single_replicate_candidates"]:
        print(
            item["accession"],
            "| Format:", item["file_format"],
            "| Replicates:", item["biological_replicates"],
            "| Size:", item["file_size"],
        )

    print("\nErrors:", result["errors"])
    print("\nContact file discovery completed")
