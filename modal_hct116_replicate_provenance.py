"""Read-only ENCODE provenance and replicate-file discovery."""

import json
from urllib.request import Request, urlopen

import modal

app = modal.App("hic-diagonal-narrowing-replicate-provenance")

image = modal.Image.debian_slim(python_version="3.12")

BASE = "https://www.encodeproject.org"

ENDPOINTS = {
    "parent_1": (
        "/files/ENCFF583BZD/?format=json&frame=object"
    ),
    "parent_2": (
        "/files/ENCFF624NZV/?format=json&frame=object"
    ),
    "experiment_hic_files": (
        "/search/?type=File"
        "&dataset=%2Fexperiments%2FENCSR477GZK%2F"
        "&file_format=hic"
        "&format=json&frame=object&limit=all"
    ),
}

FIELDS = (
    "accession",
    "status",
    "file_format",
    "output_type",
    "assembly",
    "file_size",
    "biological_replicates",
    "technical_replicates",
    "derived_from",
    "analysis_step_version",
    "preferred_default",
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
            content = response.read(10_000_001)

        if len(content) > 10_000_000:
            raise RuntimeError("Response exceeded 10 MB")

        return {
            "status": "SUCCESS",
            "data": json.loads(content),
        }

    except Exception as exc:
        return {
            "status": "FAILED",
            "error": f"{type(exc).__name__}: {exc}",
        }


def summarize_file(record):
    return {
        field: record.get(field)
        for field in FIELDS
    }


@app.function(
    image=image,
    cpu=1,
    memory=512,
    timeout=300,
)
def inspect_provenance():
    responses = {
        name: fetch(path)
        for name, path in ENDPOINTS.items()
    }

    result = {
        "experiment": "ENCSR477GZK",
        "pooled_hic_file": "ENCFF573OPJ",
        "metadata_only": True,
        "parents": {},
        "hic_files": [],
        "errors": {},
    }

    for name in ("parent_1", "parent_2"):
        response = responses[name]

        if response["status"] == "SUCCESS":
            result["parents"][name] = summarize_file(
                response["data"]
            )
        else:
            result["errors"][name] = response["error"]

    search = responses["experiment_hic_files"]

    if search["status"] == "SUCCESS":
        data = search["data"]

        result["reported_hic_file_count"] = data.get(
            "total"
        )

        for item in data.get("@graph", []):
            if item.get("file_format") != "hic":
                continue

            result["hic_files"].append(
                summarize_file(item)
            )

        result["hic_files"].sort(
            key=lambda item: item["accession"] or ""
        )

    else:
        result["errors"]["hic_search"] = search["error"]

    result["replicate_specific_hic_candidates"] = [
        item for item in result["hic_files"]
        if item.get("status") == "released"
        and len(
            item.get("biological_replicates") or []
        ) == 1
    ]

    return result


@app.local_entrypoint()
def main():
    result = inspect_provenance.remote()

    print("\nHCT116 REPLICATE PROVENANCE AUDIT")
    print("=" * 65)

    print("Experiment:", result["experiment"])
    print("Pooled matrix:", result["pooled_hic_file"])
    print("Metadata only:", result["metadata_only"])

    print("\nIMMEDIATE SOURCE FILES")

    for name, parent in result["parents"].items():
        print(f"\n{name}")

        for field, value in parent.items():
            print(f"  {field}: {value}")

    print("\nEXPERIMENT HIC FILES")
    print(
        "Reported count:",
        result.get("reported_hic_file_count"),
    )

    for item in result["hic_files"]:
        print(
            json.dumps(item, sort_keys=True)
        )

    print("\nSINGLE-REPLICATE HIC CANDIDATES")

    candidates = result[
        "replicate_specific_hic_candidates"
    ]

    print("Count:", len(candidates))

    for item in candidates:
        print(
            "  File:",
            item["accession"],
            "| Biological replicates:",
            item["biological_replicates"],
            "| Output:",
            item["output_type"],
        )

    print("\nQUERY ERRORS")
    print(result["errors"])

    print("\nProvenance audit completed")
