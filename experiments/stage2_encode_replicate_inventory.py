"""Inspect ENCODE processed-file availability for independent Hi-C replication.

This is a metadata-only Stage 2 feasibility audit. It never downloads BAM,
FASTQ, or Hi-C data and never equates replicate labels with validated
independent contact matrices. Python standard library only.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EXPERIMENT = "ENCSR477GZK"
MATRIX_FORMATS = {"hic", "cool", "mcool"}
PROCESSED_FORMATS = MATRIX_FORMATS | {"bam", "pairs", "bedpe"}
HOST = "https://www.encodeproject.org"


def accession_from_link(value):
    """Return an ENCODE file accession from an object or file URL."""
    if isinstance(value, dict):
        return value.get("accession") or accession_from_link(value.get("@id"))
    if isinstance(value, str):
        return next(
            (part for part in value.strip("/").split("/") if part.startswith("ENCFF")),
            None,
        )
    return None


def normalize_replicates(value):
    if value is None:
        return []
    if not isinstance(value, list):
        value = [value]
    return sorted({str(v) for v in value if v is not None})


def relevant(item):
    fmt = str(item.get("file_format") or "").lower()
    kind = str(item.get("output_type") or "").lower()
    return fmt in PROCESSED_FORMATS or (
        "contact matrix" in kind or "alignment" in kind
    )


def file_summary(item):
    parents = [
        a
        for value in (item.get("derived_from") or [])
        if (a := accession_from_link(value))
    ]
    return {
        "accession": item.get("accession"),
        "format": item.get("file_format"),
        "output_type": item.get("output_type"),
        "assembly": item.get("assembly"),
        "status": item.get("status"),
        "biological_replicates": normalize_replicates(
            item.get("biological_replicates")
        ),
        "technical_replicates": normalize_replicates(
            item.get("technical_replicates")
        ),
        "derived_from_file_accessions": sorted(set(parents)),
    }


def classify_matrix(item):
    """Triage metadata, *not* biological independence or data quality."""
    bio = item["biological_replicates"]
    if len(bio) == 1:
        return "single_biological_replicate_label_unverified"
    if len(bio) > 1:
        return "multiple_biological_replicate_labels"
    return "replicate_label_unavailable"


def build_report(experiment, listed_files, detailed_files):
    """Pure assembly for unit testing with offline metadata fixtures."""
    snapshots = []
    for row in listed_files:
        accession = row.get("accession")
        if not accession or not relevant(row):
            continue
        full = detailed_files.get(accession, row)
        snap = file_summary(full)
        if not snap["accession"]:
            snap["accession"] = accession
        snapshots.append(snap)

    snapshots.sort(key=lambda x: (str(x["format"]), str(x["accession"])))
    matrices = [
        {**row, "replicate_assessment": classify_matrix(row)}
        for row in snapshots
        if str(row["format"]).lower() in MATRIX_FORMATS
        or "contact matrix" in str(row["output_type"]).lower()
    ]

    replicate_groups = []
    for entry in experiment.get("replicates") or []:
        if not isinstance(entry, dict):
            replicate_groups.append({"reference": str(entry)})
            continue
        replicate_groups.append({
            "reference": entry.get("@id"),
            "biological_replicate_number": entry.get(
                "biological_replicate_number"
            ),
            "technical_replicate_number": entry.get(
                "technical_replicate_number"
            ),
        })

    return {
        "experiment_accession": EXPERIMENT,
        "experiment_status": experiment.get("status"),
        "assembly_expected": "GRCh38",
        "listed_file_count": len(listed_files),
        "listed_formats": dict(
            sorted(Counter(
                str(row.get("file_format") or "unknown")
                for row in listed_files
            ).items())
        ),
        "experiment_replicate_metadata": replicate_groups,
        "processed_file_count": len(snapshots),
        "matrix_file_count": len(matrices),
        "matrix_files": matrices,
        "processed_file_metadata": snapshots,
        "inference_limitations": [
            "A biological_replicates field does not prove independent processing.",
            "A file derived from multiple BAMs may be a pooled product.",
            "No alignments, contact matrices, or sequencing reads were downloaded.",
            "Hi-C quality, replication concordance, and copy-number state were not assessed.",
            "Use this inventory to decide whether replicate-specific contact extraction is feasible.",
        ],
    }


def fetch_json(url, attempts=3):
    req = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "hic-diagonal-narrowing-stage2-metadata-audit/1.0",
        },
    )
    for attempt in range(attempts):
        try:
            with urlopen(req, timeout=45) as response:
                return json.load(response)
        except HTTPError as exc:
            if exc.code not in (429, 500, 502, 503, 504) or attempt == attempts - 1:
                raise
        except (URLError, TimeoutError):
            if attempt == attempts - 1:
                raise
        time.sleep(attempt + 1)
    raise RuntimeError("Unreachable retry guard")


def query_encode():
    experiment = fetch_json(
        f"{HOST}/experiments/{EXPERIMENT}/?format=json"
    )
    query = urlencode({
        "type": "File",
        "dataset": f"/experiments/{EXPERIMENT}/",
        "limit": "all",
        "format": "json",
    })
    listing = fetch_json(f"{HOST}/search/?{query}")
    rows = listing.get("@graph")
    if not isinstance(rows, list):
        raise RuntimeError("ENCODE search did not return a file graph")
    total = listing.get("total")
    if isinstance(total, int) and total != len(rows):
        raise RuntimeError(
            f"Search returned {len(rows)} records but reports {total}"
        )
    targets = sorted({
        row["accession"]
        for row in rows
        if relevant(row) and row.get("accession")
    })
    if len(targets) > 100:
        raise RuntimeError("Refusing more than 100 detailed file requests")

    details = {}
    for accession in targets:
        details[accession] = fetch_json(
            f"{HOST}/files/{accession}/?format=json"
        )
    return experiment, rows, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("research_artifacts/stage2/encode_replicate_inventory_v1.json"),
    )
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite: {args.output}")

    experiment, files, details = query_encode()
    report = build_report(experiment, files, details)
    report["audited_at_utc"] = datetime.now(timezone.utc).isoformat()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("STAGE 2: ENCODE REPLICATE AVAILABILITY (METADATA ONLY)")
    print("Experiment:", EXPERIMENT)
    print("Listed files:", report["listed_file_count"])
    print("Processed files:", report["processed_file_count"])
    print("Matrix candidates:", report["matrix_file_count"])
    print("Formats:", report["listed_formats"])
    for matrix in report["matrix_files"]:
        print(
            f"  {matrix['accession']}: {matrix['format']}; "
            f"bio={matrix['biological_replicates']}; "
            f"input_files={matrix['derived_from_file_accessions']}; "
            f"triage={matrix['replicate_assessment']}"
        )
    print("No claim of independent matrix availability from labels alone.")
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
