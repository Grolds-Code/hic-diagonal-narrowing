"""Cheap, read-only audit of ENCODE HCT116 processed-contact provenance.

Focus on the two .hic files whose replicate fields were absent in the
first inventory; also include the two pooled-source .hic files and
the pooled .pairs output for comparison.

No BAM, FASTQ, pairs or .hic payload is downloaded. Only five small
ENCODE JSON metadata objects are requested. No Modal run is needed.

    uv run python -m experiments.stage2_processed_lineage
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from experiments.stage2_encode_replicate_inventory import (
    accession_from_link,
    fetch_json,
    normalize_replicates,
)

ACCESSIONS = (
    "ENCFF573OPJ",
    "ENCFF726YRX",
    "ENCFF743HCX",
    "ENCFF876NPA",
    "ENCFF498SJT",
)


def reference_name(value):
    if value is None:
        return None
    if isinstance(value, dict):
        return value.get("accession") or value.get("@id") or value.get("uuid")
    if isinstance(value, str):
        return value
    return str(value)


def summarize_metadata(item):
    """Preserve available provenance, but never infer independent replicates."""
    parents = sorted({
        accession for value in (item.get("derived_from") or [])
        if (accession := accession_from_link(value))
    })
    accession = item.get("accession")
    if not isinstance(accession, str) or not accession.startswith("ENCFF"):
        raise ValueError("Missing ENCODE file accession")
    biological = normalize_replicates(item.get("biological_replicates"))
    technical = normalize_replicates(item.get("technical_replicates"))
    return {
        "accession": accession,
        "file_format": item.get("file_format"),
        "output_type": item.get("output_type"),
        "status": item.get("status"),
        "assembly": item.get("assembly"),
        "file_size_bytes": item.get("file_size"),
        "biological_replicates": biological,
        "technical_replicates": technical,
        "derived_from_accessions": parents,
        "analysis_reference": reference_name(item.get("analysis")),
        "step_run_reference": reference_name(item.get("step_run")),
        "dataset_reference": reference_name(item.get("dataset")),
        "pipeline_reference": reference_name(item.get("pipeline")),
        "submitted_file_name": item.get("submitted_file_name"),
        "accession_labels": [str(x) for x in item.get("aliases", [])],
        "independent_biological_contact_map_established": False,
        "assessment": (
            "derived_from_multiple_replicate_BAMs"
            if len(biological) > 1
            else "single_label_provenance_unverified"
            if len(biological) == 1
            else "missing_replicate_label_provenance_unresolved"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "research_artifacts/stage2/"
            "hct116_stage2_processed_lineage_v1.json"
        ),
    )
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(
            f"Refusing to overwrite existing research output: {args.output}"
        )

    summaries = []
    for accession in ACCESSIONS:
        record = fetch_json(
            f"https://www.encodeproject.org/files/{accession}/?format=json"
        )
        result = summarize_metadata(record)
        if result["accession"] != accession:
            raise RuntimeError(f"ENCODE accession mismatch for {accession}")
        summaries.append(result)

    report = {
        "study": "HCT116 ENCSR477GZK processed Hi-C file lineage",
        "source": "https://www.encodeproject.org",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": summaries,
        "data_payloads_downloaded": False,
        "replicate_contact_validation": "not_tested",
        "limitations": [
            "File-level metadata may omit older processing ancestry.",
            "No replicate-specific contact map should be inferred from "
            "a missing derived_from or biological_replicates field.",
            "Unknown sorting of biological replicate BAMs prevents "
            "assuming trivial indexed genomic-region extraction.",
            "No Hi-C genomic contacts or W80 measures were computed.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("STAGE 2 — ENCODE PROCESSED CONTACT FILE LINEAGE")
    print("=" * 65)
    for row in summaries:
        print(
            f"{row['accession']} | {row['file_format']} | "
            f"{row['output_type']} | "
            f"bio={row['biological_replicates']} | "
            f"source={row['derived_from_accessions']}"
        )
        print(
            f"  size={row['file_size_bytes']} "
            f"analysis={row['analysis_reference']} "
            f"step_run={row['step_run_reference']} "
            f"pipeline={row['pipeline_reference']}"
        )
        print("  assessment:", row["assessment"])
    print("Contact payloads downloaded: NO")
    print("Biological replicate W80: NOT TESTED")
    print("Saved:", args.output)
    print("Step 401: METADATA AUDIT COMPLETE")


if __name__ == "__main__":
    main()
