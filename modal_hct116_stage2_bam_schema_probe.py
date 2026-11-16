"""Modal CPU-only bounded probe of replicate BAM headers, not whole BAMs.

Run: modal run modal_hct116_stage2_bam_schema_probe.py

Each of two ENCODE BAMs is read only via an explicitly limited HTTP
byte-range request (up to 2 MiB each); when byte ranges are unsupported,
the probe skips rather than streaming the whole file. This does not
measure region-specific Hi-C contacts or biological replication.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from urllib.error import HTTPError, URLError

import modal

app = modal.App("hic-stage2-bam-schema-probe")
image = (
    modal.Image.debian_slim(python_version="3.12")
    .add_local_file(
        "experiments/stage2_encode_replicate_inventory.py",
        "/root/stage2_encode_replicate_inventory.py",
    )
    .add_local_file(
        "experiments/stage2_bam_bounded_schema.py",
        "/root/stage2_bam_bounded_schema.py",
    )
)

REPLICATE_BAMS = {
    "biological_1": "ENCFF624NZV",
    "biological_2": "ENCFF583BZD",
}


@app.function(image=image, cpu=2, memory=4096, timeout=300)
def inspect_bam_headers():
    import sys

    sys.path.insert(0, "/root")
    from stage2_encode_replicate_inventory import fetch_json
    from stage2_bam_bounded_schema import MAX_BYTES, bounded_download

    output = {
        "assay": "HCT116 ENCODE ENCSR477GZK",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_type": "bounded_bam_schema_probe",
        "per_bam_max_download_bytes": MAX_BYTES,
        "full_bam_downloaded": False,
        "biological_replication_status": "not_tested",
        "files": {},
    }
    for label, accession in REPLICATE_BAMS.items():
        metadata = fetch_json(
            f"https://www.encodeproject.org/files/{accession}/?format=json"
        )
        expected = label.removeprefix("biological_")
        labels = {str(v) for v in metadata.get("biological_replicates", [])}
        summary = {
            "accession": accession,
            "file_format": metadata.get("file_format"),
            "assembly": metadata.get("assembly"),
            "file_size_bytes": metadata.get("file_size"),
            "file_status": metadata.get("status"),
            "bio_replicates": sorted(labels),
            "technical_replicates": metadata.get("technical_replicates"),
            "output_type": metadata.get("output_type"),
            "alignment_and_contact_semantics": "not_validated",
        }
        if (
            metadata.get("file_format") != "bam"
            or metadata.get("assembly") != "GRCh38"
            or labels != {expected}
        ):
            summary["probe"] = {
                "status": "metadata_mismatch_no_download",
            }
            output["files"][label] = summary
            continue
        try:
            summary["probe"] = bounded_download(
                accession, metadata.get("href")
            )
        except (HTTPError, URLError, TimeoutError, ValueError, RuntimeError) as exc:
            summary["probe"] = {
                "status": "bounded_probe_error",
                "error_class": type(exc).__name__,
                "error_message": str(exc)[:300],
            }
        output["files"][label] = summary

    output["limitations"] = [
        "A beginning-of-file alignment sample is not representative.",
        "BAM header/flags alone do not validate Hi-C contact-pair semantics.",
        "No genomic interval or Hi-C replicate W80 was computed.",
        "No BAM or .hic was downloaded in full.",
        "Independent matrices and matched controls remain necessary.",
    ]
    return output


@app.local_entrypoint()
def main(
    output: str = (
        "research_artifacts/stage2/"
        "hct116_stage2_bam_bounded_schema_v1.json"
    ),
):
    destination = Path(output)
    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite: {destination}")

    report = inspect_bam_headers.remote()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nSTAGE 2 — BOUNDED BAM SCHEMA PROBE")
    print("=" * 60)
    for label, item in report["files"].items():
        probe = item.get("probe", {})
        details = probe.get("parsed", {})
        print(
            f"{label} ({item['accession']}): {probe.get('status')}; "
            f"BAM header references={details.get('reference_count')}; "
            f"sampled alignments={details.get('complete_alignment_records_sampled')}"
        )
        print("  Initial sample flags:",
              details.get("flag_counts_in_initial_file_prefix_only"))
    print("Full BAM files downloaded: NO")
    print("Hi-C replicate W80 analysis: NOT PERFORMED")
    print("Saved:", destination)
