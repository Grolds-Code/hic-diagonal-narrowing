"""Modal CPU preflight for HCT116 Stage 2 replicate analysis.

Performs a read-only regression against the *existing* pooled contact
volume and inventories public ENCODE replicate metadata. It deliberately
does not download large BAM/HIC files or claim replicate concordance.

Run from a Modal-authenticated checkout:
    modal run modal_hct116_stage2_replicate_preflight.py
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import modal

app = modal.App("hic-stage2-replicate-preflight")

volume = modal.Volume.from_name("hic-diagonal-narrowing-data")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy==2.5.4")
    .add_local_file(
        "experiments/stage2_encode_replicate_inventory.py",
        "/root/stage2_encode_replicate_inventory.py",
    )
    .add_local_file(
        "experiments/stage2_pooled_preflight_core.py",
        "/root/stage2_pooled_preflight_core.py",
    )
    .add_local_file(
        "research_artifacts/stage1/hct116_stage1_coverage_triage_v1.csv",
        "/root/hct116_stage1_coverage_triage_v1.csv",
    )
)


@app.function(
    image=image,
    volumes={"/research": volume},
    cpu=2,
    memory=4096,
    timeout=600,
)
def stage2_preflight():
    import sys

    sys.path.insert(0, "/root")

    from stage2_encode_replicate_inventory import (
        build_report,
        query_encode,
    )
    from stage2_pooled_preflight_core import (
        SOURCE,
        audit_pooled_source,
    )

    source_path = Path(f"/research/matrices/{SOURCE}.npz")
    metadata_path = Path(f"/research/metadata/{SOURCE}.json")
    triage_path = Path("/root/hct116_stage1_coverage_triage_v1.csv")

    # This reads a known 500-bp slice already extracted from the 40-GB-class
    # ENCODE .hic product. It never changes the mounted Modal Volume.
    pooled = audit_pooled_source(
        source_path,
        metadata_path,
        triage_path,
    )

    # Small public metadata requests only: no BAM, FASTQ or HIC downloads.
    experiment, files, detailed_files = query_encode()
    inventory = build_report(experiment, files, detailed_files)

    if inventory["matrix_file_count"] < 1:
        raise RuntimeError("No ENCODE matrix metadata found")

    by_accession = {
        row["accession"]: row
        for row in inventory["processed_file_metadata"]
    }
    for accession in ("ENCFF624NZV", "ENCFF583BZD"):
        if (
            accession not in by_accession
            or str(by_accession[accession]["format"]).lower() != "bam"
        ):
            raise RuntimeError(
                "Expected separate biological-replicate BAM unavailable: "
                + accession
            )

    return {
        "study": "HCT116 chr10:20-21Mb 500bp Stage 2 preflight",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "compute_backend": "Modal CPU",
        "volume_name": "hic-diagonal-narrowing-data",
        "pooled_stage1_regression": pooled,
        "encode_replicate_inventory": inventory,
        "stage2_replication_status": "not_tested",
        "limitations": [
            "Pooled source is from both biological replicates.",
            "Validated W80 reproduction is not biological replication.",
            "R1/R2 remain exploratory and mappability-sensitive.",
            "No separate biological-replicate contact matrices were built.",
            "No large BAM, FASTQ, pairs or HIC products were downloaded.",
            "The Modal Volume and previous Stage 2 files were not modified.",
        ],
    }


@app.local_entrypoint()
def main(
    output: str = (
        "research_artifacts/stage2/"
        "hct116_stage2_modal_replicate_preflight_v1.json"
    ),
):
    destination = Path(output)
    if destination.exists():
        raise RuntimeError(
            f"Refusing to overwrite any existing Stage 2 output: {destination}"
        )

    result = stage2_preflight.remote()

    # Do not create an empty destination directory until the remote audit
    # succeeds. No result is written on the research Modal Volume.
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    pooled = result["pooled_stage1_regression"]
    inventory = result["encode_replicate_inventory"]

    print("\nSTAGE 2 — MODAL REPLICATE PREFLIGHT")
    print("=" * 70)
    print("Pooled source SHA256:", pooled["source_sha256"])
    print("Stage 1 scan rows independently reproduced:",
          pooled["stage1_rows_reproduced"])
    print("Stage 1 provisional called windows:",
          pooled["provisional_coverage_screen_calls"])
    print("Candidates:", pooled["candidates"])
    print("ENCODE processed products:", inventory["processed_file_count"])
    print("ENCODE Hi-C matrix candidates:", inventory["matrix_file_count"])
    print("Biological replication: NOT TESTED")
    print("Large sequencing/contact files downloaded: NONE")
    print("Saved local JSON:", destination)
    print("Stage 2 Modal preflight: PASS")
