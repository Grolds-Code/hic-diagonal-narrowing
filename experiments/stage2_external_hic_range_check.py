"""Bounded HTTP range preflight for four independent HCT116 GEO Hi-C maps.

Never downloads full .hic payloads. Checks only a 64-KiB prefix and the
.hic magic/version, not contact data, resolutions or library quality.

Run from an authenticated local WSL environment, no Modal needed:
    uv run python -m experiments.stage2_external_hic_range_check
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import re
import struct
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

BYTES_PER_FILE = 65536

# Exact filenames are reported by the four individual NCBI GEO sample pages.
# Directory is NCBI's standard GEO sample supplemental-file layout.
SAMPLES = {
    "GSM9287304": "GSM9287304_STRS_HCT116_WT_control_0h_1_1_inter_30.hic",
    "GSM9287305": "GSM9287305_STRS_HCT116_WT_control_0h_2_1_inter_30.hic",
    "GSM9287306": "GSM9287306_STRS_HCT116_WT_control_0h_3_1_inter_30.hic",
    "GSM9287307": "GSM9287307_STRS_HCT116_WT_control_0h_4_1_inter_30.hic",
}


def geo_file_url(accession: str, filename: str) -> str:
    if not re.fullmatch(r"GSM\d{7}", accession):
        raise ValueError("Invalid GEO sample accession")
    if not re.fullmatch(r"GSM\d{7}_[A-Za-z0-9_.-]+\.hic", filename):
        raise ValueError("Invalid GEO HIC filename")
    if not filename.startswith(accession + "_"):
        raise ValueError("GEO sample and filename accession mismatch")
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/"
        f"{accession[:-3]}nnn/{accession}/suppl/{filename}"
    )


def extract_hic_header(first_bytes: bytes) -> dict:
    if len(first_bytes) < 8 or first_bytes[:4] != b"HIC\x00":
        raise ValueError("Unexpected HIC magic/header")
    version = struct.unpack_from("<i", first_bytes, 4)[0]
    if not 1 <= version <= 30:
        raise ValueError("Implausible HIC version")
    return {"magic": "HIC", "version": version}


def check_range_response(response, max_bytes: int = BYTES_PER_FILE) -> dict:
    """Fail safely if server ignores Range; never stream a large 200 body."""
    status = getattr(response, "status", None)
    if status != 206:
        return {
            "status": "range_not_honored",
            "http_status": status,
            "payload_read": False,
        }
    content_range = response.headers.get("Content-Range", "")
    match = re.fullmatch(r"bytes 0-(\d+)/(\d+|\*)", content_range)
    if not match:
        return {
            "status": "invalid_content_range",
            "content_range": content_range,
            "payload_read": False,
        }
    end_byte = int(match.group(1))
    if end_byte >= max_bytes:
        # Compliant server may report exactly bytes 0-65535.
        if end_byte != max_bytes - 1:
            return {"status": "oversize_content_range", "payload_read": False}
    raw = response.read(max_bytes + 1)
    if len(raw) > max_bytes:
        raise RuntimeError("HTTP response exceeds maximum byte allowance")
    header = extract_hic_header(raw)
    return {
        "status": "range_header_valid",
        "http_status": status,
        "content_range": content_range,
        "compressed_bytes_read": len(raw),
        "header": header,
    }


def probe(url: str) -> dict:
    if not url.startswith(
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/"
    ) or not url.endswith(".hic"):
        raise ValueError("Refusing unexpected remote host or file")
    request = Request(
        url,
        headers={
            "Range": f"bytes=0-{BYTES_PER_FILE - 1}",
            "Accept-Encoding": "identity",
            "User-Agent": "hic-narrowing-stage2-geo-header-audit/1.0",
        },
    )
    with urlopen(request, timeout=45) as response:
        return check_range_response(response)


def main():
    destination = Path(
        "research_artifacts/stage2/"
        "hct116_stage2_external_hic_header_preflight_v1.json"
    )
    if destination.exists():
        raise SystemExit(f"Refusing to overwrite: {destination}")
    report = {
        "study": "GSE310051 untreated HCT116 independent controls",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "probe_max_bytes_per_file": BYTES_PER_FILE,
        "full_files_downloaded": False,
        "hic_zoom_resolutions_checked": False,
        "contact_counts_analyzed": False,
        "biological_replication_status": "not_tested",
        "samples": [],
    }
    for accession, filename in SAMPLES.items():
        url = geo_file_url(accession, filename)
        try:
            outcome = probe(url)
        except (HTTPError, URLError, TimeoutError, RuntimeError, ValueError) as exc:
            outcome = {
                "status": "bounded_probe_failed",
                "error_class": type(exc).__name__,
                "error_description": str(exc)[:240],
            }
        record = {
            "sample": accession,
            "matrix_filename": filename,
            "url_constructed_not_yet_independently_verified": url,
            "result": outcome,
        }
        report["samples"].append(record)
        print(accession, outcome["status"], outcome.get("header", {}))

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("Full .hic downloads: NO")
    print("500-bp zoom and local coverage: NOT YET CHECKED")
    print("Saved:", destination)


if __name__ == "__main__":
    main()
