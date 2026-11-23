"""Bounded, read-only HIC v8 header feasibility audit for four HCT116 controls.

Run: uv run python -m experiments.stage2_external_hic_resolution_audit
Does NOT examine contact blocks, occupancy, W80, or biological concordance.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import re
import struct
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from experiments.stage2_external_hic_range_check import SAMPLES, geo_file_url

PREFIX_BYTES = 262144
TOTAL_BYTE_LIMIT = 4 * PREFIX_BYTES
MAX_STRING_BYTES = PREFIX_BYTES
TARGET_END = 21_000_000
DESTINATION = Path("research_artifacts/stage2/hct116_external_hic_resolution_audit_v2.json")


class InvalidHicHeader(ValueError):
    """Invalid or truncated HIC v8 header."""


class HeaderReader:
    def __init__(self, data: bytes):
        self.data, self.offset = data, 0

    def take(self, size: int) -> bytes:
        if size < 0 or self.offset + size > len(self.data):
            raise InvalidHicHeader("Truncated HIC header")
        chunk = self.data[self.offset:self.offset + size]
        self.offset += size
        return chunk

    def i32(self) -> int:
        return struct.unpack("<i", self.take(4))[0]

    def i64(self) -> int:
        return struct.unpack("<q", self.take(8))[0]

    def cstring(self) -> str:
        end = self.data.find(b"\0", self.offset, min(len(self.data), self.offset + MAX_STRING_BYTES))
        if end < 0:
            raise InvalidHicHeader(f"Unterminated header string or prefix budget exceeded at offset {self.offset}")
        raw = self.take(end - self.offset + 1)[:-1]
        try:
            return raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise InvalidHicHeader("Invalid header string encoding") from exc


def parse_hic_v8_header(data: bytes) -> dict:
    r = HeaderReader(data)
    if r.take(4) != b"HIC\0":
        raise InvalidHicHeader("Invalid HIC magic")
    version = r.i32()
    if version != 8:
        raise InvalidHicHeader(f"Unsupported HIC version: {version}")
    if r.i64() <= 0:
        raise InvalidHicHeader("Invalid master index location")
    genome = r.cstring()
    attributes = r.i32()
    if not 0 <= attributes <= 256:
        raise InvalidHicHeader("Invalid attribute count")
    for _ in range(attributes):
        r.cstring()
        r.cstring()
    chromosome_count = r.i32()
    if not 1 <= chromosome_count <= 1000:
        raise InvalidHicHeader("Invalid chromosome count")
    chromosomes = {}
    for _ in range(chromosome_count):
        name, length = r.cstring(), r.i32()
        if not name or length <= 0 or name in chromosomes:
            raise InvalidHicHeader("Invalid or duplicate chromosome")
        chromosomes[name] = length
    bp_count = r.i32()
    if not 1 <= bp_count <= 100:
        raise InvalidHicHeader("Invalid BP resolution count")
    bp_resolutions = [r.i32() for _ in range(bp_count)]
    if any(x <= 0 for x in bp_resolutions) or len(set(bp_resolutions)) != len(bp_resolutions):
        raise InvalidHicHeader("Invalid BP resolutions")
    frag_count = r.i32()
    if not 0 <= frag_count <= 100:
        raise InvalidHicHeader("Invalid fragment resolution count")
    frag_resolutions = [r.i32() for _ in range(frag_count)]
    if any(x <= 0 for x in frag_resolutions):
        raise InvalidHicHeader("Invalid fragment resolutions")
    genome_ok = genome.lower().replace("_", "").replace("-", "") in {"hg38", "grch38"}
    chr10_name = next((c for c in ("10", "chr10") if c in chromosomes), None)
    chr10_length = chromosomes.get(chr10_name) if chr10_name else None
    target_ok = chr10_length is not None and chr10_length >= TARGET_END
    return {
        "hic_version": version,
        "genome_identifier": genome,
        "grch38_identifier_recognized": genome_ok,
        "chromosome_10_name": chr10_name,
        "chromosome_10_length": chr10_length,
        "target_interval": "chr10:20000000-21000000",
        "target_within_chromosome": target_ok,
        "bp_resolutions": bp_resolutions,
        "bp_500_present": 500 in bp_resolutions,
        "header_bytes_parsed": r.offset,
        "metadata_gate_passed": bool(genome_ok and target_ok and 500 in bp_resolutions),
        "local_500bp_contacts_checked": False,
        "local_500bp_coverage_sufficient": None,
    }


def validate_range_response(response, limit: int = PREFIX_BYTES) -> tuple[bytes, int]:
    if getattr(response, "status", None) != 206:
        raise ValueError("Range not honored; response body not read")
    content_range = response.headers.get("Content-Range", "")
    match = re.fullmatch(r"bytes 0-(\d+)/(\d+)", content_range)
    if match is None:
        raise ValueError("Invalid Content-Range; response body not read")
    end, total = map(int, match.groups())
    if end != limit - 1 or total <= end:
        raise ValueError("Unexpected Content-Range; response body not read")
    raw = response.read(limit + 1)
    if len(raw) != limit:
        raise ValueError("Short or oversized range response")
    return raw, total


def audit_one(accession: str, filename: str, *, opener=urlopen) -> dict:
    request = Request(geo_file_url(accession, filename), headers={
        "Range": f"bytes=0-{PREFIX_BYTES - 1}",
        "Accept-Encoding": "identity",
        "User-Agent": "hic-diagonal-narrowing-resolution-audit/1.0",
    })
    with opener(request, timeout=45) as response:
        payload, total = validate_range_response(response)
    result = parse_hic_v8_header(payload)
    result.update(status="metadata_checked", range_bytes_read=len(payload),
                  remote_file_size_bytes=total)
    return result


def build_report(*, opener=urlopen) -> dict:
    samples = []
    for accession, filename in SAMPLES.items():
        try:
            result = audit_one(accession, filename, opener=opener)
        except (HTTPError, URLError, TimeoutError, OSError, RuntimeError, ValueError) as exc:
            result = {
                "status": "metadata_unverified",
                "error_class": type(exc).__name__,
                "error_description": str(exc)[:200],
                "bp_500_present": None,
                "metadata_gate_passed": False,
                "local_500bp_contacts_checked": False,
                "local_500bp_coverage_sufficient": None,
            }
        samples.append({"sample": accession, "matrix_filename": filename,
                        "url": geo_file_url(accession, filename), "result": result})
    return {
        "study": "GSE310051 HCT116 WT untreated 0h",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "max_bytes_per_file": PREFIX_BYTES,
        "maximum_total_payload_bytes": TOTAL_BYTE_LIMIT,\n        "revision": "v2: expanded bounded prefix for potentially long Juicer metadata strings",
        "full_hic_files_downloaded": False,
        "contact_data_analyzed": False,
        "biological_replication_tested": False,
        "all_four_pass_metadata_gate": all(s["result"]["metadata_gate_passed"] for s in samples),
        "decision": "Metadata only; local raw contact coverage and W80 remain untested.",
        "samples": samples,
    }


def main() -> None:
    if DESTINATION.exists():
        raise SystemExit(f"Refusing to overwrite research artifact: {DESTINATION}")
    report = build_report()
    DESTINATION.parent.mkdir(parents=True, exist_ok=True)
    with DESTINATION.open("x", encoding="utf-8") as out:
        json.dump(report, out, indent=2, sort_keys=True)
        out.write("\n")
    for row in report["samples"]:
        r = row["result"]
        print(row["sample"], r["status"], "500bp:", r.get("bp_500_present"),
              "genome:", r.get("genome_identifier", "unverified"))
    print("Metadata gate all four:", report["all_four_pass_metadata_gate"])
    print("Local contacts/occupancy: NOT CHECKED")
    print("Saved:", DESTINATION)


if __name__ == "__main__":
    main()
