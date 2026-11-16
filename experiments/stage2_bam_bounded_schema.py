"""Bounded Hi-C BAM schema audit. No whole-file or contact extraction.

Decode only complete BGZF members from an explicitly range-limited HTTP
response, inspect the BAM header and a small initial alignment sample.
An initial coordinate-sorted sample is NOT representative of read quality,
contacts, chromatin structure or biological replication.
"""

from __future__ import annotations

from collections import Counter
import gzip
import io
import re
import struct
from urllib.request import Request, urlopen

MAX_BYTES = 2 * 1024 * 1024
MAX_DECOMPRESSED_BYTES = 24 * 1024 * 1024
MAX_ALIGNMENTS = 2000


def parse_bgzf_prefix(raw: bytes) -> tuple[bytes, int]:
    """Return data from complete BGZF members, discarding truncated tail."""
    decoded = bytearray()
    index = 0
    blocks = 0
    while index + 18 <= len(raw):
        if raw[index:index + 3] != b"\\x1f\\x8b\\x08":
            raise ValueError("Invalid BGZF/GZIP member magic")
        if raw[index + 3] & 4 == 0:
            raise ValueError("BGZF FEXTRA flag not set")
        extra_length = struct.unpack_from("<H", raw, index + 10)[0]
        extra_start = index + 12
        extra_end = extra_start + extra_length
        if extra_end > len(raw):
            break
        cursor = extra_start
        member_size = None
        while cursor + 4 <= extra_end:
            slen = struct.unpack_from("<H", raw, cursor + 2)[0]
            if cursor + 4 + slen > extra_end:
                raise ValueError("Invalid BGZF extra subfield")
            if raw[cursor:cursor + 2] == b"BC" and slen == 2:
                member_size = (
                    struct.unpack_from("<H", raw, cursor + 4)[0] + 1
                )
            cursor += slen + 4
        if member_size is None or member_size < 26:
            raise ValueError("Missing BGZF BC member size")
        if index + member_size > len(raw):
            break
        block = raw[index:index + member_size]
        payload = gzip.decompress(block)
        if len(decoded) + len(payload) > MAX_DECOMPRESSED_BYTES:
            raise ValueError("Bounded BGZF decompressed-data limit exceeded")
        decoded.extend(payload)
        index += member_size
        blocks += 1
    if blocks == 0:
        raise ValueError("No complete BGZF member in bounded prefix")
    return bytes(decoded), blocks


def describe_bam_prefix(raw: bytes) -> dict:
    """Parse BAM header and a bounded count of complete alignments."""
    data, blocks = parse_bgzf_prefix(raw)
    if len(data) < 12 or data[:4] != b"BAM\\x01":
        raise ValueError("BGZF prefix is not BAM")
    l_text = struct.unpack_from("<i", data, 4)[0]
    if l_text < 0 or l_text > 2_000_000 or 8 + l_text + 4 > len(data):
        raise ValueError("Truncated or oversized BAM header")
    header = data[8:8 + l_text].decode("utf-8", errors="replace")
    pos = 8 + l_text
    n_ref = struct.unpack_from("<i", data, pos)[0]
    pos += 4
    if n_ref < 1 or n_ref > 20_000:
        raise ValueError("Implausible BAM reference count")
    references = []
    for _ in range(n_ref):
        if pos + 4 > len(data):
            raise ValueError("Incomplete BAM reference list")
        name_size = struct.unpack_from("<i", data, pos)[0]
        pos += 4
        if name_size < 2 or name_size > 2048 or pos + name_size + 4 > len(data):
            raise ValueError("Invalid BAM reference name")
        name = data[pos:pos + name_size]
        if name[-1] != 0:
            raise ValueError("Unterminated BAM reference name")
        pos += name_size
        length_bp = struct.unpack_from("<i", data, pos)[0]
        pos += 4
        references.append((name[:-1].decode("utf-8", "replace"), length_bp))

    flag_counts = Counter()
    mapqs = Counter()
    parsed = 0
    while parsed < MAX_ALIGNMENTS and pos + 36 <= len(data):
        size = struct.unpack_from("<i", data, pos)[0]
        if size < 32 or size > 16 * 1024 * 1024:
            raise ValueError("Invalid BAM alignment block length")
        if pos + 4 + size > len(data):
            break
        refid, refpos, packed, flag_packed = struct.unpack_from(
            "<iiII", data, pos + 4
        )
        flag = (flag_packed >> 16) & 0xFFFF
        mapq = (packed >> 8) & 0xFF
        flag_counts["paired"] += bool(flag & 1)
        flag_counts["proper_pair"] += bool(flag & 2)
        flag_counts["unmapped"] += bool(flag & 4)
        flag_counts["mate_unmapped"] += bool(flag & 8)
        flag_counts["secondary"] += bool(flag & 0x100)
        flag_counts["duplicate"] += bool(flag & 0x400)
        flag_counts["supplementary"] += bool(flag & 0x800)
        flag_counts["chr10_alignment"] += (
            0 <= refid < len(references)
            and references[refid][0] == "chr10"
        )
        mapqs[str(mapq)] += 1
        parsed += 1
        pos += size + 4

    headers = header.splitlines()
    sorting = [
        part[3:] for line in headers if line.startswith("@HD")
        for part in line.split("\t") if part.startswith("SO:")
    ]
    return {
        "bgzf_blocks_decoded": blocks,
        "uncompressed_prefix_bytes": len(data),
        "sam_header_sort_order": sorting,
        "reference_count": len(references),
        "chr10_reference_bp": next(
            (length for name, length in references if name == "chr10"),
            None,
        ),
        "read_group_header_records": sum(
            line.startswith("@RG") for line in headers
        ),
        "program_header_records": sum(
            line.startswith("@PG") for line in headers
        ),
        "complete_alignment_records_sampled": parsed,
        "flag_counts_in_initial_file_prefix_only": dict(flag_counts),
        "mapq_counts_in_initial_file_prefix_only": dict(mapqs),
        "sample_scope": (
            "Only the beginning of the BAM file, not a random sample "
            "or a genomic-region query; not an estimator of global properties"
        ),
        "contact_pair_semantics": "not_established",
    }


def bounded_download(accession: str, href: str) -> dict:
    """Request precisely first MAX_BYTES. Never read an unbounded response."""
    if not re.fullmatch(r"ENCFF[A-Z0-9]{6}", accession):
        raise ValueError("Unexpected file accession")
    required_prefix = f"/files/{accession}/@@download/"
    if not isinstance(href, str) or not href.startswith(required_prefix):
        return {"status": "missing_or_noncanonical_download_href"}

    url = "https://www.encodeproject.org" + href
    req = Request(
        url,
        headers={
            "Range": f"bytes=0-{MAX_BYTES - 1}",
            "Accept-Encoding": "identity",
            "User-Agent": "hic-stage2-bam-bounded-schema-probe/1.0",
        },
    )
    with urlopen(req, timeout=60) as response:
        if response.status != 206:
            return {
                "status": "range_not_supported",
                "http_status": response.status,
            }
        content_range = response.headers.get("Content-Range", "")
        if not re.match(r"^bytes 0-\\d+/\\d+$", content_range):
            return {
                "status": "unrecognized_content_range",
                "content_range": content_range,
            }
        raw = response.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise RuntimeError("Server exceeded strict byte-range budget")
    return {
        "status": "bounded_prefix_obtained",
        "http_status": 206,
        "compressed_bytes": len(raw),
        "parsed": describe_bam_prefix(raw),
    }
