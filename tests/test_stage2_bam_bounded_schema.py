"""Pure offline tests of a bounded BAM header/schema probe."""

import struct
import unittest
from unittest.mock import patch
import zlib

from experiments.stage2_bam_bounded_schema import (
    MAX_BYTES,
    bounded_download,
    describe_bam_prefix,
    parse_bgzf_prefix,
)


def _bgzf_member(data):
    """A single valid BGZF gzip member for deterministic tiny tests."""
    z = zlib.compressobj(level=6, wbits=-15)
    compressed = z.compress(data) + z.flush()
    header = (
        b"\x1f\x8b\x08\x04"
        + b"\x00\x00\x00\x00"
        + b"\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<H", 2)
        + b"\x00\x00"
    )
    trailer = struct.pack(
        "<II", zlib.crc32(data) & 0xFFFFFFFF, len(data)
    )
    block = bytearray(header + compressed + trailer)
    struct.pack_into("<H", block, 16, len(block) - 1)
    return bytes(block)


def _bam_fixture():
    header = b"@HD\tVN:1.6\tSO:coordinate\n@PG\tID:probe\n"
    reference = b"chr10\x00"
    bam = (
        b"BAM\x01"
        + struct.pack("<i", len(header))
        + header
        + struct.pack("<i", 1)
        + struct.pack("<i", len(reference))
        + reference
        + struct.pack("<i", 135534747)
    )
    name = b"r1\x00"
    bin_mq_nl = (60 << 8) | len(name)
    flag_nc = (0x41 << 16) | 0
    # BAM fixed alignment core is 32 bytes after block_size.
    core = struct.pack(
        "<iiIIiiii",
        0, 100, bin_mq_nl, flag_nc,
        0, 0, 200, 0,
    )
    return bam + struct.pack("<i", len(core) + len(name)) + core + name


class BoundedBamSchemaTests(unittest.TestCase):
    def test_valid_bgzf_and_bam_header(self):
        block = _bgzf_member(_bam_fixture())
        self.assertEqual(parse_bgzf_prefix(block)[1], 1)
        result = describe_bam_prefix(block)
        self.assertEqual(result["chr10_reference_bp"], 135534747)
        self.assertEqual(result["sam_header_sort_order"], ["coordinate"])
        self.assertEqual(result["reference_count"], 1)
        self.assertEqual(result["complete_alignment_records_sampled"], 1)
        self.assertEqual(
            result["flag_counts_in_initial_file_prefix_only"]["paired"],
            1,
        )
        self.assertEqual(result["contact_pair_semantics"], "not_established")

    def test_truncated_final_bgzf_member_is_discarded(self):
        block = _bgzf_member(_bam_fixture())
        parsed, count = parse_bgzf_prefix(block + block[:7])
        self.assertEqual(count, 1)
        self.assertEqual(parsed, _bam_fixture())

    def test_invalid_bgzf_rejected(self):
        with self.assertRaises(ValueError):
            describe_bam_prefix(b"not a bam")

    def test_missing_range_support_never_reads_response_body(self):
        class FakeResponse:
            status = 200
            def __init__(self):
                self.read_called = False
            def __enter__(self):
                return self
            def __exit__(self, *args):
                return False
            def read(self, *args):
                self.read_called = True
                raise AssertionError("Should never read non-range body")

        response = FakeResponse()
        with patch(
            "experiments.stage2_bam_bounded_schema.urlopen",
            return_value=response,
        ):
            result = bounded_download(
                "ENCFF624NZV",
                "/files/ENCFF624NZV/@@download/ENCFF624NZV.bam",
            )
        self.assertEqual(result["status"], "range_not_supported")
        self.assertFalse(response.read_called)

    def test_unrecognized_download_path_blocks_request(self):
        with patch(
            "experiments.stage2_bam_bounded_schema.urlopen",
        ) as mocked:
            result = bounded_download(
                "ENCFF624NZV",
                "https://other.example.com/large.bam",
            )
        self.assertEqual(
            result["status"], "missing_or_noncanonical_download_href"
        )
        mocked.assert_not_called()

    def test_max_download_budget_is_bounded(self):
        self.assertLessEqual(MAX_BYTES, 2 * 1024 * 1024)


if __name__ == "__main__":
    unittest.main()
