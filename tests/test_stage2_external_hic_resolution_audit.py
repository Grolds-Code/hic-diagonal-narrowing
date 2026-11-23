"""Offline regression tests for strictly bounded Stage 2 HIC v8 metadata audit."""
import io
import struct
import unittest

from experiments.stage2_external_hic_resolution_audit import (
    InvalidHicHeader, PREFIX_BYTES, SAMPLES, audit_one, build_report,
    parse_hic_v8_header, validate_range_response,
)


def i(value):
    return struct.pack("<i", value)


def header(*, version=8, genome="hg38", resolutions=(1000, 500, 250),
           chromosome_length=133797422):
    return (b"HIC\x00" + i(version) + struct.pack("<q", 90000000)
            + genome.encode() + b"\x00" + i(1) + b"key\x00value\x00"
            + i(2) + b"All\x00" + i(1000000000) + b"10\x00" + i(chromosome_length)
            + i(len(resolutions)) + b"".join(i(x) for x in resolutions) + i(0))


class FakeResponse:
    def __init__(self, status=206, data=None, content_range=None):
        self.status = status
        self.headers = {"Content-Range": content_range or
                        f"bytes 0-{PREFIX_BYTES-1}/2100000000"}
        self.body = io.BytesIO(data if data is not None else
                               header().ljust(PREFIX_BYTES, b"\x00"))
        self.was_read = False

    def read(self, n):
        self.was_read = True
        return self.body.read(n)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class Stage2ResolutionAuditTests(unittest.TestCase):
    def test_500_metadata_pass(self):
        result = parse_hic_v8_header(header())
        self.assertTrue(result["metadata_gate_passed"])
        self.assertIsNone(result["local_500bp_coverage_sufficient"])

    def test_missing_500_fails_gate(self):
        self.assertFalse(parse_hic_v8_header(
            header(resolutions=(1000, 2000)))["metadata_gate_passed"])

    def test_wrong_genome_fails_gate(self):
        self.assertFalse(parse_hic_v8_header(
            header(genome="hg19"))["metadata_gate_passed"])

    def test_short_chromosome_fails_gate(self):
        self.assertFalse(parse_hic_v8_header(
            header(chromosome_length=19000000))["metadata_gate_passed"])

    def test_bad_version_rejected(self):
        with self.assertRaises(InvalidHicHeader):
            parse_hic_v8_header(header(version=9))

    def test_truncated_header_rejected(self):
        with self.assertRaises(InvalidHicHeader):
            parse_hic_v8_header(header()[:30])

    def test_invalid_counts_rejected(self):
        with self.assertRaises(InvalidHicHeader):
            parse_hic_v8_header(b"HIC\0" + i(8) + struct.pack("<q", 12345)
                                + b"hg38\0" + i(999999))

    def test_200_response_body_never_read(self):
        response = FakeResponse(status=200)
        with self.assertRaises(ValueError):
            validate_range_response(response)
        self.assertFalse(response.was_read)

    def test_incorrect_range_never_read(self):
        response = FakeResponse(content_range="bytes 0-999999/2000000000")
        with self.assertRaises(ValueError):
            validate_range_response(response)
        self.assertFalse(response.was_read)

    def test_short_payload_rejected(self):
        with self.assertRaises(ValueError):
            validate_range_response(FakeResponse(data=b"too short"))

    def test_exact_read_limit(self):
        raw, total = validate_range_response(FakeResponse())
        self.assertEqual(len(raw), PREFIX_BYTES)
        self.assertEqual(total, 2100000000)

    def test_audit_one(self):
        accession, filename = next(iter(SAMPLES.items()))
        result = audit_one(accession, filename,
                           opener=lambda *a, **k: FakeResponse())
        self.assertTrue(result["metadata_gate_passed"])
        self.assertEqual(result["range_bytes_read"], PREFIX_BYTES)

    def test_four_samples_no_biological_claim(self):
        report = build_report(opener=lambda *a, **k: FakeResponse())
        self.assertEqual(len(report["samples"]), 4)
        self.assertFalse(report["biological_replication_tested"])
        self.assertFalse(report["contact_data_analyzed"])


if __name__ == "__main__":
    unittest.main()
