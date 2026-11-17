"""Offline tests for strict no-full-download GEO .hic header preflight."""

import io
import struct
import unittest
from unittest.mock import patch

from experiments.stage2_external_hic_range_check import (
    BYTES_PER_FILE,
    SAMPLES,
    check_range_response,
    extract_hic_header,
    geo_file_url,
    probe,
)


class FakeResponse:
    def __init__(self, status, content_range="", data=b""):
        self.status = status
        self.headers = {"Content-Range": content_range}
        self.body = io.BytesIO(data)
        self.did_read = False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self, n):
        self.did_read = True
        return self.body.read(n)


class ExternalHicPreflightTests(unittest.TestCase):
    def test_geo_paths_match_sample_accessions(self):
        self.assertEqual(len(SAMPLES), 4)
        for accession, filename in SAMPLES.items():
            self.assertIn(
                f"/{accession[:-3]}nnn/{accession}/suppl/",
                geo_file_url(accession, filename),
            )

    def test_invalid_file_cross_sample_rejected(self):
        with self.assertRaises(ValueError):
            geo_file_url("GSM9287304", SAMPLES["GSM9287307"])

    def test_parse_hic_v9(self):
        header = b"HIC\x00" + struct.pack("<i", 9) + b"0" * 12
        self.assertEqual(extract_hic_header(header)["version"], 9)

    def test_invalid_hic_header_rejected(self):
        with self.assertRaises(ValueError):
            extract_hic_header(b"notHICa")

    def test_server_without_range_never_reads_payload(self):
        response = FakeResponse(200, data=b"x" * 256)
        with patch(
            "experiments.stage2_external_hic_range_check.urlopen",
            return_value=response,
        ):
            result = probe(
                geo_file_url("GSM9287304", SAMPLES["GSM9287304"])
            )
        self.assertEqual(result["status"], "range_not_honored")
        self.assertFalse(response.did_read)

    def test_range_success_checks_hic_magic(self):
        header = b"HIC\x00" + struct.pack("<i", 8) + b"x" * 32
        response = FakeResponse(
            206,
            content_range="bytes 0-39/2000000000",
            data=header,
        )
        result = check_range_response(response)
        self.assertEqual(result["status"], "range_header_valid")
        self.assertEqual(result["header"]["version"], 8)

    def test_oversized_content_range_skips_read(self):
        response = FakeResponse(
            206,
            content_range="bytes 0-999999/2000000000",
            data=b"x" * 256,
        )
        result = check_range_response(response)
        self.assertEqual(result["status"], "oversize_content_range")
        self.assertFalse(response.did_read)

    def test_size_budget(self):
        self.assertEqual(BYTES_PER_FILE, 65536)


if __name__ == "__main__":
    unittest.main()
