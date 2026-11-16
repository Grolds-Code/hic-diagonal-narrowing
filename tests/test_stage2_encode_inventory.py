"""Offline tests for the Stage 2 ENCODE inventory; no external requests."""

import unittest

from experiments.stage2_encode_replicate_inventory import (
    accession_from_link,
    build_report,
    classify_matrix,
    file_summary,
)


class Stage2InventoryTests(unittest.TestCase):
    def test_parent_accessions(self):
        self.assertEqual(
            accession_from_link("/files/ENCFF583BZD/"),
            "ENCFF583BZD",
        )
        self.assertEqual(
            accession_from_link(
                {"accession": "ENCFF624NZV"}
            ),
            "ENCFF624NZV",
        )
        self.assertIsNone(accession_from_link(None))

    def test_replicate_labels_do_not_imply_independence(self):
        pooled = file_summary({
            "accession": "ENCFF573OPJ",
            "file_format": "hic",
            "output_type": "mapping quality thresholded contact matrix",
            "biological_replicates": [2, 1],
            "derived_from": [
                "/files/ENCFF583BZD/",
                {"accession": "ENCFF624NZV"},
            ],
        })
        self.assertEqual(
            pooled["derived_from_file_accessions"],
            ["ENCFF583BZD", "ENCFF624NZV"],
        )
        self.assertEqual(
            classify_matrix(pooled),
            "multiple_biological_replicate_labels",
        )
        self.assertNotIn("validated", classify_matrix(pooled))

    def test_report_preserves_distinct_processed_products(self):
        experiment = {
            "status": "released",
            "replicates": [
                {"biological_replicate_number": 1, "technical_replicate_number": 1},
                {"biological_replicate_number": 2, "technical_replicate_number": 1},
            ],
        }
        listed = [
            {
                "accession": "ENCFF573OPJ",
                "file_format": "hic",
                "output_type": "contact matrix",
            },
            {
                "accession": "ENCFF583BZD",
                "file_format": "bam",
                "output_type": "alignments",
            },
            {
                "accession": "ENCFFRAW",
                "file_format": "fastq",
                "output_type": "reads",
            },
        ]
        detailed = {
            "ENCFF573OPJ": {
                **listed[0],
                "biological_replicates": [1, 2],
                "derived_from": ["/files/ENCFF583BZD/"],
            },
            "ENCFF583BZD": {
                **listed[1],
                "biological_replicates": [1],
            },
        }
        report = build_report(experiment, listed, detailed)
        self.assertEqual(report["listed_file_count"], 3)
        self.assertEqual(report["processed_file_count"], 2)
        self.assertEqual(report["matrix_file_count"], 1)
        self.assertEqual(
            report["matrix_files"][0]["biological_replicates"],
            ["1", "2"],
        )
        self.assertEqual(
            report["matrix_files"][0]["replicate_assessment"],
            "multiple_biological_replicate_labels",
        )

    def test_unlabeled_matrix_is_not_declared_replicated(self):
        self.assertEqual(
            classify_matrix({
                "biological_replicates": [],
            }),
            "replicate_label_unavailable",
        )


if __name__ == "__main__":
    unittest.main()
