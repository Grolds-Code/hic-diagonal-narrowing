"""Pure tests for Stage 2 ENCODE processed-contact lineage safeguards."""

import unittest

from experiments.stage2_processed_lineage import summarize_metadata


class ProcessedLineageTests(unittest.TestCase):
    def test_missing_replicate_data_never_implies_independence(self):
        result = summarize_metadata({
            "accession": "ENCFF743HCX",
            "file_format": "hic",
            "output_type": "contact matrix",
            "biological_replicates": [],
            "derived_from": [],
        })
        self.assertFalse(
            result["independent_biological_contact_map_established"]
        )
        self.assertEqual(
            result["assessment"],
            "missing_replicate_label_provenance_unresolved",
        )

    def test_known_pooled_source_and_no_unjustified_validation(self):
        result = summarize_metadata({
            "accession": "ENCFF573OPJ",
            "file_format": "hic",
            "output_type": "mapping quality thresholded contact matrix",
            "biological_replicates": [2, 1],
            "derived_from": [
                "/files/ENCFF624NZV/",
                {"accession": "ENCFF583BZD"},
            ],
            "file_size": 100_000_000,
        })
        self.assertEqual(result["biological_replicates"], ["1", "2"])
        self.assertEqual(result["derived_from_accessions"],
                         ["ENCFF583BZD", "ENCFF624NZV"])
        self.assertFalse(
            result["independent_biological_contact_map_established"]
        )

    def test_single_biological_label_still_requires_validation(self):
        result = summarize_metadata({
            "accession": "ENCFFXXXXXX",
            "file_format": "hic",
            "biological_replicates": [1],
        })
        self.assertEqual(result["assessment"],
                         "single_label_provenance_unverified")
        self.assertFalse(
            result["independent_biological_contact_map_established"]
        )

    def test_missing_accession_rejected(self):
        with self.assertRaises(ValueError):
            summarize_metadata({"file_format": "hic"})


if __name__ == "__main__":
    unittest.main()
