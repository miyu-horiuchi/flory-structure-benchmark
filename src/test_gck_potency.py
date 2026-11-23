import csv
import unittest
from pathlib import Path


class GckPotencyTests(unittest.TestCase):
    def test_excludes_misannotated_protein_kinase_measurements(self):
        path = Path(__file__).resolve().parents[1] / "gck-kinetics" / "data" / "gck_potency_all.csv"
        with path.open(encoding="utf-8", newline="") as handle:
            records = {(row["source_type"], row["source_record"]) for row in csv.DictReader(handle)}
        excluded = {
            ("chembl", "18938924"),
            ("chembl", "19147864"),
            ("chembl", "22423052"),
            ("chembl", "22843734"),
            ("chembl", "24885021"),
            ("bindingdb", "51238608"),
            ("bindingdb", "51272807"),
            ("bindingdb", "51355123"),
            ("bindingdb", "51382841"),
            ("bindingdb", "51496798"),
        }
        self.assertFalse(records & excluded)
