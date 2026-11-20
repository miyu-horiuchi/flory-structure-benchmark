import copy
import json
import unittest
from pathlib import Path

from prepare_data import combine_data


class PersonalComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        root = Path(__file__).resolve().parent
        self.original = json.loads((root / "benchmark.data.json").read_text())
        self.personal = json.loads((root.parent / "personal-af3/results.json").read_text())

    def test_preserves_every_original_field(self) -> None:
        before = copy.deepcopy(self.original)
        result = combine_data(self.original, self.personal)
        self.assertEqual(self.original, before)
        self.assertEqual({k: v for k, v in result.items() if k != "personal_af3"}, before)
        self.assertEqual(len(result["models"]), 7)
        self.assertEqual(sum(r["tier"] == "primary" for r in result["samples"]), 140)

    def test_personal_counts_and_missing_glucose_oxygen(self) -> None:
        result = combine_data(self.original, self.personal)["personal_af3"]
        self.assertEqual(sum(c["clean_hits"] for c in result["cells"]), 13)
        self.assertEqual(sum(c["scored_samples"] for c in result["cells"]), 15)
        for row in (r for r in result["samples"] if r["target"] == "GCK"):
            self.assertEqual(row["evaluation_status"], "evaluation_error")
            self.assertIsNone(row.get("ligand_rmsd_receptor_aligned_a"))
            self.assertEqual(row["ligand_atom_inventory"]["missing_reference_atom_names"], ["O1"])

    def test_retains_aurora_chemical_mismatch(self) -> None:
        for row in (r for r in self.personal["samples"] if r["target"] == "AURKA"):
            chemistry = row["chemical_input_comparison"]
            self.assertEqual(chemistry["server_formal_charge"], 0)
            self.assertEqual(chemistry["original_formal_charge"], 1)
            self.assertFalse(chemistry["canonical_isomeric_smiles_equal"])
            self.assertTrue(row["input_phosphorylation_retained"])

    def test_rejects_missing_or_duplicate_samples(self) -> None:
        for samples in [self.personal["samples"][:-1], [self.personal["samples"][0]] * 20]:
            with self.subTest(samples=len(samples)), self.assertRaises(ValueError):
                combine_data(self.original, {**self.personal, "samples": samples})

    def test_rejects_missing_notice(self) -> None:
        with self.assertRaisesRegex(ValueError, "terms notice"):
            combine_data(self.original, {**self.personal, "notice": ""})

    def test_rejects_changed_ranking_selection(self) -> None:
        self.personal["cells"][0]["selected_sample"] = 4
        with self.assertRaisesRegex(ValueError, "native confidence"):
            combine_data(self.original, self.personal)


if __name__ == "__main__":
    unittest.main()
