"""Check snapshot counts, sample identity boundaries, units and positive joins."""

from decimal import Decimal

from build_tables import ANTOINE, BORZILLERI, MOLLICA, csv_read, json_read


def main():
    rows = csv_read("gck_kinetics_joined.csv")
    found = {row["compound_key"]: row for row in rows}
    assert len(rows) == len(found) == 22
    baseline = [r for r in rows if r["baseline"] == "yes"]
    assert {r["compound_key"] for r in baseline} == {f"borzilleri:{i}" for i in range(1, 18)}
    assert all("50 mM glucose" in r["kinetics_conditions"] for r in baseline)
    for key in ("3", "11"):
        assert found[f"borzilleri:{key}"]["chembl_id"] == "CHEMBL3580738"
    assert found["borzilleri:3"]["koff_s-1"] != found["borzilleri:11"]["koff_s-1"]
    assert found["borzilleri:8"]["smiles"] == found["borzilleri:6"]["smiles"]
    assert not found["borzilleri:8"]["chembl_id"]
    assert found["mollica:3"]["source_doi_or_patent"] == ANTOINE
    assert found["mollica:4"]["potency_type"] == "EC50"
    assert Decimal(found["mollica:4"]["potency_value"]) == 45
    for key in ("mollica:6", "mollica:7"):
        assert not found[key]["smiles"] and "racemate" in found[key]["structure_status"]
    for row in rows:
        assert Decimal(row["potency_value"]) > 0 and row["potency_units"]
        assert Decimal(row["kon_M-1s-1"]) > 0 and Decimal(row["koff_s-1"]) > 0
        assert Decimal(row["residence_time_min"]) > 0
        assert bool(row["KD_SPR_nM"]) == (row["source_doi_or_patent"] == BORZILLERI)
    observations = csv_read("gck_kinetics_observations.csv")
    assert len(observations) == 28 and {r["compound_key"] for r in observations} == set(found)
    assert sum(r["source_doi_or_patent"] == MOLLICA for r in observations) == 7
    potency = csv_read("gck_potency_all.csv")
    assert len(potency) == 2953
    assert not any(r["smiles"] for r in potency if r["compound_id"] in ("mollica:6", "mollica:7"))
    assert len(csv_read("gck_kinetics_candidates.csv")) == 2
    assert len(json_read("chembl_gck_activities.json")) == 3222
    bindingdb = csv_read("bindingdb_gck_activities.csv")
    assert len(bindingdb) == 1484
    assert all(not r["kon (M-1-s-1)"].strip() and not r["koff (s-1)"].strip() for r in bindingdb)
    print("Validated 22 joined samples, 28 kinetic observations, 2 exclusions and 2,953 potency records")


if __name__ == "__main__":
    main()
