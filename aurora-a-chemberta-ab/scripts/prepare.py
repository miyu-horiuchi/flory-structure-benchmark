import csv
import hashlib
import math
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from common import (
    DATA, GROUPED_LABEL, GROUPED_PLANS, LABELS, RESULTS, ROOT, RUN, prespec_sha, read_json,
    sha256, verify_manifest, write_json,
)


def rows(file):
    with (DATA / file).open(newline="") as handle:
        return list(csv.DictReader(handle))


def source_rows():
    assay = read_json(DATA / "chembl_assay_CHEMBL3991630.json")
    document = read_json(DATA / "chembl_document_CHEMBL3991601.json")
    assert assay["assay_chembl_id"] == "CHEMBL3991630"
    assert assay["assay_type"] == "B" and assay["target_chembl_id"] == "CHEMBL4722"
    assert assay["variant_sequence"] is None and "cell lysate" in assay["description"]
    assert assay["document_chembl_id"] == document["document_chembl_id"] == "CHEMBL3991601"
    assert document["doi"] == "10.1126/science.aan4368"
    activities = [
        row for row in rows("chembl_aurka_exact_match_activities.csv")
        if row["assay_chembl_id"] == assay["assay_chembl_id"] and row["standard_type"] == "Kd"
    ]
    for row in activities:
        assert row["document_chembl_id"] == document["document_chembl_id"]
        assert row["assay_type"] == "B" and row["standard_units"] == "nM"
        assert row["numeric_uncensored"] == "True" and row["standard_relation"] == "="
    assert len(activities) == len({row["inchikey"] for row in activities}) == 15
    return rows("Lake2018_panel.csv"), {row["inchikey"]: row for row in activities}


def pkd(raw):
    value = float(raw)
    assert math.isfinite(value) and value > 0
    return 9 - math.log10(value)


def grouped_folds(scaffolds, label):
    count, size = (3, 5) if label == "primary" else (4, 6)
    folds = [{"id": f"G{i:02}", "scaffold_ids": [], "members": {label: []}} for i in range(1, count + 1)]
    eligible = [scaffold for scaffold in scaffolds if scaffold["members"][label]]
    for scaffold in sorted(eligible, key=lambda s: (-len(s["members"][label]), s["scaffold"])):
        fold = min(folds, key=lambda f: (len(f["members"][label]), f["id"]))
        fold["scaffold_ids"].append(scaffold["id"])
        fold["members"][label].extend(scaffold["members"][label])
    for fold in folds:
        fold["members"][label].sort()
        assert len(fold["members"][label]) == size
    return folds, eligible


def prepare():
    seal = prespec_sha()
    manifest = DATA / "SHA256SUMS"
    source_hashes = verify_manifest(manifest, DATA)
    panel, activities = source_rows()
    labels = {name: {} for name in LABELS}
    compounds = {}
    scaffold_members = defaultdict(list)
    for row in panel:
        name, smiles, key = (row[field] for field in ("compound_name", "canonical_isomeric_smiles", "inchikey"))
        mol = Chem.MolFromSmiles(smiles)
        assert name not in compounds and mol is not None and Chem.MolToInchiKey(mol) == key
        assert Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) == smiles
        scaffold = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol), canonical=True, isomericSmiles=False)
        scaffold_members[scaffold].append(name)
        alpha = float(row["alpha_mean"])
        assert math.isfinite(alpha) and alpha > 0
        compounds[name] = {"smiles": smiles, "inchikey": key, "alpha": alpha, "log10_alpha": math.log10(alpha)}
        labels["coupled"][name] = pkd(row["kd_minus_tpx2_nM"])
        if key in activities:
            measured = activities[key]
            assert measured["compound_name"] == name
            labels["primary"][name] = pkd(measured["standard_value"])

    assert len(compounds) == len({c["inchikey"] for c in compounds.values()}) == 24
    assert len(labels["primary"]) == 15
    folds = []
    for index, (scaffold, members) in enumerate(sorted(scaffold_members.items()), start=1):
        folds.append({"id": f"S{index:02d}", "scaffold": scaffold, "members": {
            label: sorted(set(members) & values.keys()) for label, values in labels.items()
        }})
    assert len(folds) == 23 and sum(bool(f["members"]["primary"]) for f in folds) == 14
    spec = (ROOT / "PRESPEC.md").read_text()
    for name in labels["primary"]:
        assert f"| {name} | {compounds[name]['inchikey']} |" in spec
    for fold in folds:
        members = "; ".join(fold["members"]["coupled"])
        primary = "yes" if fold["members"]["primary"] else "no"
        assert f"| {fold['id']} | {members} | {primary} |" in spec
    data = {
        "prespec_sha256": seal, "source_manifest_sha256": sha256(manifest), "source_sha256": source_hashes,
        "compounds": compounds, "labels": labels, "folds": folds, "dynamics_columns": ["log10_alpha"],
        "geometry_columns": [],
    }
    if GROUPED_LABEL:
        data["folds"], data["scaffolds"] = grouped_folds(folds, GROUPED_LABEL)
        data["labels"] = {GROUPED_LABEL: labels[GROUPED_LABEL]}
        data["exploratory"] = {
            "added_utc": "2026-11-23", "post_hoc": True, "plan_sha256": sha256(GROUPED_PLANS[GROUPED_LABEL]),
            "assignment": "Decreasing scaffold size, then canonical scaffold SMILES; assign intact to least-loaded fold, breaking ties by fold ID.",
        }
    return data


if __name__ == "__main__":
    data = prepare()
    if GROUPED_LABEL:
        RUN.mkdir(parents=True, exist_ok=True)
        RESULTS.mkdir(parents=True, exist_ok=True)
    write_json(RUN / "prepared.json", data)
    for label, values in data["labels"].items():
        print(f"{label}: n={len(values)}; folds={sum(bool(f['members'][label]) for f in data['folds'])}")
    print(f"Prepared from sealed prespec {data['prespec_sha256']}")
