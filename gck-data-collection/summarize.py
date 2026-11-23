import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from collect import ALLOSTERIC, DEFAULT_OUTPUT_DIR, csv_rows, load_json, save_json


LIGAND_CLASSES = {
    "GLC": "substrate",
    "BGC": "substrate",
    "AGS": "ATP analogue",
    "ANP": "ATP analogue",
    "F6P": "GKRP effector",
    "F6R": "GKRP effector",
    "MG": "cofactor ion",
    "DMS": "solvent",
    "EDO": "cryoprotectant/buffer",
    "GOL": "cryoprotectant/buffer",
    "IOD": "ion",
    "K": "ion",
    "NA": "ion",
    "SO4": "ion",
}
POTENCY_TYPES = {"EC50", "AC50", "IC50", "Ki", "Kd"}
ACTIVATION_TYPES = {"EC50", "AC50"}


def names(values):
    return ";".join(sorted(set(values)))


def numeric_potency(row):
    return (row["standard_type"] in POTENCY_TYPES and row["standard_value"] is not None
            and row["standard_units"] == "nM" and float(row["standard_value"]) > 0)


def assay_class(row):
    assay = row["assay_description"].lower()
    if row["standard_type"] in ACTIVATION_TYPES and any(k in assay for k in ("activat", "induction")):
        return "activation"
    if assay.startswith(("inhibition of recombinant human pancreatic glucokinase by g-6-p",
                         "inhibition of recombinant human n-terminal gst-fused glucokinase",
                         "inhibition of hk4 (unknown origin)")):
        return "direct enzyme inhibition assay"
    return ""


def main(output_dir):
    candidates = load_json(output_dir / "pdb_candidates.json")
    human = {id_.split("_")[0] for id_ in candidates["hits"]["human_uniprot"]}
    polymers = {
        id_.split("_")[0]: candidates["polymers"][id_]
        for id_ in candidates["hits"]["human_uniprot"] + candidates["hits"]["xenopus_uniprot"] + ["9XYT_1"]
    }
    pdb = load_json(output_dir / "pdb_raw.json")
    chembl = load_json(output_dir / "chembl_activities_raw.json")
    structures = load_json(output_dir / "structure_matches.json")

    ligands_by_entry = defaultdict(set)
    entries_by_ligand = defaultdict(set)
    for eid, ligand in pdb["ligand_entities"].items():
        pdb_id = eid.split("_")[0]
        comp_id = ligand["pdbx_entity_nonpoly"]["comp_id"]
        ligands_by_entry[pdb_id].add(comp_id)
        entries_by_ligand[comp_id].add(pdb_id)
    assert set(entries_by_ligand) == ALLOSTERIC | LIGAND_CLASSES.keys()
    assert set(polymers) == set(pdb["entries"])

    bindingdb = defaultdict(list)
    entry_rows = []
    for pdb_id, entry in sorted(pdb["entries"].items()):
        polymer = polymers[pdb_id]
        entity = polymer["rcsb_polymer_entity"]
        entry_info = entry["rcsb_entry_info"]
        ligands = sorted(ligands_by_entry[pdb_id])
        annotations = entry.get("rcsb_binding_affinity", [])
        for aff in annotations:
            assert aff["comp_id"] in ligands
            if aff["comp_id"] in ALLOSTERIC:
                bindingdb[aff["comp_id"]].append({"pdb_id": pdb_id, **aff})
        entry_rows.append({
            "pdb_id": pdb_id,
            "organism": names(x["scientific_name"] for x in polymer.get("rcsb_entity_source_organism", [])),
            "gck_uniprot": names(polymer["rcsb_polymer_entity_container_identifiers"].get("uniprot_ids", [])),
            "gck_entity": polymer["rcsb_id"],
            "gck_description": entity["pdbx_description"],
            "gck_fragment": entity.get("pdbx_fragment", ""),
            "gck_depositor_mutations": entity.get("pdbx_mutation", ""),
            "method": names(exp["method"] for exp in entry["exptl"]),
            "resolution_angstrom": names(map(str, entry_info.get("resolution_combined", []))),
            "deposited_models": entry_info["deposited_model_count"],
            "gck_polymer_copies": entity.get("pdbx_number_of_molecules", ""),
            "title": entry["struct"]["title"],
            "doi": entry.get("rcsb_primary_citation", {}).get("pdbx_database_id_DOI", ""),
            "nonpolymer_ligands": json.dumps([
                {"comp_id": l, "name": " ".join(pdb["ligands"][l]["chem_comp"]["name"].split()),
                 "class": "allosteric-site drug-like ligand" if l in ALLOSTERIC else LIGAND_CLASSES[l]}
                for l in ligands
            ]),
            "allosteric_ligand_ids": names(set(ligands) & ALLOSTERIC),
            "rcsb_binding_affinity": json.dumps(annotations),
        })
    csv_rows(output_dir / "pdb_gck_entries.csv", list(entry_rows[0]), entry_rows)

    potency_by_molecule = defaultdict(list)
    for act in chembl:
        if numeric_potency(act):
            potency_by_molecule[act["molecule_chembl_id"]].append(act)

    allosteric_rows, join_rows = [], []
    for ligand in sorted(ALLOSTERIC):
        comp = pdb["ligands"][ligand]
        desc = comp["rcsb_chem_comp_descriptor"]
        eids = entries_by_ligand[ligand]
        human_eids = eids & human
        matches = structures[ligand]
        exact = {m["molecule_chembl_id"] for m in matches["chembl"]["molecules"]}
        uni = {s["compoundId"] for comp_match in matches["unichem"].get("compounds", [])
               for s in comp_match.get("sources", []) if s["id"] == 1}
        assert exact == uni, (ligand, exact, uni)
        activity = [row for id_ in exact for row in potency_by_molecule[id_]]
        best = min(activity, default={}, key=lambda row: (
            assay_class(row) != "activation", row["standard_relation"] != "=",
            float(row["standard_value"]), row["activity_id"]))
        bdb = bindingdb[ligand]
        ligand_role = ("activator" if ligand != "0H5"
                       else "putative inhibitor (PDB deposition; functional assay unverified)")
        row = {
            "comp_id": ligand,
            "name": " ".join(comp["chem_comp"]["name"].split()),
            "ligand_role": ligand_role,
            "smiles": desc["SMILES_stereo"],
            "inchikey": desc["InChIKey"],
            "pdb_ids": names(eids),
            "human_gck_pdb_ids": names(human_eids),
            "human_gck_entry_count": len(human_eids),
        }
        allosteric_rows.append(row)
        connectivity_id = "CHEMBL2346884" if ligand == "1FY" else ""
        join_rows.append({
            **row,
            "chembl_ids_exact_inchikey": names(exact),
            "chembl_ids_unichem": names(uni),
            "chembl_best_type": best.get("standard_type", ""),
            "chembl_best_value": best.get("standard_value", ""),
            "chembl_best_units": best.get("standard_units", ""),
            "chembl_best_relation": best.get("standard_relation", ""),
            "chembl_best_activity_id": best.get("activity_id", ""),
            "chembl_best_assay_id": best.get("assay_chembl_id", ""),
            "chembl_best_document_id": best.get("document_chembl_id", ""),
            "chembl_all_potency_activity_ids": names(str(x["activity_id"]) for x in activity),
            "chembl_connectivity_only_id_excluded_from_strict_join": connectivity_id,
            "chembl_connectivity_only_note": (
                "Racemate, compound 30, EC50 1600 nM, activity 12717683; PDB is the (2S) enantiomer"
                if connectivity_id else ""
            ),
            "bindingdb_pdb_ids": names(a["pdb_id"] for a in bdb),
            "bindingdb_types": names(a["type"] for a in bdb),
            "bindingdb_via_pdb_annotations": json.dumps(bdb),
            "bindingdb_via_pdb_record_count": len(bdb),
            "potency_sources": names(
                source for source, data in (("BindingDB-via-PDB", bdb), ("ChEMBL", activity)) if data),
        })
    csv_rows(output_dir / "pdb_gck_activators.csv", list(allosteric_rows[0]), allosteric_rows)
    csv_rows(output_dir / "pdb_chembl_join.csv", list(join_rows[0]), join_rows)

    def molecule_count(rows):
        return len({r["molecule_chembl_id"] for r in rows})

    all_potencies = [row for row in chembl if numeric_potency(row)]
    joined = [row for row in join_rows if row["potency_sources"]]
    metrics = {
        "pdb_entries": len(entry_rows),
        "pdb_human_entries": len(human),
        "pdb_methods": dict(Counter(r["method"] for r in entry_rows)),
        "pdb_models_per_entry": dict(Counter(r["deposited_models"] for r in entry_rows)),
        "pdb_allosteric_ligands": len(allosteric_rows),
        "pdb_activator_ligands": sum(r["ligand_role"] == "activator" for r in allosteric_rows),
        "pdb_allosteric_entries": len({e for row in allosteric_rows for e in row["pdb_ids"].split(";")}),
        "bindingdb_affinity_ligands": sum(bool(r["bindingdb_types"]) for r in join_rows),
        "bindingdb_entry_annotations": sum(len(rows) for rows in bindingdb.values()),
        "chembl_exact_potency_ligands": sum(bool(r["chembl_best_type"]) for r in join_rows),
        "pdb_union_potency_ligands": len(joined),
        "pdb_union_potency_entries": len({e for row in joined for e in row["pdb_ids"].split(";")}),
        "pdb_union_multiple_entry_ligands": [r["comp_id"] for r in joined if r["human_gck_entry_count"] > 1],
        "pdb_union_activation_ligands": sum(r["chembl_best_type"] in ACTIVATION_TYPES or
                                           bool(ACTIVATION_TYPES & set(r["bindingdb_types"].split(";")))
                                           for r in joined),
        "chembl_rows": len(chembl),
        "chembl_numeric_potency_compounds": molecule_count(all_potencies),
        "chembl_numeric_ec50_ac50_compounds": molecule_count(
            row for row in all_potencies if row["standard_type"] in ACTIVATION_TYPES),
        "chembl_numeric_ic50_ki_compounds": molecule_count(
            row for row in all_potencies if row["standard_type"] in {"IC50", "Ki"}),
        "chembl_numeric_four_endpoint_compounds": molecule_count(
            row for row in all_potencies if row["standard_type"] != "Kd"),
        "chembl_numeric_kd_compounds": molecule_count(
            row for row in all_potencies if row["standard_type"] == "Kd"),
        "chembl_activation_assay_compounds": molecule_count(
            row for row in all_potencies if assay_class(row) == "activation"),
        "chembl_direct_inhibition_assay_compounds": molecule_count(
            row for row in all_potencies if assay_class(row) == "direct enzyme inhibition assay"),
    }
    save_json(output_dir / "part1_metrics.json", metrics)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                        help="directory for snapshots and CSVs (default: %(default)s)")
    args = parser.parse_args()
    main(args.output_dir)
