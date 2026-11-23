import argparse
import csv
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests


DEFAULT_OUTPUT_DIR = Path("gck-output")
RCSB = "https://data.rcsb.org/rest/v1/core/"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
CHEMBL = "https://www.ebi.ac.uk"


def get(url, **params):
    response = requests.get(url, params=params, timeout=60)
    if response.status_code == 204:
        return {}
    response.raise_for_status()
    return response.json()


def indexed_search(attribute, operator, value):
    query = {
        "query": {"type": "terminal", "service": "text",
                  "parameters": {"attribute": attribute, "operator": operator, "value": value}},
        "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True},
    }
    result = get(SEARCH, json=json.dumps(query))
    hits = {hit["identifier"] for hit in result.get("result_set", [])}
    assert len(hits) == result.get("total_count", 0)
    return hits


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def load_json(path):
    return json.loads(path.read_text())


def csv_rows(path, columns, rows):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parallel_get(url, ids):
    ids = sorted(ids)
    with ThreadPoolExecutor(max_workers=12) as pool:
        return dict(zip(ids, pool.map(lambda id_: get(url + id_.replace("_", "/")), ids)))


def fetch_pdb(output_dir):
    reference = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
    description = "rcsb_polymer_entity.pdbx_description"
    hits = {
        "human_uniprot": sorted(indexed_search(reference, "exact_match", "P35557")),
        "xenopus_uniprot": sorted(indexed_search(reference, "exact_match", "Q91753")),
        "glucokinase": sorted(indexed_search(description, "contains_phrase", "glucokinase")),
        "hexokinase4": sorted(indexed_search(description, "contains_phrase", "hexokinase-4")),
    }
    candidates = set().union(*hits.values())
    polymers = parallel_get(RCSB + "polymer_entity/", candidates)
    save_json(output_dir / "pdb_candidates.json", {"hits": hits, "polymers": polymers})


def fetch_entries(output_dir, ids):
    entries = parallel_get(RCSB + "entry/", {id_.split("_")[0] for id_ in ids})
    ligand_entities = parallel_get(RCSB + "nonpolymer_entity/", {
        f"{entry_id}_{entity_id}"
        for entry_id, entry in entries.items()
        for entity_id in entry["rcsb_entry_container_identifiers"].get("non_polymer_entity_ids", [])
    })
    ligands = parallel_get(RCSB + "chemcomp/", {
        entity["pdbx_entity_nonpoly"]["comp_id"] for entity in ligand_entities.values()
    })
    save_json(output_dir / "pdb_raw.json",
              {"entries": entries, "ligand_entities": ligand_entities, "ligands": ligands})


def fetch_chembl(output_dir):
    path = "/chembl/api/data/activity.json"
    first = get(CHEMBL + path, target_chembl_id="CHEMBL3820", limit=1000, order_by="activity_id")
    total = first["page_meta"]["total_count"]
    page_size = first["page_meta"]["limit"]
    with ThreadPoolExecutor(max_workers=5) as pool:
        pages = pool.map(lambda offset: get(CHEMBL + path, target_chembl_id="CHEMBL3820",
                                            limit=page_size, offset=offset, order_by="activity_id"),
                         range(page_size, total, page_size))
        activities = first["activities"]
        for page in pages:
            activities.extend(page["activities"])
    assert len({a["activity_id"] for a in activities}) == len(activities) == total
    save_json(output_dir / "chembl_activities_raw.json", activities)
    columns = ("activity_id molecule_chembl_id parent_molecule_chembl_id molecule_pref_name canonical_smiles "
               "standard_type standard_value standard_units standard_relation standard_flag pchembl_value "
               "assay_chembl_id assay_type assay_description assay_variant_accession assay_variant_mutation "
               "document_chembl_id document_journal document_year data_validity_comment potential_duplicate "
               "activity_comment target_chembl_id").split()
    csv_rows(output_dir / "chembl_gck_activities.csv", columns, activities)


def fetch_structure_matches(output_dir):
    ligands = load_json(output_dir / "pdb_raw.json")["ligands"]
    ids = sorted(ALLOSTERIC)

    def lookup(id_):
        key = ligands[id_]["rcsb_chem_comp_descriptor"]["InChIKey"]
        result = get(CHEMBL + "/chembl/api/data/molecule.json",
                     molecule_structures__standard_inchi_key=key)
        response = requests.post(CHEMBL + "/unichem/api/v1/compounds",
                                 json={"type": "sourceID", "compound": id_, "sourceID": 3},
                                 timeout=30)
        response.raise_for_status()
        return {"chembl": result, "unichem": response.json()}

    with ThreadPoolExecutor(max_workers=8) as pool:
        save_json(output_dir / "structure_matches.json", dict(zip(ids, pool.map(lookup, ids))))


ALLOSTERIC = set(("0H4 0H5 0H6 1FW 1FX 1FY 1J9 1JD 2N8 3LZ 4DC 8WJ 8WM 926 AJB AJI B84 "
                 "G1S G2T HKG HKM LOI MRK S41 TK1 VO1 VO2").split())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                        help="directory for snapshots and CSVs (default: %(default)s)")
    parser.add_argument("stage", choices=["candidates", "entries", "chembl", "matches"])
    parser.add_argument("ids", nargs="*")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.stage == "candidates":
        fetch_pdb(args.output_dir)
    elif args.stage == "entries":
        hits = load_json(args.output_dir / "pdb_candidates.json")["hits"] if not args.ids else {}
        fetch_entries(args.output_dir, args.ids or hits["human_uniprot"] + hits["xenopus_uniprot"] + ["9XYT_1"])
    elif args.stage == "chembl":
        fetch_chembl(args.output_dir)
    else:
        fetch_structure_matches(args.output_dir)
