"""Collect the full ChEMBL activity dump and the exact BindingDB UniProt subset."""

import argparse
import csv
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import io
import json
from pathlib import Path
from urllib.parse import urljoin
import zipfile

import requests


DATA_DIR = Path(__file__).resolve().parent.parent / "data"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/"
BDB = "https://www.bindingdb.org"
BDB_ARCHIVE = "BindingDB_All_202611_tsv.zip"
RATES = ("kon (M-1-s-1)", "koff (s-1)")


def get(url, **params):
    response = requests.get(url, params=params, timeout=90)
    response.raise_for_status()
    return response.json()


def write_csv(name, records, fields):
    with (DATA_DIR / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def chembl():
    target = get(CHEMBL + "target/CHEMBL3820.json")
    assert target["organism"] == "Homo sapiens"
    assert {x["accession"] for x in target["target_components"]} == {"P35557"}
    activities = []
    page = get(CHEMBL + "activity.json", target_chembl_id="CHEMBL3820",
               limit=1000, order_by="activity_id")
    total = page["page_meta"]["total_count"]
    while True:
        activities.extend(page["activities"])
        if not page["page_meta"]["next"]:
            break
        page = get(urljoin(CHEMBL, page["page_meta"]["next"]))
    assert len(activities) == len({a["activity_id"] for a in activities}) == total
    docs = sorted({a["document_chembl_id"] for a in activities if a["document_chembl_id"]})
    with ThreadPoolExecutor(max_workers=8) as pool:
        metadata = dict(zip(docs, pool.map(lambda doc: get(CHEMBL + f"document/{doc}.json"), docs)))
    for a in activities:
        doc = metadata.get(a["document_chembl_id"], {})
        a.update(doi=doc.get("doi"), pubmed_id=doc.get("pubmed_id"),
                 document_title=doc.get("title"))
    (DATA_DIR / "chembl_gck_activities.json").write_text(
        json.dumps(activities, indent=2) + "\n", encoding="utf-8")
    print(f"ChEMBL: {total} activities; {len(docs)} documents")
    print(Counter(a["standard_type"] for a in activities))


def compounds():
    activities = json.loads((DATA_DIR / "chembl_gck_activities.json").read_text(encoding="utf-8"))
    dois = {"10.1039/C4MD00027G", "10.1021/acs.jmedchem.6b00632",
            "10.1021/acs.jmedchem.5b01367"}
    records = sorted({a["record_id"] for a in activities if a["doi"] in dois})
    with ThreadPoolExecutor(max_workers=8) as pool:
        data = list(pool.map(lambda rec: get(CHEMBL + f"compound_record/{rec}.json"), records))
    (DATA_DIR / "chembl_source_compounds.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"ChEMBL: {len(data)} source compound records")


def bindingdb():
    affinity = get(BDB + "/rest/getLigandsByUniprots", uniprot="P35557", cutoff=1000000)
    (DATA_DIR / "bindingdb_gck_rest.json").write_text(
        json.dumps(affinity, indent=2) + "\n", encoding="utf-8")
    hits = affinity["getLindsByUniprotsResponse"]["affinities"]
    print(f"BindingDB affinity REST: {len(hits)} hits; columns: {sorted({k for x in hits for k in x})}")
    path = DATA_DIR / BDB_ARCHIVE
    if not zipfile.is_zipfile(path):
        partial = path.with_suffix(".part")
        with requests.get(BDB + f"/rwd/bind/downloads/{BDB_ARCHIVE}",
                          stream=True, timeout=120) as response:
            response.raise_for_status()
            with partial.open("wb") as handle:
                for chunk in response.iter_content(2**20):
                    handle.write(chunk)
        partial.replace(path)
    with zipfile.ZipFile(path) as archive:
        member, = [x for x in archive.namelist() if x.endswith(".tsv")]
        with archive.open(member) as binary:
            handle = io.TextIOWrapper(binary, encoding="utf-8", newline="")
            fields = next(csv.reader([next(handle)], delimiter="\t"))
            accession_cols = [i for i, name in enumerate(fields)
                              if name.startswith("UniProt (SwissProt) Primary ID of Target Chain")]
            assert accession_cols and all(r in fields for r in RATES), fields[:30]
            rows = []
            for line in handle:
                if "P35557" not in line:
                    continue
                row, = csv.reader([line], delimiter="\t")
                if any(i < len(row) and row[i].strip() == "P35557" for i in accession_cols):
                    rows.append(dict(zip(fields, row)))
    write_csv("bindingdb_gck_activities.csv", rows, fields)
    print(f"BindingDB TSV: {len(rows)} exact P35557 rows; "
          f"kon={sum(bool(r.get(RATES[0], '').strip()) for r in rows)}, "
          f"koff={sum(bool(r.get(RATES[1], '').strip()) for r in rows)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("database", choices=["chembl", "compounds", "bindingdb"])
    globals()[parser.parse_args().database]()
