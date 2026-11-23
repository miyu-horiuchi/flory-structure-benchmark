"""Join human GCK measurements by source sample, not only by chemical registry ID."""

import csv
import json
from collections import defaultdict
from decimal import Decimal
import re

from fetch_databases import DATA_DIR, write_csv


RAW_INPUTS = {
    "chembl_gck_activities.json": "chembl",
    "chembl_source_compounds.json": "compounds",
    "bindingdb_gck_activities.csv": "bindingdb",
}
BORZILLERI = "10.1039/C4MD00027G"
MOLLICA = "10.1021/acs.jmedchem.6b00632"
ANTOINE = "10.1021/bi900374c"
GKRP = "10.1021/acs.jmedchem.5b01367"
JOIN = """compound_name smiles potency_type potency_value potency_units potency_assay_conditions
kinetics_type kon_M-1s-1 koff_s-1 residence_time_min KD_SPR_nM kinetics_method kinetics_conditions
source_doi_or_patent source_type extraction_confidence notes compound_key chembl_id potency_source
potency_relation t_half_s baseline structure_status record_id""".split()
POTENCY = """compound_name smiles potency_type potency_value potency_relation potency_units
potency_assay_conditions source_doi_or_patent source_type source_record compound_id notes""".split()
BORZ_POTENCY = ("G6PDH/NADP-coupled activation of recombinant human GCK; matrix of 22 "
                "activator and 16 glucose concentrations; EC50 is half-maximal Km reduction (SI Fig. S2)")
BORZ_SPR = ("Biacore 3000, 25 °C, biotinylated recombinant human GCK (E. coli) captured on "
            "streptavidin; 50 mM glucose, 10 mM HEPES pH 7.4, 150 mM NaCl, 0.5 mM TCEP, "
            "0.05% P20, 1 mg/mL BSA, 1% DMSO; three ligand concentrations in duplicate; global 1:1 fit")
ANT_BUFFER = ("50 mM HEPES/Na, 100 mM NaCl, 5 mM MgCl2, 2 mM TCEP, pH 7.1, 5% DMSO; "
              "human pancreatic GCK isoform 1 residues 12–C terminus, N-terminal His6/TEV")
MOL_NAMES = {"1": "Piragliatin (AKY5)", "2": "TAFMT (MWI119)", "3": "LY2121260 (SGN71)",
             "4": "MWI121", "5": "S-44520", "6": "S-49513 (racemate)", "7": "S-49164 (racemate)"}


def csv_read(name):
    with (DATA_DIR / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def json_read(name):
    return json.loads((DATA_DIR / name).read_text(encoding="utf-8"))


def number(value):
    return Decimal(str(value))


def rate(value, scale):
    return str((number(value) * number(scale)).normalize())


def divide(value, denominator):
    return f"{number(value) / number(denominator):.6g}"


def format_potency(name, smiles, kind, value, units, conditions, source, source_type,
                   record, compound, relation="=", notes=""):
    return dict(zip(POTENCY, (name, smiles, kind, value, relation, units, conditions,
                              source, source_type, record, compound, notes)))


def chem_potency(a):
    return format_potency(a["molecule_pref_name"] or a["molecule_chembl_id"],
                          a["canonical_smiles"], a["standard_type"], a["standard_value"],
                          a["standard_units"], a["assay_description"],
                          a["doi"] or a["document_chembl_id"], "chembl", a["activity_id"],
                          a["molecule_chembl_id"], a["standard_relation"])


def join_potency(row):
    return {key: row[key] for key in POTENCY if key.startswith("potency_")} | {
        "potency_source": f'{row["source_doi_or_patent"]}; {row["source_type"]} {row["source_record"]}'}


def group_activities(activities):
    groups = defaultdict(dict)
    for a in activities:
        if a["doi"] in (BORZILLERI, MOLLICA):
            group = groups[(a["doi"], a["record_id"])]
            assert a["assay_chembl_id"] not in group
            group[a["assay_chembl_id"]] = a
    return groups


def build_kinetics(activities, sources):
    observations = []
    added_potencies = []
    baseline_connectivity = next(a["canonical_smiles"] for a in activities if a["record_id"] == 2481285)
    for (doi, record_id), records in group_activities(activities).items():
        key = sources[record_id]["compound_key"].split(";")[0].split(",")[0]
        is_baseline = doi == BORZILLERI
        if is_baseline:
            p = records["CHEMBL3583541"]
            kd = records["CHEMBL3583540"]
            on = records["CHEMBL3583538"]
            off = records["CHEMBL3583539"]
            half = records["CHEMBL3582523"]
            assert (on["units"], off["units"], half["units"]) == ("10'5/M/s", "10'-3/s", "s")
            name = p["molecule_pref_name"] or f"Borzilleri compound {key}"
            if key == "17":
                name = "TAFMT (Borzilleri compound 17)"
            smiles = p["canonical_smiles"]
            chembl_id = p["molecule_chembl_id"]
            structure = "ChEMBL canonical SMILES"
            notes = "Original Table 1 unavailable; values recovered by ChEMBL compound-record ID; assay conditions from paper SI."
            confidence = "high"
            source_type = "chembl"
            if key in ("6", "15"):
                expected = {"6": ("376000", "0.00657", "105.5", "17.5"),
                            "15": ("379000", "0.0727", "9.7", "189.2")}[key]
                actual = (rate(on["value"], "1e5"), rate(off["value"], "1e-3"),
                          half["value"], kd["standard_value"])
                assert tuple(map(number, expected)) == tuple(map(number, actual))
                source_type = "paper_SI"
                notes += " All SPR values cross-checked in Fig. S1; SI Table S1 also confirms koff."
            elif key == "7":
                assert number(rate(off["value"], "1e-3")) == number("0.00701")
                notes += " Human koff cross-checked in SI Table S1."
            if key in ("3", "11"):
                polarity = "(-)" if key == "3" else "(+)"
                name += f" {polarity}"
                structure = "connectivity only; absolute configuration unverified"
                notes += " ChEMBL merges compounds 3 and 11 under one non-stereo SMILES; source compound names are opposite optical isomers."
                confidence = "med"
            elif key == "5":
                structure = "racemate; connectivity only"
            elif key == "6":
                structure = "racemate; connectivity only"
            elif key == "8":
                smiles = baseline_connectivity
                chembl_id = ""
                structure = "connectivity only; absolute configuration unverified"
                notes += (" Withheld ChEMBL SMILES: its fluorothiazole/cyclopropylsulfone structure "
                          "contradicts the source name and paper SI. Compound 8 is the isolated (-) "
                          "enantiomer of the scaffold of compound 6; only canonical connectivity is shown.")
                confidence = "med"
            potency = chem_potency(p)
            potency["potency_assay_conditions"] = BORZ_POTENCY
            observations.append(dict(
                compound_name=name, smiles=smiles, **join_potency(potency),
                kinetics_type="kon;koff;t_half", **{"kon_M-1s-1": rate(on["value"], "1e5")},
                koff_s_1=rate(off["value"], "1e-3"),
                residence_time_min=divide(1, number(off["value"]) * number("0.06")),
                t_half_s=half["value"], KD_SPR_nM=kd["standard_value"], kinetics_method="SPR",
                kinetics_conditions=BORZ_SPR,
                source_doi_or_patent=doi, source_type=source_type, extraction_confidence=confidence,
                notes=notes + " Residence time calculated as 1/koff; dissociation half-life is listed separately.",
                compound_key=f"borzilleri:{key}", chembl_id=chembl_id,
                record_id=record_id, baseline="yes", structure_status=structure))
            continue

        kd = records.get("CHEMBL4428660") or records["CHEMBL4428661"]
        on = records.get("CHEMBL4428665") or records["CHEMBL4428662"]
        off = records.get("CHEMBL4428667") or records["CHEMBL4428663"]
        dwell = records.get("CHEMBL4428668") or records["CHEMBL4428664"]
        assert (on["units"], off["units"], dwell["units"]) == ("10'6/s/M", "/s", "s")
        if key == "4":
            p, = [a for a in activities if a["molecule_chembl_id"] == kd["molecule_chembl_id"]
                  and a["doi"] == "10.1016/j.bmcl.2009.04.107" and a["standard_type"] == "EC50"]
        else:
            p = kd
        direct = key == "2"
        conditions = ("His6-tagged human pancreatic GCK isoform 1, E. coli BL21(DE3); "
                      "glucose present (concentration not verified)" if direct else
                      "His6-tagged human pancreatic GCK isoform 1, E. coli BL21(DE3); "
                      "TAFMT displacement; glucose concentration not verified")
        potency = chem_potency(p)
        if p is kd:
            potency["potency_assay_conditions"] = conditions
        potency["potency_assay_conditions"] += "; buffer and temperature not verified"
        notes = ("Kinetic values from ChEMBL; compound identifiers, affinities and stereochemical "
                 "scope cross-checked against SI CSV/PDF; conditions not imputed from Antoine et al.")
        is_racemate = key in ("6", "7")
        if is_racemate:
            notes += (" Kinetics and Kd refer to the racemate, not its individual enantiomers. "
                      "Single-isomer ChEMBL SMILES withheld; both configurations are listed in the SI CSV.")
        observations.append(dict(
            compound_name=MOL_NAMES[key], smiles="" if is_racemate else kd["canonical_smiles"],
            **join_potency(potency), kinetics_type="kon;koff;residence_time",
            **{"kon_M-1s-1": rate(on["value"], "1e6")},
            koff_s_1=off["value"], residence_time_min=divide(dwell["value"], 60),
            kinetics_method="stopped-flow fluorometry" + (" (direct binding)" if direct else " (competition)"),
            kinetics_conditions=conditions + "; buffer and temperature not verified",
            source_doi_or_patent=doi, source_type="chembl", extraction_confidence="med",
            notes=notes, record_id=record_id,
            compound_key={"1": "borzilleri:12", "2": "borzilleri:17"}.get(key, f"mollica:{key}"),
            chembl_id=kd["molecule_chembl_id"], baseline="yes" if key in ("1", "2") else "no",
            structure_status="racemate; no single-isomer SMILES" if is_racemate else "ChEMBL canonical SMILES"))

    mappings = {r["compound_key"]: r for r in observations if r["source_doi_or_patent"] == MOLLICA}
    for key, figure, kd, on, off, atp in [
        ("borzilleri:17", "5A; 6", 70, 1130000, "0.43", False),
        ("borzilleri:17", "5A; 8", 40, 1670000, "0.62", True),
        ("mollica:3", "5B; 7", 400, 940000, "0.37", False),
        ("mollica:3", "5B; 9", 300, 1090000, "0.30", True),
    ]:
        identity = mappings[key]
        condition = "100 mM glucose" + (" + 4 mM AMP-PNP" if atp else "") + "; " + ANT_BUFFER
        affinity_method = "equilibrium direct fluorescence" if key == "borzilleri:17" else "equilibrium competition fluorescence"
        potency = format_potency(identity["compound_name"], identity["smiles"], "Kd", kd, "nM",
                                 f"{affinity_method}; room temperature; {condition}", ANTOINE,
                                 "paper_figure", f"Fig. {figure.split(';')[0]}", identity["chembl_id"],
                                 "≈", "Estimated equilibrium Kd from fitted curves, reported in caption/text")
        added_potencies.append(potency)
        measurements = dict(compound_name=identity["compound_name"], smiles=identity["smiles"],
                            compound_key=key, chembl_id=identity["chembl_id"],
                            baseline=identity["baseline"], structure_status=identity["structure_status"])
        observations.append(dict(
            **measurements, **join_potency(potency), kinetics_type="kon;koff",
            **{"kon_M-1s-1": on}, koff_s_1=off, residence_time_min=divide(1, 60 * number(off)),
            kinetics_method="stopped-flow fluorometry; association fit and displacement-derived koff",
            kinetics_conditions=f"25 °C; {condition}", source_doi_or_patent=ANTOINE,
            source_type="paper_figure", extraction_confidence="high",
            record_id=f"Figs. {figure}",
            notes=("Rates stated numerically in Results linked to figures; not digitized. "
                   "Residence time calculated as 1/koff.")))
    for r in observations:
        r["koff_s-1"] = r.pop("koff_s_1")
    return observations, added_potencies


def potency_rows(activities, supplemental, observations):
    accepted = {"EC50", "AC50", "SC1.5", "Kd", "Ki", "IC50"}
    rows = []
    excluded_measurements = defaultdict(set)
    for a in activities:
        if a["standard_type"] not in accepted or a["standard_value"] is None:
            continue
        if re.search(r"\b(?:MBP|myelin basic protein) as substrate\b",
                     a["assay_description"] or "", re.IGNORECASE):
            if a["pubmed_id"] and a["standard_units"] == "nM":
                key = (str(a["pubmed_id"]), a["molecule_chembl_id"], a["standard_type"],
                       a["standard_relation"])
                excluded_measurements[key].add(number(a["standard_value"]))
            continue
        rows.append(chem_potency(a))
    baseline = {r["record_id"]: r for r in observations if r["source_doi_or_patent"] == BORZILLERI}
    mollica = {r["record_id"]: r for r in observations if r["source_doi_or_patent"] == MOLLICA}
    activities_by_id = {a["activity_id"]: a for a in activities}
    for row in rows:
        a = activities_by_id[row["source_record"]]
        if a["doi"] == BORZILLERI:
            identity = baseline[a["record_id"]]
            row.update(compound_name=identity["compound_name"], smiles=identity["smiles"],
                       compound_id=identity["compound_key"], notes=identity["structure_status"],
                       potency_assay_conditions=BORZ_SPR if row["potency_type"] == "Kd" else BORZ_POTENCY)
        elif a["doi"] == MOLLICA:
            identity = mollica[a["record_id"]]
            row.update(compound_name=identity["compound_name"], smiles=identity["smiles"],
                       compound_id=identity["compound_key"], notes=identity["structure_status"],
                       potency_assay_conditions=identity["kinetics_conditions"])
    for a in csv_read("bindingdb_gck_activities.csv"):
        for kind in ("EC50", "IC50", "Kd", "Ki"):
            raw = a.get(f"{kind} (nM)", "").strip()
            if not raw:
                continue
            match = re.fullmatch(r"([<>=~≤≥]*)(.*)", raw)
            relation, value = match[1] or "=", match[2].strip()
            if a["Curation/DataSource"] == "ChEMBL" and a["PMID"]:
                ligand = a["ChEMBL ID of Ligand"] or a["BindingDB Ligand Name"].split("::", 1)[0]
                excluded_values = excluded_measurements.get((a["PMID"], ligand, kind, relation))
                if excluded_values:
                    measured = number(value)
                    if measured in excluded_values or (value.isdigit() and len(excluded_values) == 1
                                                       and round(next(iter(excluded_values))) == measured):
                        continue
            ph, temp = a["pH"].strip(), a["Temp (C)"].strip()
            conditions = "; ".join(x for x in [f"pH {ph}" if ph else "",
                                              f"{temp} °C" if temp else "",
                                              "glucose concentration not available in TSV"] if x)
            monomer = "BDBM" + a["BindingDB MonomerID"]
            rows.append(format_potency(a["BindingDB Ligand Name"] or monomer, a["Ligand SMILES"],
                                       kind, value, "nM", conditions,
                                       a["Article DOI"] or a["Patent Number"] or a["BindingDB Entry DOI"],
                                       "bindingdb", a["BindingDB Reactant_set_id"], monomer,
                                       relation, "BindingDB source SMILES; not recanonicalized"))
    for example, value in ((1, "3.6"), (5, "0.074"), (8, "0.16"), (17, "0.35"), (18, "0.43")):
        rows.append(format_potency(
            f"US8492405B2 Example {example}", "", "EC50", value, "µM",
            "Human GST-tagged liver GCK, 5 mM final glucose, 37 °C, Mant-ATP/Mant-ADP fluorescence",
            "US8492405B2", "patent", f"Experimental Example 1, Table 1, Example {example}", "",
            notes="Example structure not independently normalized; not used for a kinetics join."))
    return rows + supplemental


def candidates(activities, sources):
    rows = []
    fields = ["compound_name", "smiles", "source", "chembl_id", "available_evidence", "reason"]
    for a in activities:
        if a["doi"] != GKRP or a["standard_type"] != "Kd":
            continue
        entry = sources[a["record_id"]]
        rows.append(dict(zip(fields, (f'GKRP disrupter {entry["compound_key"]}', a["canonical_smiles"],
                                      GKRP, a["molecule_chembl_id"],
                                      f'human GCK SPR Kd {a["standard_relation"]} {a["standard_value"]} nM',
                                      "No human GCK association/dissociation rates verified; affinity alone is not kinetics."))))
    write_csv("gck_kinetics_candidates.csv", rows, fields)


def main():
    missing = [name for name in RAW_INPUTS if not (DATA_DIR / name).is_file()]
    if missing:
        commands = "\n".join(f"  python fetch_databases.py {RAW_INPUTS[name]}" for name in missing)
        raise SystemExit(
            f"Missing inputs in {DATA_DIR}: {', '.join(missing)}.\n"
            f"Run from gck-kinetics/scripts first:\n{commands}"
        )
    activities = json_read("chembl_gck_activities.json")
    sources = {r["record_id"]: r for r in json_read("chembl_source_compounds.json")}
    observations, supplemental = build_kinetics(activities, sources)
    ordered_dois = [BORZILLERI, ANTOINE, MOLLICA]
    observations.sort(key=lambda r: (
        not (r["baseline"] == "yes"), int(r["compound_key"].split(":")[1]),
        ordered_dois.index(r["source_doi_or_patent"]), str(r["record_id"])))
    joined = {}
    for row in observations:
        joined.setdefault(row["compound_key"], row)
    write_csv("gck_kinetics_joined.csv", joined.values(), JOIN)
    write_csv("gck_kinetics_observations.csv", observations, JOIN)
    potency = potency_rows(activities, supplemental, observations)
    write_csv("gck_potency_all.csv", potency, POTENCY)
    candidates(activities, sources)
    print(f"Joined: {len(joined)} samples; {sum(r['baseline'] == 'yes' for r in joined.values())} baseline; "
          f"{sum(r['baseline'] == 'no' for r in joined.values())} new")
    print(f"Measurements: {len(observations)} kinetic observations; {len(potency)} potency records")


if __name__ == "__main__":
    main()
