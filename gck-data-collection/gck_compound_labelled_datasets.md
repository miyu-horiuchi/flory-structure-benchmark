# Compound-labelled GCK protein-state measurements: dataset survey and potency joins

Scope: human glucokinase (GCK / hexokinase IV, UniProt P35557). Searched RCSB PDB (Search + Data API, accession / entity-name / 70%-identity sequence clusters), ChEMBL (CHEMBL3820), BindingDB annotations served through RCSB, PRIDE, BMRB, SASBDB, Zenodo / figshare / OSF / MDverse, PDB-Dev, and Europe PMC full text (HDX, HDX-MS, pulsed proteolysis, DSF/thermal shift, SAXS, ITC, NMR, deposited MD), plus the data-availability statements of the hits. No training was done.

## Dataset table

| # | Dataset / source | Method (protein-state readout) | Per-compound? | Distinct GCK compounds | Compounds with potency labels | Join size | Access |
|---|---|---|---|---:|---:|---:|---|
| 1 | RCSB PDB human GCK crystal set | X-ray crystal structures (single-model, static) | yes | 26 activators + 1 putative inhibitor (0H5) across 35 human entries | 22 (ChEMBL EC50; 19 of them also BindingDB-via-PDB) | **22 compounds / 25 deposits** | RCSB; `pdb_gck_entries.csv`, `pdb_chembl_join.csv` |
| 2 | Liu et al. 2012, JBC (Pfizer), PMC3340189 | SAXS Rg glucose-titration = solution conformational-state fraction, per GKA | yes | 4 GKAs (compounds A–D) | 4 (EC50 691, 174, 92, 343 nM, quoted in the paper with α/β activator parameters) | 4 (paper-internal; no deposited curves) | PMC3340189; SAXS curves never deposited (no SASBDB entry) |
| 3 | Assemat et al. 2015, Anal Biochem | 19F NMR FAXS: R2 / spin-echo signal of a fluorinated spy ligand, per screened molecule | yes | 8-molecule library (+ TAFMT spy ligand, + reference GKAs) | 8 (signal correlates with EC50, EC50s are Servier in-house) | 8 in-paper; 0 publicly resolvable (structures undisclosed, no deposition) | doi:10.1016/j.ab.2015.02.006 |
| 4 | Bowler et al. 2013, Biochemistry (Miller lab) | ITC thermodynamics + transient kinetics of activator binding | yes | 3 (GKA22, "compound A", RO-28-1675) | 3 (published EC50/kinetics) | 3 in-paper; no deposited thermograms | PMC article; no repository deposit |
| 5 | Roman et al. 2025, Diabetes (dorzagliatin vs MK-0941), PMC12278794 | X-ray + enzyme kinetics + islet perifusion + 2 ns MD | yes | 2 | 2 (both clinical GKAs with published EC50) | 2 | PMC12278794; the dorzagliatin complex is not yet in the public PDB snapshot (full-text search returns no entry), MK-0941 pose taken from 1V4S/3IDH-based modelling |
| 6 | Sternisha et al. 2018 + Gordon et al. 2023 (Methods Enzymol) HDX-MS / pulsed proteolysis | HDX-MS, pulsed proteolysis (true dynamics) | **no** | 0 (unliganded WT and activating *variants*, glucose only) | n/a | 0 | PRIDE has no matching project (PXD061668, PXD037264 are unrelated GCK proteomics) |
| 7 | Larion / Miller methyl-TROSY and 13C NMR relaxation series | NMR millisecond dynamics (17 Ile methyls, 3 Trp) | **no** | 0 small-molecule activators (glucose-bound / apo) | n/a | 0 | BMRB has no activator-bound GCK deposition |
| 8 | Copenhagen MD deposit `10.17894/ucph.4e38c597-fae2-4654-a48a-8e714263c5a1` | All-atom MD trajectories | **no** | 0 (WT + 5 disease variants, no ligand) | n/a | 0 | ERDA / UCPH |
| 9 | PDB NMR, cryo-EM or multi-model GCK ensembles | — | — | 0 entries exist | n/a | 0 | RCSB |
| 10 | Zenodo / figshare / OSF / MDverse ligand-bound GCK MD | MD trajectories with a GKA bound | — | none found (Zenodo API returned HTTP 403; figshare/OSF/MDverse searches returned no GCK trajectory set) | n/a | 0 | — |

## Largest joinable cohort

**22 distinct GCK activators across 25 human crystal deposits, every one carrying an EC50 label** (ChEMBL CHEMBL3820 exact InChIKey/UniChem identity; 19 also annotated via BindingDB on the PDB entry). Only MRK (1V4S, 3F9M, 3ID8) and S41 (3S41, 4DHY) have more than one deposit. Unjoined drug-like ligands: 0H4, 0H5, 1FY, 8WJ, TK1 (1FY matches only the racemate CHEMBL2346884, not the crystallized (2S) enantiomer).

Caveat that matters for the Model A vs Model B design: this cohort is *static structure*, not dynamics. The largest cohort with an actual per-compound protein-state/dynamics readout is **4 compounds** (Liu 2012 SAXS, each with EC50), or 8 if the undisclosed Servier 19F NMR library were obtainable. Compound-labelled GCK dynamics data at the scale needed to test the dynamics hypothesis does not exist in public repositories; the blocker reported earlier is confirmed, and the only path at n > 20 is crystal-structure features, or generating dynamics in-house from the 25 liganded deposits.

## Supporting files

`pdb_gck_entries.csv` (37 deposits: 35 human, 1 Xenopus, 1 ancestral; method, resolution, mutations, all nonpolymer entities, binding-affinity annotations), `pdb_gck_activators.csv`, `chembl_gck_activities.csv` (3,222 activities, 1,135 compounds with numeric EC50/AC50/IC50/Ki), `pdb_chembl_join.csv` (per-ligand identity maps, PDB IDs, 53 BindingDB annotation rows), `part1_summary.md`, and the collection scripts `collect.py` / `summarize.py`.
