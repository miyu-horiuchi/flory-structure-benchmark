# Personal AlphaFold 3 Server comparison

**Personal, non-commercial exploratory use only. Shown in the same comparison table as the other seven models, but not company research or a fully input-matched eight-model leaderboard. The original seven-model data and AF3's separate usage restrictions are unchanged.**

**This information is subject to AlphaFold Server Output Terms of Use found at https://alphafoldserver.com/output-terms**

Read the [preserved full output terms](output-terms.txt). These restrictions apply to this comparison and its derived metrics, including the [results/provenance JSON](results.json). No commercial/company use, automated docking or interaction-screening systems, or training similar biomolecular structure-prediction models. No Google endorsement is implied. Personal GitHub hosting does not change these restrictions or make this public page private.

## Results

Twenty native predictions were generated: five per target with seed 42, coordinate-free protein sequences, the frozen shared unpaired MSAs and templates disabled. Partner chains and native PTMs were retained. The frozen version-3 evaluator and experimental references were reused.

| Target | Selected ligand RMSD after receptor alignment | Strict pass / generated | Complete / scorable |
|---|---:|---:|---:|
| PPARγ | 1.40 Å | 5/5 | 5/5 |
| GCK | Not scorable: glucose O1 missing | 0/5 | 0/5 |
| AURKA | 0.93 Å | 5/5 | 5/5 |
| ERK2 | 0.48 Å | 3/5 | 5/5 |

Overall, 13/20 generated predictions passed, or 13/15 complete/scorable predictions. The five GCK outputs are evaluation errors, not five measured poor poses. The highest native `ranking_score` selected the displayed prediction before geometry screening; ties use the lower sample index. Sample 0 was selected for all four targets.

## Important limitations

- **GCK:** every exported glucose ligand contains 11 heavy atoms instead of 12, with O1 missing. No atom was reconstructed and no partial-ligand RMSD was substituted.
- **AURKA:** the server's CCD VX6 definition is neutral; the earlier benchmark used a +1 VX-680 form. Canonical isomeric SMILES differ. Native phospho-Thr288 and TPX2 were retained, but this is not a chemically identical input comparison.
- **Partner placement:** selected PPARγ SRC-1 partner Cα RMSD is 24.23 Å after receptor alignment, despite a good ligand pose; selected AURKA TPX2 is 1.80 Å.
- **ERK2:** two non-selected samples have severe ligand–protein clashes.
- These retrospective public cases may overlap training. One seed and five samples do not establish generalisation, functional-state accuracy, dynamics, kinetics, thermodynamics or fine-tuning suitability.

The JSON retains the four cells, all twenty sample records, native confidence, evaluation errors, atom inventories, chemical-input comparisons, CIF hashes, original job names, model version strings, sequence/MSA/reference hashes and input/output manifests. Coordinate arrays are omitted from this display export; raw predicted structures are not redistributed on this site. Original analysis-file hashes and the frozen private evidence-archive hash identify the source material; the original seven-model dataset is unchanged.

## Citation

Abramson et al. (2024), *Accurate structure prediction of biomolecular interactions with AlphaFold 3*, Nature 630, 493–500. https://doi.org/10.1038/s41586-024-07487-w
