# Aurora A ChemBERTa A/B prespec

The primary panel (n=15) is underpowered for the 0.3 Δrho decision bar; a null is uninformative, not evidence of no effect. The Lake-coupled run is a sensitivity analysis only, never the headline.

## Endpoints

- **Primary (independent):** `CHEMBL3991630`, `CHEMBL3991601`, Klaeger et al., *The target landscape of clinical kinase drugs*, Science (2017), DOI `10.1126/science.aan4368`, PMID 29191878; not Lake et al. (2018). Assay type B (binding), in vitro subcellular Kinobeads epsilon competition with mass-spectrometry readout on human cell lysates (K-562, COLO 205, SK-N-BE(2), MV-4-11 listed by ChEMBL). Endogenous human Aurora A; no engineered variant or purified construct, phosphorylation or TPX2 occupancy specified. The assay reports apparent Kd in nM. Select one positive, uncensored Kd per panel compound by full InChIKey. Score `pKd = 9 - log10(Kd_nM)`; larger is stronger apparent binding. Aurora A therefore uses a different endpoint from the other targets.
- **Secondary (coupled) sensitivity:** Lake et al. PNAS (2018), DOI `10.1073/pnas.1811158115`, SI Fig. S9B/C, phosphorylated Aurora A without TPX2, 24 positive mean Kd values in nM; score `9 - log10(Kd_nM)`. The label is part of the same global fit used to derive Model B's cooperativity factor, so this comparison is circular and not evidence of independent predictive gain.

### Primary panel (full InChIKey)

| Compound | InChIKey |
|---|---|
| AMG-900 | IVUGFMLRJOCGAS-UHFFFAOYSA-N |
| AT-9283 | LOLPPWBBNUVNQZ-UHFFFAOYSA-N |
| BI-847325 | FLBNLJLONKAPLR-DQSJHHFOSA-N |
| CYC-116 | GPSZYOIFQZPWEJ-UHFFFAOYSA-N |
| ENMD-2076 | BLQYVHBZHAISJM-CMDGGOBGSA-N |
| GSK-1070916 | QTBWCSQGBMPECM-UHFFFAOYSA-N |
| KW-2449 | YYLKKYCXAOBSRM-JXMROGBWSA-N |
| MK-5108 | LCVIRAZGMYMNNT-VVONHTQRSA-N |
| MLN-8054 | HHFBDROWDBDFBR-UHFFFAOYSA-N |
| MLN-8237 | ZLHFILGSQDJULK-UHFFFAOYSA-N |
| PF-03814735 | RYYNGWLOYLRZLK-RBUKOAKNSA-N |
| SNS-314 | FAYAUAZLLLJJGH-UHFFFAOYSA-N |
| TAK-901 | WKDACQVEJIVHMZ-UHFFFAOYSA-N |
| VX-680 | GCIKSSRWRFVXBI-UHFFFAOYSA-N |
| danusertib | XKFTZKGMDDZMJI-HSZRJFAPSA-N |

## Features and folds

Canonical isomeric SMILES from the frozen source panel. No geometry columns survive: no Aurora A pocket geometry table is available, no gate is relaxed, and no predicted complexes are used. Model A is ligand-only, matching PPARγ. PXD062561 offers ligand-resolved HDX for 0/24 compounds, so Model B adds exactly one feature: `log10(published α)` (Lake SI Table S2). All 24 have a positive published mean, including SNS-314; retain its mean but flag the poorly constrained fit (only a lower bound α > 10 was established). Do not use Kd, IC50 or HDX as features. Standardise this one column using the training-fold mean and population SD, constant SD replaced by 1; Model A uses exactly zero in that same slot **after** preprocessing. Head shape, initialisation and trainable-parameter count are identical.

Leave-one-typed-Bemis-Murcko-scaffold-out, RDKit `GetScaffoldForMol` then canonical non-isomeric SMILES; sorted scaffold SMILES determine fold IDs shared across endpoints. No validation fold. Hold out all compounds with the named scaffold; train on the complement within each endpoint. The two endpoints use 14 and 23 folds respectively; folds marked `yes` are present for the primary.

| Fold | Held-out members (coupled panel) | Primary |
|---|---|---|
| S01 | ENMD-2076 | yes |
| S02 | CCT-129202 | no |
| S03 | TAK-901 | yes |
| S04 | ZM-447439 | no |
| S05 | VX-680 | yes |
| S06 | SNS-314 | yes |
| S07 | AT-9283 | yes |
| S08 | danusertib | yes |
| S09 | PHA-680632 | no |
| S10 | KW-2449 | yes |
| S11 | MK-8745 | no |
| S12 | JNJ-7706621 | no |
| S13 | hesperadin | no |
| S14 | BI-847325 | yes |
| S15 | CYC-116 | yes |
| S16 | CCT-137690 | no |
| S17 | PF-03814735 | yes |
| S18 | GSK-1070916 | yes |
| S19 | TAE-684 | no |
| S20 | MLN-8054; MLN-8237 | yes |
| S21 | MK-5108 | yes |
| S22 | AMG-900 | yes |
| S23 | reversine | no |

## Model and training

Frozen reference protocol SHA256 `062b06c15a7e7eeec35aad09f5dc39ba1275b60a202246c74b48d1885847eec1`. Load `seyonec/ChemBERTa-zinc-base-v1` offline from revision `761d6a18cf99db371e0b43baf3e2d21b3e865a20` and verify file hashes and no missing encoder weights. Tokenize with no augmentation, max length 256, truncation, dynamic padding and attention masks. Attention-mask-weighted mean of the final hidden states includes special tokens. Encoder dropout zero. Head: `Linear(769, 32)` → GELU → `Linear(32, 1)`, no dropout.

Headline LoRA: rank 8, alpha 16, dropout 0, bias none; query and value in every self-attention layer, all base parameters frozen. Supporting variant: frozen encoder, same head, no adapters. Per fold/label/variant/arm/seed reset Python, NumPy and Torch seeds to 20261123, 20261124, 20261125, 20261126 or 20261127, and use the same seed for independently shuffled batch order. Execution order: fold, label (primary then coupled), variant (LoRA then frozen), arm (A then B), then seed. Deterministic CPU float32. Optimiser: AdamW over trainable parameters, lr 1e-4, betas 0.9/0.999, eps 1e-8, weight decay 0.01. MSE on labels standardised by training-fold mean and population SD. Batch 8, 20 epochs, global grad-norm clip 1.0, constant LR, zero warmup, no accumulation, no early stopping, final epoch only. Invert predictions to pKd before metrics. A/B optimiser step counts match in each fold and seed.

## Metrics and decision

Headline aggregation: arithmetic mean of the five seed-specific OOF predictions per compound; pooled OOF predictions mix models trained on different subsets. Spearman uses average ranks, undefined for constants or fewer than two compounds. Report pooled and descriptive per-fold Spearman and pKd RMSE per arm; all five individual-seed pooled rhos and their mean, population SD, min and max. Supporting pairwise accuracy: skip tied true pairs, credit half for tied predictions. Top-5 enrichment (pooled only): fraction of the predicted top five in the actual top five divided by random expectation `5/n`; ties at the cutoff get fractional membership, making the result invariant to compound order. Supporting metrics do not determine a winner.

Per label/variant pooled comparison, reset the bootstrap RNG seed to 20261123, draw 1,500 resamples of n compound indices with replacement, use the same indices for A/B/truth and rerank within each sample. Use 2.5th/97.5th percentiles on defined `rho_B - rho_A` values; report the undefined fraction. Model B wins a cell only if pooled Δrho ≥ 0.3 **and** the lower paired 95% bound > 0. Otherwise report a null/sub-threshold result as uninformative at this n, not evidence of no effect. The coupled result cannot support an independent predictive claim regardless of whether it clears the numerical bar.

## Controls and interpretation

The numerical bar is necessary, not sufficient for a cooperativity-specific gain. Run these controls for both endpoints and both model variants, with the same folds and seed-mean OOF aggregate. For each trained control use the identical pKd scaling, head, optimiser, step count, RNG resets and batch order as its corresponding A/B fit.

1. **Fold mean:** predict the training-fold mean pKd for every held-out compound, a non-trained deterministic baseline so all five seeds would be identical.
2. **Row shuffle:** ten independent permutations of the single log10(α) column across the sorted compounds eligible for each endpoint, drawn once per endpoint and draw and then fixed across folds, variants and training seeds. Standardise using each training fold only.
3. **Gaussian sham:** ten independent standard normal draws of the same column shape as the shuffled block, fixed across folds, variants and training seeds and standardised using each training fold only.
4. **Cooperativity only:** use the true log10(α) column with 768 exact-zero embedding slots in place of the ChemBERTa embedding, preserving the head shape and initialisation. The bypassed encoder and adapters cannot affect prediction; use the corresponding variant's seeded head and train with the same optimizer schedule.

For sham draw 1–10, use NumPy `default_rng(SeedSequence([20261123, label_index, block_index, draw]))`, where label order is primary (0), coupled (1), and block order is shuffle (0), Gaussian (1). Controls follow A/B within each fold/endpoint/variant, in order: features-only, shuffles 1–10, Gaussian draws 1–10; fold mean requires no fit. Report pooled rho and pKd RMSE for the fold-mean and features-only controls and every sham draw. Compute each sham's Δrho versus the same real Model A; report counts out of 10 with defined Δrho at or above the real Δrho, and report undefined draws explicitly. Any sham reproducing a positive real Δrho means the gain is not attributable to the cooperativity signal. The numerical bar on a coupled comparison cannot overcome its circularity.

Report Model A's per-fold population prediction SD versus observed pKd SD and the ratio when defined. Call a non-singleton fold near-constant if prediction SD is at most 10% of its held-out label SD, provided label SD is nonzero; singleton and constant-label folds cannot diagnose this. Report its pooled rank correlation and RMSE against fold-mean-only predictions (not only against the truth); classify rank tracking at rho ≥ 0.8. Pooled agreement with training means and uninformative singleton folds limit interpretation of a large Δrho.
