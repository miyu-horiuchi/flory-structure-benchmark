# Aurora A ChemBERTa A/B

**The sealed decision rule did not resolve a Model B benefit at n=15.** The primary LoRA Δρ was 0.0733; the rule requires Δρ ≥0.3 and a positive lower 95% confidence bound. This result is not evidence against cooperativity or dynamics.

13 of 14 occupied primary folds contain one compound. The fold-mean-only predictor scores ρ -0.9991: on singleton folds, holding out a higher pKd makes the training mean lower by construction. This arithmetic property is expected and does not diagnose model failure or anti-prediction.

Model A’s rank correlation with the training-fold mean is 0.5170 for LoRA and 0.5796 for frozen. Sham draws at or above the real Δρ: 0/10 shuffled and 1/10 Gaussian for LoRA; 8/10 shuffled and 9/10 Gaussian for frozen. The positive LoRA gain is not attributable to cooperativity.

In the protocol owner’s power analysis, the unchanged decision rule fired 38% of the time at n=26 with a planted true Δρ of exactly 0.3. Regrouping folds changed this by less than one percentage point. The CI-excludes-zero condition alone had 31% power at n=26 and needed about 104 compounds for 80% power. At n=15 the sealed result is unresolved; grouped k-fold is a comparability analysis, not a rescue.

**Sealed primary:** The n=15 panel is underpowered for the 0.3 Δrho bar. A null is uninformative, not evidence of no effect.

Model A is a ligand-only baseline, matching PPARγ: no geometry columns survive. PXD062561 has 0/24 ligand-resolved HDX coverage. Model B therefore adds only training-fold-standardised log10(cooperativity α); no substitutes or predicted complexes were used. SNS-314 remains included despite its poorly constrained α (only a lower bound >10 could be established).

The primary endpoint is pKd from Klaeger et al. (2017), human cell-lysate Kinobeads binding assay CHEMBL3991630, not Lake data. Aurora A is scored on a different endpoint from the other three targets. The Lake no-TPX2 pKd sensitivity is **circular**: it shares the global fit that defined α, and is not headline or independent evidence.

Prespec SHA256: `39fb48b381d695970166c2120299d640af577a3bf03d7306d4727fb23fafc03e`.

The sealed prespec names fold-first execution. The requester asked to finish the primary endpoint before starting the Lake-coupled sensitivity. Only dispatch order across independent endpoint fits changed; the sealed prespec and its SHA256 above are unchanged.

The first primary check stopped with an `IndexError`: a local loop variable shadowed the feature lookup in `check.py`. The checker-only correction did not change any fit or computed metric.

Metrics use seed-mean out-of-fold predictions (pooled across models trained on different subsets). The RMSE scale is pKd. A fold with one compound has undefined rho; a two-compound rho is an ordering sign, not a stable estimate.

| Endpoint | Variant | n / folds | ρ A | ρ B | RMSE A | RMSE B | Δρ | 95% paired CI | Undefined bootstrap | Shuffle / Gaussian at or above Δρ | Decision |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|
| Independent | lora | 15 / 14 | -0.5147 | -0.4415 | 0.9253 | 0.9230 | 0.0733 | [0.0000, 0.2582] | 0.00% | 0/10 / 1/10 | Sham reproduced the gain; not attributable to cooperativity |
| Independent | frozen | 15 / 14 | -0.5773 | -0.5773 | 0.9238 | 0.9240 | 0.0000 | [0.0000, 0.0000] | 0.00% | 8/10 / 9/10 | Null/sub-threshold; uninformative |
| Lake coupled (sensitivity) | lora | 24 / 23 | -0.1374 | -0.0791 | 0.9378 | 0.9336 | 0.0583 | [-0.0446, 0.1742] | 0.00% | 0/10 / 1/10 | Sham reproduced the gain; not attributable to cooperativity |
| Lake coupled (sensitivity) | frozen | 24 / 23 | -0.1817 | -0.1643 | 0.9304 | 0.9293 | 0.0174 | [-0.0706, 0.1077] | 0.00% | 3/10 / 3/10 | Sham reproduced the gain; not attributable to cooperativity |

## Controls

Sham Δρ is relative to the same real Model A. Near-constant means prediction SD ≤10% of label SD on an evaluable fold. Singleton or constant-label folds cannot diagnose near-constancy. Rank tracking means ρ ≥0.8 against the fold-mean predictor.

| Endpoint | Variant | Fold-mean ρ / RMSE | Cooperativity-only ρ / RMSE | A vs fold mean ρ / RMSE | A tracks fold means by rank | Near-constant / evaluable folds |
|---|---|---|---|---|---|---|
| Independent | lora | -0.9991 / 0.9311 | -0.9759 / 0.9292 | 0.5170 / 0.1146 | No | 1/1 |
| Independent | frozen | -0.9991 / 0.9311 | -0.9955 / 0.9358 | 0.5796 / 0.0919 | No | 1/1 |
| Lake coupled (sensitivity) | lora | -0.9928 / 0.9509 | -0.9296 / 0.9475 | 0.1396 / 0.1294 | No | 0/1 |
| Lake coupled (sensitivity) | frozen | -0.9928 / 0.9509 | -0.9809 / 0.9541 | 0.1931 / 0.0917 | No | 0/1 |

| Endpoint | Variant | Sham block | Draw | ρ | RMSE | Δρ vs A |
|---|---|---|---:|---:|---:|---:|
| Independent | lora | shuffle | 1 | -0.5183 | 0.9239 | -0.0036 |
| Independent | lora | shuffle | 2 | -0.5147 | 0.9257 | 0.0000 |
| Independent | lora | shuffle | 3 | -0.5147 | 0.9279 | 0.0000 |
| Independent | lora | shuffle | 4 | -0.5147 | 0.9264 | 0.0000 |
| Independent | lora | shuffle | 5 | -0.5255 | 0.9248 | -0.0107 |
| Independent | lora | shuffle | 6 | -0.4719 | 0.9260 | 0.0429 |
| Independent | lora | shuffle | 7 | -0.4915 | 0.9261 | 0.0232 |
| Independent | lora | shuffle | 8 | -0.4558 | 0.9237 | 0.0590 |
| Independent | lora | shuffle | 9 | -0.4683 | 0.9246 | 0.0465 |
| Independent | lora | shuffle | 10 | -0.5147 | 0.9245 | 0.0000 |
| Independent | lora | gaussian | 1 | -0.5147 | 0.9264 | 0.0000 |
| Independent | lora | gaussian | 2 | -0.5183 | 0.9224 | -0.0036 |
| Independent | lora | gaussian | 3 | -0.5147 | 0.9261 | 0.0000 |
| Independent | lora | gaussian | 4 | -0.5183 | 0.9279 | -0.0036 |
| Independent | lora | gaussian | 5 | -0.5147 | 0.9248 | 0.0000 |
| Independent | lora | gaussian | 6 | -0.5147 | 0.9235 | 0.0000 |
| Independent | lora | gaussian | 7 | -0.4272 | 0.9213 | 0.0876 |
| Independent | lora | gaussian | 8 | -0.4754 | 0.9238 | 0.0393 |
| Independent | lora | gaussian | 9 | -0.5183 | 0.9267 | -0.0036 |
| Independent | lora | gaussian | 10 | -0.5147 | 0.9224 | 0.0000 |
| Independent | frozen | shuffle | 1 | -0.5702 | 0.9194 | 0.0071 |
| Independent | frozen | shuffle | 2 | -0.5934 | 0.9234 | -0.0161 |
| Independent | frozen | shuffle | 3 | -0.5773 | 0.9255 | 0.0000 |
| Independent | frozen | shuffle | 4 | -0.5916 | 0.9227 | -0.0143 |
| Independent | frozen | shuffle | 5 | -0.5773 | 0.9249 | 0.0000 |
| Independent | frozen | shuffle | 6 | -0.5773 | 0.9252 | 0.0000 |
| Independent | frozen | shuffle | 7 | -0.5773 | 0.9246 | 0.0000 |
| Independent | frozen | shuffle | 8 | -0.5719 | 0.9244 | 0.0054 |
| Independent | frozen | shuffle | 9 | -0.5773 | 0.9244 | 0.0000 |
| Independent | frozen | shuffle | 10 | -0.5451 | 0.9209 | 0.0322 |
| Independent | frozen | gaussian | 1 | -0.5773 | 0.9251 | 0.0000 |
| Independent | frozen | gaussian | 2 | -0.5451 | 0.9169 | 0.0322 |
| Independent | frozen | gaussian | 3 | -0.5934 | 0.9234 | -0.0161 |
| Independent | frozen | gaussian | 4 | -0.5773 | 0.9245 | 0.0000 |
| Independent | frozen | gaussian | 5 | -0.5451 | 0.9214 | 0.0322 |
| Independent | frozen | gaussian | 6 | -0.5308 | 0.9183 | 0.0465 |
| Independent | frozen | gaussian | 7 | -0.5719 | 0.9233 | 0.0054 |
| Independent | frozen | gaussian | 8 | -0.5773 | 0.9241 | 0.0000 |
| Independent | frozen | gaussian | 9 | -0.5773 | 0.9242 | 0.0000 |
| Independent | frozen | gaussian | 10 | -0.5308 | 0.9165 | 0.0465 |
| Lake coupled (sensitivity) | lora | shuffle | 1 | -0.1574 | 0.9398 | -0.0200 |
| Lake coupled (sensitivity) | lora | shuffle | 2 | -0.1617 | 0.9393 | -0.0243 |
| Lake coupled (sensitivity) | lora | shuffle | 3 | -0.0939 | 0.9368 | 0.0435 |
| Lake coupled (sensitivity) | lora | shuffle | 4 | -0.1209 | 0.9364 | 0.0165 |
| Lake coupled (sensitivity) | lora | shuffle | 5 | -0.1522 | 0.9396 | -0.0148 |
| Lake coupled (sensitivity) | lora | shuffle | 6 | -0.1922 | 0.9401 | -0.0548 |
| Lake coupled (sensitivity) | lora | shuffle | 7 | -0.1600 | 0.9387 | -0.0226 |
| Lake coupled (sensitivity) | lora | shuffle | 8 | -0.0878 | 0.9356 | 0.0496 |
| Lake coupled (sensitivity) | lora | shuffle | 9 | -0.1122 | 0.9380 | 0.0252 |
| Lake coupled (sensitivity) | lora | shuffle | 10 | -0.1330 | 0.9378 | 0.0043 |
| Lake coupled (sensitivity) | lora | gaussian | 1 | -0.1583 | 0.9396 | -0.0209 |
| Lake coupled (sensitivity) | lora | gaussian | 2 | -0.2052 | 0.9387 | -0.0678 |
| Lake coupled (sensitivity) | lora | gaussian | 3 | -0.1270 | 0.9381 | 0.0104 |
| Lake coupled (sensitivity) | lora | gaussian | 4 | -0.1652 | 0.9391 | -0.0278 |
| Lake coupled (sensitivity) | lora | gaussian | 5 | -0.1600 | 0.9394 | -0.0226 |
| Lake coupled (sensitivity) | lora | gaussian | 6 | -0.1600 | 0.9385 | -0.0226 |
| Lake coupled (sensitivity) | lora | gaussian | 7 | -0.0835 | 0.9368 | 0.0539 |
| Lake coupled (sensitivity) | lora | gaussian | 8 | -0.1400 | 0.9388 | -0.0026 |
| Lake coupled (sensitivity) | lora | gaussian | 9 | -0.1487 | 0.9403 | -0.0113 |
| Lake coupled (sensitivity) | lora | gaussian | 10 | -0.0791 | 0.9364 | 0.0583 |
| Lake coupled (sensitivity) | frozen | shuffle | 1 | -0.1991 | 0.9318 | -0.0174 |
| Lake coupled (sensitivity) | frozen | shuffle | 2 | -0.2157 | 0.9315 | -0.0339 |
| Lake coupled (sensitivity) | frozen | shuffle | 3 | -0.1948 | 0.9311 | -0.0130 |
| Lake coupled (sensitivity) | frozen | shuffle | 4 | -0.1452 | 0.9258 | 0.0365 |
| Lake coupled (sensitivity) | frozen | shuffle | 5 | -0.2104 | 0.9326 | -0.0287 |
| Lake coupled (sensitivity) | frozen | shuffle | 6 | -0.1800 | 0.9323 | 0.0017 |
| Lake coupled (sensitivity) | frozen | shuffle | 7 | -0.1426 | 0.9301 | 0.0391 |
| Lake coupled (sensitivity) | frozen | shuffle | 8 | -0.0783 | 0.9252 | 0.1035 |
| Lake coupled (sensitivity) | frozen | shuffle | 9 | -0.1748 | 0.9283 | 0.0070 |
| Lake coupled (sensitivity) | frozen | shuffle | 10 | -0.2052 | 0.9321 | -0.0235 |
| Lake coupled (sensitivity) | frozen | gaussian | 1 | -0.1670 | 0.9316 | 0.0148 |
| Lake coupled (sensitivity) | frozen | gaussian | 2 | -0.1965 | 0.9327 | -0.0148 |
| Lake coupled (sensitivity) | frozen | gaussian | 3 | -0.1443 | 0.9288 | 0.0374 |
| Lake coupled (sensitivity) | frozen | gaussian | 4 | -0.1600 | 0.9298 | 0.0217 |
| Lake coupled (sensitivity) | frozen | gaussian | 5 | -0.1826 | 0.9301 | -0.0009 |
| Lake coupled (sensitivity) | frozen | gaussian | 6 | -0.1557 | 0.9292 | 0.0261 |
| Lake coupled (sensitivity) | frozen | gaussian | 7 | -0.1809 | 0.9311 | 0.0009 |
| Lake coupled (sensitivity) | frozen | gaussian | 8 | -0.1783 | 0.9310 | 0.0035 |
| Lake coupled (sensitivity) | frozen | gaussian | 9 | -0.1852 | 0.9326 | -0.0035 |
| Lake coupled (sensitivity) | frozen | gaussian | 10 | -0.1870 | 0.9305 | -0.0052 |

## Supporting pooled and individual-seed metrics

| Endpoint | Variant | Arm | Pairwise accuracy | Top-5 enrichment | Seed pooled ρ (20261123–27) | Mean | SD | Min | Max |
|---|---|---|---:|---:|---|---:|---:|---:|---:|
| Independent | lora | A | 0.3654 | 0.0000 | -0.3914, -0.4879, -0.3735, -0.4450, -0.5541 | -0.4504 | 0.0657 | -0.5541 | -0.3735 |
| Independent | lora | B | 0.3846 | 0.0000 | -0.3003, -0.4844, -0.3610, -0.4218, -0.5612 | -0.4257 | 0.0914 | -0.5612 | -0.3003 |
| Independent | frozen | A | 0.2981 | 0.0000 | -0.5702, -0.5416, -0.5630, -0.4772, -0.4415 | -0.5187 | 0.0506 | -0.5702 | -0.4415 |
| Independent | frozen | B | 0.2981 | 0.0000 | -0.5880, -0.6023, -0.5630, -0.4665, -0.4415 | -0.5323 | 0.0656 | -0.6023 | -0.4415 |
| Lake coupled (sensitivity) | lora | A | 0.4457 | 0.0000 | -0.0278, -0.0913, -0.2852, -0.0852, -0.1948 | -0.1369 | 0.0916 | -0.2852 | -0.0278 |
| Lake coupled (sensitivity) | lora | B | 0.4710 | 0.9600 | -0.0200, -0.0226, -0.2096, -0.0287, -0.1496 | -0.0861 | 0.0787 | -0.2096 | -0.0200 |
| Lake coupled (sensitivity) | frozen | A | 0.4384 | 0.0000 | -0.2122, -0.2209, -0.2296, -0.0348, -0.1357 | -0.1666 | 0.0739 | -0.2296 | -0.0348 |
| Lake coupled (sensitivity) | frozen | B | 0.4420 | 0.0000 | -0.2339, -0.2365, -0.1939, -0.0504, -0.1530 | -0.1736 | 0.0687 | -0.2365 | -0.0504 |

## Scaffold membership (23 total typed scaffolds)

| Fold | Primary members | Coupled members |
|---|---|---|
| S01 | ENMD-2076 | ENMD-2076 |
| S02 | — | CCT-129202 |
| S03 | TAK-901 | TAK-901 |
| S04 | — | ZM-447439 |
| S05 | VX-680 | VX-680 |
| S06 | SNS-314 | SNS-314 |
| S07 | AT-9283 | AT-9283 |
| S08 | danusertib | danusertib |
| S09 | — | PHA-680632 |
| S10 | KW-2449 | KW-2449 |
| S11 | — | MK-8745 |
| S12 | — | JNJ-7706621 |
| S13 | — | hesperadin |
| S14 | BI-847325 | BI-847325 |
| S15 | CYC-116 | CYC-116 |
| S16 | — | CCT-137690 |
| S17 | PF-03814735 | PF-03814735 |
| S18 | GSK-1070916 | GSK-1070916 |
| S19 | — | TAE-684 |
| S20 | MLN-8054; MLN-8237 | MLN-8054; MLN-8237 |
| S21 | MK-5108 | MK-5108 |
| S22 | AMG-900 | AMG-900 |
| S23 | — | reversine |

## Per-fold: Independent, lora

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S01 | 1 | undefined | undefined | 0.6228 | 0.6098 | 0.0000 | 0.0000 | undefined |
| S03 | 1 | undefined | undefined | 0.2867 | 0.2807 | 0.0000 | 0.0000 | undefined |
| S05 | 1 | undefined | undefined | 0.5901 | 0.5890 | 0.0000 | 0.0000 | undefined |
| S06 | 1 | undefined | undefined | 0.9024 | 0.9186 | 0.0000 | 0.0000 | undefined |
| S07 | 1 | undefined | undefined | 0.2956 | 0.3089 | 0.0000 | 0.0000 | undefined |
| S08 | 1 | undefined | undefined | 0.2074 | 0.2063 | 0.0000 | 0.0000 | undefined |
| S10 | 1 | undefined | undefined | 1.3245 | 1.3135 | 0.0000 | 0.0000 | undefined |
| S14 | 1 | undefined | undefined | 0.5559 | 0.5727 | 0.0000 | 0.0000 | undefined |
| S15 | 1 | undefined | undefined | 0.8327 | 0.8203 | 0.0000 | 0.0000 | undefined |
| S17 | 1 | undefined | undefined | 0.0736 | 0.0613 | 0.0000 | 0.0000 | undefined |
| S18 | 1 | undefined | undefined | 0.7049 | 0.7171 | 0.0000 | 0.0000 | undefined |
| S20 | 2 | 1.0000 | 1.0000 | 0.8806 | 0.8777 | 0.0058 | 0.1990 | 0.0293 |
| S21 | 1 | undefined | undefined | 2.4374 | 2.4295 | 0.0000 | 0.0000 | undefined |
| S22 | 1 | undefined | undefined | 0.5733 | 0.5589 | 0.0000 | 0.0000 | undefined |

## Per-fold: Independent, frozen

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S01 | 1 | undefined | undefined | 0.6665 | 0.6647 | 0.0000 | 0.0000 | undefined |
| S03 | 1 | undefined | undefined | 0.2511 | 0.2493 | 0.0000 | 0.0000 | undefined |
| S05 | 1 | undefined | undefined | 0.5597 | 0.5603 | 0.0000 | 0.0000 | undefined |
| S06 | 1 | undefined | undefined | 0.9380 | 0.9406 | 0.0000 | 0.0000 | undefined |
| S07 | 1 | undefined | undefined | 0.2693 | 0.2675 | 0.0000 | 0.0000 | undefined |
| S08 | 1 | undefined | undefined | 0.2075 | 0.2086 | 0.0000 | 0.0000 | undefined |
| S10 | 1 | undefined | undefined | 1.3764 | 1.3768 | 0.0000 | 0.0000 | undefined |
| S14 | 1 | undefined | undefined | 0.5467 | 0.5498 | 0.0000 | 0.0000 | undefined |
| S15 | 1 | undefined | undefined | 0.8369 | 0.8369 | 0.0000 | 0.0000 | undefined |
| S17 | 1 | undefined | undefined | 0.0614 | 0.0591 | 0.0000 | 0.0000 | undefined |
| S18 | 1 | undefined | undefined | 0.6911 | 0.6928 | 0.0000 | 0.0000 | undefined |
| S20 | 2 | 1.0000 | 1.0000 | 0.8926 | 0.8928 | 0.0187 | 0.1990 | 0.0938 |
| S21 | 1 | undefined | undefined | 2.4078 | 2.4075 | 0.0000 | 0.0000 | undefined |
| S22 | 1 | undefined | undefined | 0.4670 | 0.4645 | 0.0000 | 0.0000 | undefined |

## Per-fold: Lake coupled (sensitivity), lora

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S01 | 1 | undefined | undefined | 0.1328 | 0.1491 | 0.0000 | 0.0000 | undefined |
| S02 | 1 | undefined | undefined | 0.4239 | 0.4159 | 0.0000 | 0.0000 | undefined |
| S03 | 1 | undefined | undefined | 0.0532 | 0.0597 | 0.0000 | 0.0000 | undefined |
| S04 | 1 | undefined | undefined | 1.6968 | 1.7098 | 0.0000 | 0.0000 | undefined |
| S05 | 1 | undefined | undefined | 0.0936 | 0.0875 | 0.0000 | 0.0000 | undefined |
| S06 | 1 | undefined | undefined | 0.1277 | 0.1516 | 0.0000 | 0.0000 | undefined |
| S07 | 1 | undefined | undefined | 0.5867 | 0.6029 | 0.0000 | 0.0000 | undefined |
| S08 | 1 | undefined | undefined | 0.1592 | 0.1535 | 0.0000 | 0.0000 | undefined |
| S09 | 1 | undefined | undefined | 0.0334 | 0.0214 | 0.0000 | 0.0000 | undefined |
| S10 | 1 | undefined | undefined | 0.7936 | 0.7772 | 0.0000 | 0.0000 | undefined |
| S11 | 1 | undefined | undefined | 1.2229 | 1.2001 | 0.0000 | 0.0000 | undefined |
| S12 | 1 | undefined | undefined | 0.7004 | 0.6842 | 0.0000 | 0.0000 | undefined |
| S13 | 1 | undefined | undefined | 0.4379 | 0.4275 | 0.0000 | 0.0000 | undefined |
| S14 | 1 | undefined | undefined | 0.1641 | 0.1877 | 0.0000 | 0.0000 | undefined |
| S15 | 1 | undefined | undefined | 0.5466 | 0.5286 | 0.0000 | 0.0000 | undefined |
| S16 | 1 | undefined | undefined | 0.1318 | 0.1402 | 0.0000 | 0.0000 | undefined |
| S17 | 1 | undefined | undefined | 0.3496 | 0.3643 | 0.0000 | 0.0000 | undefined |
| S18 | 1 | undefined | undefined | 1.2223 | 1.2438 | 0.0000 | 0.0000 | undefined |
| S19 | 1 | undefined | undefined | 0.4790 | 0.4725 | 0.0000 | 0.0000 | undefined |
| S20 | 2 | -1.0000 | -1.0000 | 0.9584 | 0.9491 | 0.0688 | 0.2852 | 0.2412 |
| S21 | 1 | undefined | undefined | 2.6522 | 2.6373 | 0.0000 | 0.0000 | undefined |
| S22 | 1 | undefined | undefined | 1.4890 | 1.4704 | 0.0000 | 0.0000 | undefined |
| S23 | 1 | undefined | undefined | 1.2445 | 1.2359 | 0.0000 | 0.0000 | undefined |

## Per-fold: Lake coupled (sensitivity), frozen

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S01 | 1 | undefined | undefined | 0.0261 | 0.0317 | 0.0000 | 0.0000 | undefined |
| S02 | 1 | undefined | undefined | 0.4027 | 0.3971 | 0.0000 | 0.0000 | undefined |
| S03 | 1 | undefined | undefined | 0.0072 | 0.0112 | 0.0000 | 0.0000 | undefined |
| S04 | 1 | undefined | undefined | 1.7052 | 1.7113 | 0.0000 | 0.0000 | undefined |
| S05 | 1 | undefined | undefined | 0.0069 | 0.0066 | 0.0000 | 0.0000 | undefined |
| S06 | 1 | undefined | undefined | 0.0869 | 0.0939 | 0.0000 | 0.0000 | undefined |
| S07 | 1 | undefined | undefined | 0.5549 | 0.5562 | 0.0000 | 0.0000 | undefined |
| S08 | 1 | undefined | undefined | 0.1624 | 0.1635 | 0.0000 | 0.0000 | undefined |
| S09 | 1 | undefined | undefined | 0.0225 | 0.0155 | 0.0000 | 0.0000 | undefined |
| S10 | 1 | undefined | undefined | 0.8069 | 0.8024 | 0.0000 | 0.0000 | undefined |
| S11 | 1 | undefined | undefined | 1.2461 | 1.2403 | 0.0000 | 0.0000 | undefined |
| S12 | 1 | undefined | undefined | 0.6877 | 0.6840 | 0.0000 | 0.0000 | undefined |
| S13 | 1 | undefined | undefined | 0.3714 | 0.3694 | 0.0000 | 0.0000 | undefined |
| S14 | 1 | undefined | undefined | 0.1723 | 0.1787 | 0.0000 | 0.0000 | undefined |
| S15 | 1 | undefined | undefined | 0.5370 | 0.5339 | 0.0000 | 0.0000 | undefined |
| S16 | 1 | undefined | undefined | 0.1084 | 0.1145 | 0.0000 | 0.0000 | undefined |
| S17 | 1 | undefined | undefined | 0.3449 | 0.3509 | 0.0000 | 0.0000 | undefined |
| S18 | 1 | undefined | undefined | 1.2280 | 1.2332 | 0.0000 | 0.0000 | undefined |
| S19 | 1 | undefined | undefined | 0.4755 | 0.4705 | 0.0000 | 0.0000 | undefined |
| S20 | 2 | -1.0000 | -1.0000 | 0.9692 | 0.9668 | 0.0629 | 0.2852 | 0.2206 |
| S21 | 1 | undefined | undefined | 2.6621 | 2.6573 | 0.0000 | 0.0000 | undefined |
| S22 | 1 | undefined | undefined | 1.3038 | 1.2986 | 0.0000 | 0.0000 | undefined |
| S23 | 1 | undefined | undefined | 1.3045 | 1.3020 | 0.0000 | 0.0000 | undefined |

## Post-hoc exploratory grouped-scaffold sensitivity: independent primary

**Post-hoc exploratory sensitivity, added 2026-11-23 UTC.** The sealed leave-one-scaffold-out primary remains the headline regardless of this result.

3 folds with 5, 5, 5 held-out compounds. Decreasing scaffold size, then canonical scaffold SMILES; assign intact to least-loaded fold, breaking ties by fold ID.

**Sealed primary:** The n=15 panel is underpowered for the 0.3 Δrho bar. A null is uninformative, not evidence of no effect.

Model A is a ligand-only baseline, matching PPARγ: no geometry columns survive. PXD062561 has 0/24 ligand-resolved HDX coverage. Model B therefore adds only training-fold-standardised log10(cooperativity α); no substitutes or predicted complexes were used. SNS-314 remains included despite its poorly constrained α (only a lower bound >10 could be established).

The primary endpoint is pKd from Klaeger et al. (2017), human cell-lysate Kinobeads binding assay CHEMBL3991630, not Lake data. Aurora A is scored on a different endpoint from the other three targets. The Lake no-TPX2 pKd sensitivity is **circular**: it shares the global fit that defined α, and is not headline or independent evidence.

Prespec SHA256: `39fb48b381d695970166c2120299d640af577a3bf03d7306d4727fb23fafc03e`.

The sealed prespec names fold-first execution. The requester asked to finish the primary endpoint before starting the Lake-coupled sensitivity. Only dispatch order across independent endpoint fits changed; the sealed prespec and its SHA256 above are unchanged.

The first primary check stopped with an `IndexError`: a local loop variable shadowed the feature lookup in `check.py`. The checker-only correction did not change any fit or computed metric.

Metrics use seed-mean out-of-fold predictions (pooled across models trained on different subsets). The RMSE scale is pKd. A fold with one compound has undefined rho; a two-compound rho is an ordering sign, not a stable estimate.

| Endpoint | Variant | n / folds | ρ A | ρ B | RMSE A | RMSE B | Δρ | 95% paired CI | Undefined bootstrap | Shuffle / Gaussian at or above Δρ | Decision |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|
| Independent (exploratory grouped) | lora | 15 / 3 | -0.0322 | 0.0107 | 0.8924 | 0.8921 | 0.0429 | [0.0000, 0.1989] | 0.00% | 2/10 / 2/10 | Sham reproduced the gain; not attributable to cooperativity |
| Independent (exploratory grouped) | frozen | 15 / 3 | -0.0608 | -0.0822 | 0.8968 | 0.8975 | -0.0214 | [-0.1273, 0.0000] | 0.00% | 9/10 / 9/10 | Null/sub-threshold; uninformative |

### Controls

Sham Δρ is relative to the same real Model A. Near-constant means prediction SD ≤10% of label SD on an evaluable fold. Singleton or constant-label folds cannot diagnose near-constancy. Rank tracking means ρ ≥0.8 against the fold-mean predictor.

| Endpoint | Variant | Fold-mean ρ / RMSE | Cooperativity-only ρ / RMSE | A vs fold mean ρ / RMSE | A tracks fold means by rank | Near-constant / evaluable folds |
|---|---|---|---|---|---|---|
| Independent (exploratory grouped) | lora | -0.2459 / 0.9167 | -0.1251 / 0.9157 | 0.8126 / 0.1032 | Yes | 1/3 |
| Independent (exploratory grouped) | frozen | -0.2459 / 0.9167 | -0.3360 / 0.9225 | 0.8126 / 0.0857 | Yes | 1/3 |

| Endpoint | Variant | Sham block | Draw | ρ | RMSE | Δρ vs A |
|---|---|---|---:|---:|---:|---:|
| Independent (exploratory grouped) | lora | shuffle | 1 | -0.0322 | 0.8910 | 0.0000 |
| Independent (exploratory grouped) | lora | shuffle | 2 | -0.0322 | 0.8928 | 0.0000 |
| Independent (exploratory grouped) | lora | shuffle | 3 | -0.0322 | 0.8946 | 0.0000 |
| Independent (exploratory grouped) | lora | shuffle | 4 | -0.0107 | 0.8927 | 0.0214 |
| Independent (exploratory grouped) | lora | shuffle | 5 | -0.0322 | 0.8947 | 0.0000 |
| Independent (exploratory grouped) | lora | shuffle | 6 | -0.0536 | 0.8912 | -0.0214 |
| Independent (exploratory grouped) | lora | shuffle | 7 | -0.0322 | 0.8945 | 0.0000 |
| Independent (exploratory grouped) | lora | shuffle | 8 | 0.0107 | 0.8904 | 0.0429 |
| Independent (exploratory grouped) | lora | shuffle | 9 | -0.0536 | 0.8929 | -0.0214 |
| Independent (exploratory grouped) | lora | shuffle | 10 | 0.0107 | 0.8919 | 0.0429 |
| Independent (exploratory grouped) | lora | gaussian | 1 | -0.0107 | 0.8926 | 0.0214 |
| Independent (exploratory grouped) | lora | gaussian | 2 | -0.0322 | 0.8902 | 0.0000 |
| Independent (exploratory grouped) | lora | gaussian | 3 | -0.0322 | 0.8941 | 0.0000 |
| Independent (exploratory grouped) | lora | gaussian | 4 | -0.0322 | 0.8939 | 0.0000 |
| Independent (exploratory grouped) | lora | gaussian | 5 | 0.0107 | 0.8915 | 0.0429 |
| Independent (exploratory grouped) | lora | gaussian | 6 | -0.0322 | 0.8909 | 0.0000 |
| Independent (exploratory grouped) | lora | gaussian | 7 | 0.0036 | 0.8890 | 0.0357 |
| Independent (exploratory grouped) | lora | gaussian | 8 | 0.0214 | 0.8902 | 0.0536 |
| Independent (exploratory grouped) | lora | gaussian | 9 | -0.0536 | 0.8956 | -0.0214 |
| Independent (exploratory grouped) | lora | gaussian | 10 | -0.0536 | 0.8901 | -0.0214 |
| Independent (exploratory grouped) | frozen | shuffle | 1 | -0.0608 | 0.8927 | 0.0000 |
| Independent (exploratory grouped) | frozen | shuffle | 2 | -0.0822 | 0.8971 | -0.0214 |
| Independent (exploratory grouped) | frozen | shuffle | 3 | -0.0661 | 0.8974 | -0.0054 |
| Independent (exploratory grouped) | frozen | shuffle | 4 | -0.0447 | 0.8939 | 0.0161 |
| Independent (exploratory grouped) | frozen | shuffle | 5 | -0.0822 | 0.8994 | -0.0214 |
| Independent (exploratory grouped) | frozen | shuffle | 6 | -0.0751 | 0.8975 | -0.0143 |
| Independent (exploratory grouped) | frozen | shuffle | 7 | -0.0822 | 0.8981 | -0.0214 |
| Independent (exploratory grouped) | frozen | shuffle | 8 | -0.0751 | 0.8982 | -0.0143 |
| Independent (exploratory grouped) | frozen | shuffle | 9 | -0.0751 | 0.8987 | -0.0143 |
| Independent (exploratory grouped) | frozen | shuffle | 10 | -0.0965 | 0.8944 | -0.0357 |
| Independent (exploratory grouped) | frozen | gaussian | 1 | -0.0751 | 0.8980 | -0.0143 |
| Independent (exploratory grouped) | frozen | gaussian | 2 | -0.0125 | 0.8903 | 0.0483 |
| Independent (exploratory grouped) | frozen | gaussian | 3 | -0.0608 | 0.8981 | 0.0000 |
| Independent (exploratory grouped) | frozen | gaussian | 4 | -0.0608 | 0.8962 | 0.0000 |
| Independent (exploratory grouped) | frozen | gaussian | 5 | -0.0894 | 0.8943 | -0.0286 |
| Independent (exploratory grouped) | frozen | gaussian | 6 | -0.0357 | 0.8925 | 0.0250 |
| Independent (exploratory grouped) | frozen | gaussian | 7 | -0.0608 | 0.8968 | 0.0000 |
| Independent (exploratory grouped) | frozen | gaussian | 8 | -0.0751 | 0.8970 | -0.0143 |
| Independent (exploratory grouped) | frozen | gaussian | 9 | -0.0608 | 0.8971 | 0.0000 |
| Independent (exploratory grouped) | frozen | gaussian | 10 | -0.0661 | 0.8902 | -0.0054 |

### Supporting pooled and individual-seed metrics

| Endpoint | Variant | Arm | Pairwise accuracy | Top-5 enrichment | Seed pooled ρ (20261123–27) | Mean | SD | Min | Max |
|---|---|---|---:|---:|---|---:|---:|---:|---:|
| Independent (exploratory grouped) | lora | A | 0.4808 | 0.6000 | -0.0054, 0.0018, 0.0447, 0.0268, -0.1055 | -0.0075 | 0.0521 | -0.1055 | 0.0447 |
| Independent (exploratory grouped) | lora | B | 0.4904 | 0.6000 | 0.0769, -0.0036, 0.0250, 0.0268, -0.0661 | 0.0118 | 0.0468 | -0.0661 | 0.0769 |
| Independent (exploratory grouped) | frozen | A | 0.4904 | 0.6000 | -0.0661, -0.1037, 0.0322, -0.0286, 0.0572 | -0.0218 | 0.0598 | -0.1037 | 0.0572 |
| Independent (exploratory grouped) | frozen | B | 0.4808 | 0.6000 | -0.0661, -0.0947, 0.0322, -0.0286, 0.0250 | -0.0265 | 0.0497 | -0.0947 | 0.0322 |

### Scaffold membership (14 total typed scaffolds)

| Fold | Primary members | Coupled members |
|---|---|---|
| S01 | ENMD-2076 | ENMD-2076 |
| S03 | TAK-901 | TAK-901 |
| S05 | VX-680 | VX-680 |
| S06 | SNS-314 | SNS-314 |
| S07 | AT-9283 | AT-9283 |
| S08 | danusertib | danusertib |
| S10 | KW-2449 | KW-2449 |
| S14 | BI-847325 | BI-847325 |
| S15 | CYC-116 | CYC-116 |
| S17 | PF-03814735 | PF-03814735 |
| S18 | GSK-1070916 | GSK-1070916 |
| S20 | MLN-8054; MLN-8237 | MLN-8054; MLN-8237 |
| S21 | MK-5108 | MK-5108 |
| S22 | AMG-900 | AMG-900 |

### Grouped fold membership

| Fold | Typed scaffold IDs | Held-out compounds |
|---|---|---|
| G01 | S20, S07, S14, S18 | AT-9283; BI-847325; GSK-1070916; MLN-8054; MLN-8237 |
| G02 | S01, S05, S08, S15, S21 | CYC-116; ENMD-2076; MK-5108; VX-680; danusertib |
| G03 | S03, S06, S10, S17, S22 | AMG-900; KW-2449; PF-03814735; SNS-314; TAK-901 |

### Per-fold: Independent (exploratory grouped), lora

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G01 | 5 | -0.1000 | -0.1000 | 0.6206 | 0.6250 | 0.1210 | 0.6005 | 0.2014 |
| G02 | 5 | 0.5000 | 0.5000 | 1.1627 | 1.1609 | 0.0567 | 1.0961 | 0.0517 |
| G03 | 5 | 0.4000 | 0.4000 | 0.8075 | 0.8057 | 0.1093 | 0.6624 | 0.1651 |

### Per-fold: Independent (exploratory grouped), frozen

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G01 | 5 | -0.1000 | -0.1000 | 0.6222 | 0.6220 | 0.1054 | 0.6005 | 0.1755 |
| G02 | 5 | 0.7000 | 0.7000 | 1.1582 | 1.1597 | 0.0424 | 1.0961 | 0.0387 |
| G03 | 5 | 0.3000 | 0.3000 | 0.8270 | 0.8274 | 0.0880 | 0.6624 | 0.1328 |

## Post-hoc exploratory grouped-scaffold sensitivity: circular Lake-coupled

**Post-hoc exploratory sensitivity, added 2026-11-23 UTC.** The sealed leave-one-scaffold-out primary remains the headline regardless of this result.

4 folds with 6, 6, 6, 6 held-out compounds. Four Lake folds achieve exact equality; five cannot evenly divide 24 compounds. Decreasing scaffold size, then canonical scaffold SMILES; assign intact to least-loaded fold, breaking ties by fold ID.

**Sealed primary:** The n=15 panel is underpowered for the 0.3 Δrho bar. A null is uninformative, not evidence of no effect.

Model A is a ligand-only baseline, matching PPARγ: no geometry columns survive. PXD062561 has 0/24 ligand-resolved HDX coverage. Model B therefore adds only training-fold-standardised log10(cooperativity α); no substitutes or predicted complexes were used. SNS-314 remains included despite its poorly constrained α (only a lower bound >10 could be established).

The primary endpoint is pKd from Klaeger et al. (2017), human cell-lysate Kinobeads binding assay CHEMBL3991630, not Lake data. Aurora A is scored on a different endpoint from the other three targets. The Lake no-TPX2 pKd sensitivity is **circular**: it shares the global fit that defined α, and is not headline or independent evidence.

Prespec SHA256: `39fb48b381d695970166c2120299d640af577a3bf03d7306d4727fb23fafc03e`.

The sealed prespec names fold-first execution. The requester asked to finish the primary endpoint before starting the Lake-coupled sensitivity. Only dispatch order across independent endpoint fits changed; the sealed prespec and its SHA256 above are unchanged.

The first primary check stopped with an `IndexError`: a local loop variable shadowed the feature lookup in `check.py`. The checker-only correction did not change any fit or computed metric.

Metrics use seed-mean out-of-fold predictions (pooled across models trained on different subsets). The RMSE scale is pKd. A fold with one compound has undefined rho; a two-compound rho is an ordering sign, not a stable estimate.

| Endpoint | Variant | n / folds | ρ A | ρ B | RMSE A | RMSE B | Δρ | 95% paired CI | Undefined bootstrap | Shuffle / Gaussian at or above Δρ | Decision |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|
| Lake coupled (exploratory grouped) | lora | 24 / 4 | -0.0609 | -0.0496 | 0.9379 | 0.9349 | 0.0113 | [-0.0384, 0.0787] | 0.00% | 3/10 / 1/10 | Sham reproduced the gain; not attributable to cooperativity |
| Lake coupled (exploratory grouped) | frozen | 24 / 4 | -0.0774 | -0.0678 | 0.9352 | 0.9344 | 0.0096 | [0.0000, 0.0439] | 0.00% | 2/10 / 2/10 | Sham reproduced the gain; not attributable to cooperativity |

### Controls

Sham Δρ is relative to the same real Model A. Near-constant means prediction SD ≤10% of label SD on an evaluable fold. Singleton or constant-label folds cannot diagnose near-constancy. Rank tracking means ρ ≥0.8 against the fold-mean predictor.

| Endpoint | Variant | Fold-mean ρ / RMSE | Cooperativity-only ρ / RMSE | A vs fold mean ρ / RMSE | A tracks fold means by rank | Near-constant / evaluable folds |
|---|---|---|---|---|---|---|
| Lake coupled (exploratory grouped) | lora | -0.4092 / 0.9597 | -0.2757 / 0.9571 | 0.5491 / 0.1461 | No | 0/4 |
| Lake coupled (exploratory grouped) | frozen | -0.4092 / 0.9597 | -0.2357 / 0.9621 | 0.6622 / 0.1031 | No | 1/4 |

| Endpoint | Variant | Sham block | Draw | ρ | RMSE | Δρ vs A |
|---|---|---|---:|---:|---:|---:|
| Lake coupled (exploratory grouped) | lora | shuffle | 1 | -0.0661 | 0.9410 | -0.0052 |
| Lake coupled (exploratory grouped) | lora | shuffle | 2 | -0.0600 | 0.9392 | 0.0009 |
| Lake coupled (exploratory grouped) | lora | shuffle | 3 | -0.0296 | 0.9377 | 0.0313 |
| Lake coupled (exploratory grouped) | lora | shuffle | 4 | -0.0522 | 0.9366 | 0.0087 |
| Lake coupled (exploratory grouped) | lora | shuffle | 5 | -0.0600 | 0.9402 | 0.0009 |
| Lake coupled (exploratory grouped) | lora | shuffle | 6 | -0.0678 | 0.9384 | -0.0070 |
| Lake coupled (exploratory grouped) | lora | shuffle | 7 | -0.0843 | 0.9405 | -0.0235 |
| Lake coupled (exploratory grouped) | lora | shuffle | 8 | -0.0383 | 0.9359 | 0.0226 |
| Lake coupled (exploratory grouped) | lora | shuffle | 9 | -0.0678 | 0.9390 | -0.0070 |
| Lake coupled (exploratory grouped) | lora | shuffle | 10 | -0.0426 | 0.9368 | 0.0183 |
| Lake coupled (exploratory grouped) | lora | gaussian | 1 | -0.0661 | 0.9391 | -0.0052 |
| Lake coupled (exploratory grouped) | lora | gaussian | 2 | -0.0757 | 0.9387 | -0.0148 |
| Lake coupled (exploratory grouped) | lora | gaussian | 3 | -0.0678 | 0.9383 | -0.0070 |
| Lake coupled (exploratory grouped) | lora | gaussian | 4 | -0.0678 | 0.9398 | -0.0070 |
| Lake coupled (exploratory grouped) | lora | gaussian | 5 | -0.0704 | 0.9397 | -0.0096 |
| Lake coupled (exploratory grouped) | lora | gaussian | 6 | -0.0730 | 0.9379 | -0.0122 |
| Lake coupled (exploratory grouped) | lora | gaussian | 7 | -0.0339 | 0.9376 | 0.0270 |
| Lake coupled (exploratory grouped) | lora | gaussian | 8 | -0.0626 | 0.9389 | -0.0017 |
| Lake coupled (exploratory grouped) | lora | gaussian | 9 | -0.0670 | 0.9401 | -0.0061 |
| Lake coupled (exploratory grouped) | lora | gaussian | 10 | -0.0765 | 0.9371 | -0.0157 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 1 | -0.0826 | 0.9369 | -0.0052 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 2 | -0.0722 | 0.9362 | 0.0052 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 3 | -0.0826 | 0.9368 | -0.0052 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 4 | -0.0565 | 0.9311 | 0.0209 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 5 | -0.0826 | 0.9374 | -0.0052 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 6 | -0.0713 | 0.9361 | 0.0061 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 7 | -0.0774 | 0.9369 | 0.0000 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 8 | -0.0426 | 0.9306 | 0.0348 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 9 | -0.0687 | 0.9340 | 0.0087 |
| Lake coupled (exploratory grouped) | frozen | shuffle | 10 | -0.0774 | 0.9355 | 0.0000 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 1 | -0.0765 | 0.9365 | 0.0009 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 2 | -0.0835 | 0.9370 | -0.0061 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 3 | -0.0557 | 0.9338 | 0.0217 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 4 | -0.0730 | 0.9345 | 0.0043 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 5 | -0.0713 | 0.9351 | 0.0061 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 6 | -0.0574 | 0.9333 | 0.0200 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 7 | -0.0722 | 0.9360 | 0.0052 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 8 | -0.0774 | 0.9357 | 0.0000 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 9 | -0.0774 | 0.9370 | 0.0000 |
| Lake coupled (exploratory grouped) | frozen | gaussian | 10 | -0.0722 | 0.9352 | 0.0052 |

### Supporting pooled and individual-seed metrics

| Endpoint | Variant | Arm | Pairwise accuracy | Top-5 enrichment | Seed pooled ρ (20261123–27) | Mean | SD | Min | Max |
|---|---|---|---:|---:|---|---:|---:|---:|---:|
| Lake coupled (exploratory grouped) | lora | A | 0.4819 | 0.0000 | 0.0070, -0.0391, -0.0965, -0.0313, -0.1635 | -0.0647 | 0.0595 | -0.1635 | 0.0070 |
| Lake coupled (exploratory grouped) | lora | B | 0.4819 | 0.0000 | 0.0104, -0.0304, -0.0652, -0.0061, -0.1557 | -0.0494 | 0.0589 | -0.1557 | 0.0104 |
| Lake coupled (exploratory grouped) | frozen | A | 0.4746 | 0.0000 | -0.0696, -0.1461, -0.1617, 0.0035, -0.1209 | -0.0990 | 0.0600 | -0.1617 | 0.0035 |
| Lake coupled (exploratory grouped) | frozen | B | 0.4819 | 0.0000 | -0.0696, -0.1704, -0.1165, 0.0183, -0.1000 | -0.0877 | 0.0623 | -0.1704 | 0.0183 |

### Scaffold membership (23 total typed scaffolds)

| Fold | Primary members | Coupled members |
|---|---|---|
| S01 | ENMD-2076 | ENMD-2076 |
| S02 | — | CCT-129202 |
| S03 | TAK-901 | TAK-901 |
| S04 | — | ZM-447439 |
| S05 | VX-680 | VX-680 |
| S06 | SNS-314 | SNS-314 |
| S07 | AT-9283 | AT-9283 |
| S08 | danusertib | danusertib |
| S09 | — | PHA-680632 |
| S10 | KW-2449 | KW-2449 |
| S11 | — | MK-8745 |
| S12 | — | JNJ-7706621 |
| S13 | — | hesperadin |
| S14 | BI-847325 | BI-847325 |
| S15 | CYC-116 | CYC-116 |
| S16 | — | CCT-137690 |
| S17 | PF-03814735 | PF-03814735 |
| S18 | GSK-1070916 | GSK-1070916 |
| S19 | — | TAE-684 |
| S20 | MLN-8054; MLN-8237 | MLN-8054; MLN-8237 |
| S21 | MK-5108 | MK-5108 |
| S22 | AMG-900 | AMG-900 |
| S23 | — | reversine |

### Grouped fold membership

| Fold | Typed scaffold IDs | Held-out compounds |
|---|---|---|
| G01 | S20, S07, S11, S15, S19 | AT-9283; CYC-116; MK-8745; MLN-8054; MLN-8237; TAE-684 |
| G02 | S01, S04, S08, S12, S16, S21 | CCT-137690; ENMD-2076; JNJ-7706621; MK-5108; ZM-447439; danusertib |
| G03 | S02, S05, S09, S13, S17, S22 | AMG-900; CCT-129202; PF-03814735; PHA-680632; VX-680; hesperadin |
| G04 | S03, S06, S10, S14, S18, S23 | BI-847325; GSK-1070916; KW-2449; SNS-314; TAK-901; reversine |

### Per-fold: Lake coupled (exploratory grouped), lora

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G01 | 6 | 0.6000 | 0.6000 | 0.8074 | 0.7999 | 0.1360 | 0.6967 | 0.1952 |
| G02 | 6 | 0.5429 | 0.5429 | 1.2953 | 1.2920 | 0.1374 | 1.3162 | 0.1044 |
| G03 | 6 | -0.1429 | -0.1429 | 0.6689 | 0.6610 | 0.1479 | 0.5371 | 0.2754 |
| G04 | 6 | 0.6571 | 0.6000 | 0.8613 | 0.8662 | 0.1369 | 0.5541 | 0.2470 |

### Per-fold: Lake coupled (exploratory grouped), frozen

| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G01 | 6 | 0.6000 | 0.6000 | 0.8296 | 0.8275 | 0.1183 | 0.6967 | 0.1698 |
| G02 | 6 | 0.5429 | 0.5429 | 1.2940 | 1.2929 | 0.0867 | 1.3162 | 0.0659 |
| G03 | 6 | 0.0286 | 0.0286 | 0.5886 | 0.5864 | 0.0883 | 0.5371 | 0.1644 |
| G04 | 6 | 0.6000 | 0.6571 | 0.8886 | 0.8898 | 0.1070 | 0.5541 | 0.1932 |
