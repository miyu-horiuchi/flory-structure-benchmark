from statistics import mean, pstdev

import numpy as np
from scipy.stats import rankdata

from common import (
    ARMS, BOOTSTRAP_N, GROUPED_OUTPUTS, GROUPED_PLANS, LABELS, RESULTS, RUN, SEEDS, SHAMS, VARIANTS,
    read_json, read_prepared, run_scope, sha256, write_json,
)


def rank_correlations(x, y):
    x_rank, y_rank = rankdata(x, axis=-1), rankdata(y, axis=-1)
    x_center = x_rank - x_rank.mean(axis=-1, keepdims=True)
    y_center = y_rank - y_rank.mean(axis=-1, keepdims=True)
    numer = np.sum(x_center * y_center, axis=-1)
    denom = np.sqrt(np.sum(x_center ** 2, axis=-1) * np.sum(y_center ** 2, axis=-1))
    return np.divide(numer, denom, out=np.full(numer.shape, np.nan), where=denom != 0)


def spearman(y, p):
    rho = float(rank_correlations(y, p)) if len(y) >= 2 else np.nan
    return rho if np.isfinite(rho) else None


def paired_bootstrap(y, a, b):
    n = len(y)
    samples = np.random.default_rng(SEEDS[0]).integers(0, n, size=(BOOTSTRAP_N, n))
    delta = rank_correlations(y[samples], b[samples]) - rank_correlations(y[samples], a[samples])
    finite = delta[np.isfinite(delta)]
    ci = np.percentile(finite, [2.5, 97.5]).tolist() if len(finite) else [None, None]
    return {
        "ci95": ci, "undefined_fraction": (BOOTSTRAP_N - len(finite)) / BOOTSTRAP_N,
        "resamples": BOOTSTRAP_N, "seed": SEEDS[0],
    }


def pairwise_accuracy(y, p):
    i, j = np.triu_indices(len(y), k=1)
    true, pred = np.sign(y[i] - y[j]), np.sign(p[i] - p[j])
    valid = true != 0
    return float(np.mean((true[valid] == pred[valid]) + 0.5 * (pred[valid] == 0))) if valid.any() else None


def top_k_membership(values, k):
    cutoff = np.sort(values)[-k]
    above, tied = values > cutoff, values == cutoff
    return above + tied * ((k - above.sum()) / tied.sum())


def top5_enrichment(y, p):
    k = 5
    if len(y) < k:
        return None
    overlap = np.dot(top_k_membership(y, k), top_k_membership(p, k))
    return float(overlap * len(y) / k ** 2)


def metrics(y, p, *, pooled=False):
    result = {
        "spearman": spearman(y, p), "rmse_pkd": float(np.sqrt(np.mean((y - p) ** 2))),
        "pairwise_accuracy": pairwise_accuracy(y, p),
    }
    if pooled:
        result["top5_enrichment"] = top5_enrichment(y, p)
    return result


def seed_summary(values):
    defined = [rho for rho in values if rho is not None]
    funcs = {"mean": mean, "sd": pstdev, "min": min, "max": max}
    return {key: func(defined) if defined else None for key, func in funcs.items()} | {"n_defined": len(defined)}


def reconstruct(prepared, training, label, variant, arm, seed, members):
    folds = {f["id"]: f["members"][label] for f in prepared["folds"] if f["members"][label]}
    predictions = {}
    tested_folds = set()
    for run in training["runs"]:
        if (run["label"], run["variant"], run["arm"], run["seed"]) != (label, variant, arm, seed):
            continue
        assert run["fold"] not in tested_folds and set(run["oof"]) == set(folds[run["fold"]])
        assert not set(predictions).intersection(run["oof"])
        tested_folds.add(run["fold"])
        predictions.update(run["oof"])
    assert tested_folds == set(folds) and set(predictions) == set(members)
    return [predictions[name] for name in members]


def mean_oof(prepared, training, label, variant, arm, members):
    predictions = np.array([
        reconstruct(prepared, training, label, variant, arm, seed, members) for seed in SEEDS
    ], dtype=np.float64)
    return predictions, predictions.mean(axis=0)


def control_metrics(prepared, training, label, variant, members, folds, y, a, b, rho_a, delta):
    truth = prepared["labels"][label]
    fold_diagnostics = []
    fold_means = np.empty_like(y)
    index = {name: i for i, name in enumerate(members)}
    for fold in folds:
        held = fold["members"][label]
        indices = [index[name] for name in held]
        fold_means[indices] = mean(truth[name] for name in members if name not in held)
        label_sd, pred_sd = float(y[indices].std()), float(a[indices].std())
        ratio = pred_sd / label_sd if label_sd else None
        fold_diagnostics.append({
            "fold": fold["id"], "n": len(held), "label_sd_pkd": label_sd, "prediction_sd_pkd": pred_sd,
            "sd_ratio": ratio, "near_constant": None if ratio is None else ratio <= 0.1,
        })
    tracks = metrics(fold_means, a)
    track_rho = tracks["spearman"]
    tracks["tracks_rank"] = None if track_rho is None else track_rho >= 0.8
    _, features_only = mean_oof(prepared, training, label, variant, "features_only", members)
    controls = {
        "fold_mean": metrics(y, fold_means, pooled=True),
        "features_only": metrics(y, features_only, pooled=True),
        "shams": {},
        "model_a_diagnostics": {"agreement_with_fold_mean": tracks, "per_fold": fold_diagnostics},
        "model_b_agreement_with_fold_mean": metrics(fold_means, b),
    }
    for block, arms in SHAMS.items():
        draws = []
        for i, arm in enumerate(arms, 1):
            _, prediction = mean_oof(prepared, training, label, variant, arm, members)
            point = metrics(y, prediction, pooled=True)
            sham_rho = point["spearman"]
            sham_delta = None if rho_a is None or sham_rho is None else sham_rho - rho_a
            draws.append({"draw": i, **point, "delta_spearman": sham_delta,
                          "at_or_above_real_delta": None if sham_delta is None or delta is None else sham_delta >= delta})
        controls["shams"][block] = {
            "draws": draws,
            "n_at_or_above_real_delta": (None if delta is None else
                                         sum(draw["at_or_above_real_delta"] is True for draw in draws)),
            "n_undefined_deltas": sum(draw["delta_spearman"] is None for draw in draws),
        }
    return controls


def interpretation(label, delta, clears, controls):
    reproduces = delta is not None and delta > 0 and any(
        block["n_at_or_above_real_delta"] for block in controls["shams"].values()
    )
    if reproduces:
        return "Sham reproduced the gain; not attributable to cooperativity"
    if not clears:
        return "Null/sub-threshold; uninformative"
    if label == "coupled":
        return "Clears numerical bar; circular, no independent claim"
    return "Clears numerical bar; tested shams did not reproduce the gain"


def evaluate(prepared, training, labels=LABELS, suffix=""):
    assert training["prespec_sha256"] == prepared["prespec_sha256"]
    assert training["prepared_sha256"] == sha256(RUN / "prepared.json")
    scaffolds = prepared.get("scaffolds", prepared["folds"])
    results = {
        "prespec_sha256": prepared["prespec_sha256"], "source_manifest_sha256": prepared["source_manifest_sha256"],
        "training_sha256": sha256(RUN / f"training{suffix}.json"),
        "model_revision": training["model"]["revision"], "seed_order": SEEDS,
        "geometry_columns": prepared["geometry_columns"], "model_b_dynamics_columns": prepared["dynamics_columns"],
        "scaffold_count": len(scaffolds), "scaffolds": scaffolds,
        "hdx": {"accession": "PXD062561", "covered": 0, "n": 24},
        "sns314": "Retained; mean alpha 122 poorly constrained, lower bound alpha > 10 (Lake SI Fig. S10).",
        "labels": {},
    }
    if "exploratory" in prepared:
        results["exploratory"] = prepared["exploratory"]
        results["grouped_folds"] = prepared["folds"]
    for label in labels:
        truth = prepared["labels"][label]
        members = sorted(truth)
        y = np.array([truth[name] for name in members])
        folds = [fold for fold in prepared["folds"] if fold["members"][label]]
        entry = {
            "n": len(members), "fold_count": len(folds), "variants": {},
            "label": "pKd = 9 - log10(Kd_nM); higher is stronger binding",
            "endpoint": ("CHEMBL3991630 (Klaeger et al., 2017)" if label == "primary"
                         else "Lake et al. (2018), Fig. S9 no-TPX2 Kd"),
            "independent": label == "primary",
            "caveat": ("n=15 is underpowered; a null is uninformative, not evidence of no effect."
                       if label == "primary" else
                       "Circular sensitivity: pKd was obtained in the same global fit used to derive alpha."),
        }
        for variant in VARIANTS:
            arms, pooled = {}, {}
            for arm in ARMS:
                predictions, pooled[arm] = mean_oof(prepared, training, label, variant, arm, members)
                seed_points = [{"seed": seed, **metrics(y, point, pooled=True)}
                               for seed, point in zip(SEEDS, predictions)]
                index = {name: j for j, name in enumerate(members)}
                by_fold = []
                for fold in folds:
                    names = fold["members"][label]
                    indices = [index[name] for name in names]
                    by_fold.append({"fold": fold["id"], "members": names, "n": len(names),
                                    **metrics(y[indices], pooled[arm][indices])})
                arms[arm] = {
                    "pooled": metrics(y, pooled[arm], pooled=True), "per_fold": by_fold,
                    "individual_seeds": seed_points,
                    "seed_spearman_summary": seed_summary([point["spearman"] for point in seed_points]),
                    "pooled_oof_predictions": dict(zip(members, pooled[arm].tolist())),
                }
            a, b = (arms[arm]["pooled"]["spearman"] for arm in ARMS)
            delta = None if a is None or b is None else b - a
            boot = paired_bootstrap(y, pooled["A"], pooled["B"])
            lower = boot["ci95"][0]
            clears = delta is not None and delta >= 0.3 and lower is not None and lower > 0
            controls = control_metrics(prepared, training, label, variant, members, folds, y,
                                       pooled["A"], pooled["B"], a, delta)
            entry["variants"][variant] = {
                "arms": arms, "delta_spearman": delta, "bootstrap": boot, "clears_numerical_bar": clears,
                "controls": controls, "decision": interpretation(label, delta, clears, controls),
            }
        results["labels"][label] = entry
    if labels == LABELS:
        sensitivities = {}
        for label, directory in GROUPED_OUTPUTS.items():
            grouped = RESULTS / directory
            if not (grouped / "checks.json").is_file():
                continue
            check = read_json(grouped / "checks.json")
            grouped_metrics = grouped / "metrics.json"
            assert check["status"] == "passed" and check["metric_results_sha256"] == sha256(grouped_metrics)
            sensitivity = read_json(grouped_metrics)
            for key in ("prespec_sha256", "source_manifest_sha256", "model_revision"):
                assert sensitivity[key] == results[key]
            assert list(sensitivity["labels"]) == [label]
            assert sensitivity["exploratory"]["plan_sha256"] == sha256(GROUPED_PLANS[label])
            sensitivities[label] = sensitivity
        if sensitivities:
            results["exploratory_grouped"] = sensitivities
    return results


def fmt(x):
    return "undefined" if x is None else f"{x:.4f}"


def results_markdown(results):
    primary = results["labels"].get("primary")
    exploratory = results.get("exploratory")
    power = []
    headline = (f"**Primary LoRA interpretation:** {primary['variants']['lora']['decision']}. "
                "The numerical bar is necessary, not sufficient."
                if primary else "**Lake-coupled sensitivity only.** This report is not independent evidence.")
    if exploratory:
        headline = (f"**Post-hoc exploratory sensitivity, added {exploratory['added_utc']} UTC.** "
                    "The sealed leave-one-scaffold-out primary remains the headline regardless of this result.")
        sizes = [len(next(iter(fold["members"].values()))) for fold in results["grouped_folds"]]
        baseline = (f"{len(sizes)} folds with {', '.join(map(str, sizes))} held-out compounds. "
                    + ("Four Lake folds achieve exact equality; five cannot evenly divide 24 compounds. " if not primary else "")
                    + exploratory["assignment"])
    elif primary:
        lora = primary["variants"]["lora"]
        frozen = primary["variants"]["frozen"]
        lora_controls, frozen_controls = (cell["controls"] for cell in (lora, frozen))
        if not lora["clears_numerical_bar"]:
            headline = (f"**The sealed decision rule did not resolve a Model B benefit at n={primary['n']}.** "
                        f"The primary LoRA Δρ was {fmt(lora['delta_spearman'])}; the rule requires Δρ ≥0.3 and a "
                        "positive lower 95% confidence bound. This result is not evidence against cooperativity or dynamics.")
        singleton_count = sum(fold["n"] == 1 for fold in lora["arms"]["A"]["per_fold"])
        fold_mean = lora_controls["fold_mean"]
        baseline = (f"{singleton_count} of {primary['fold_count']} occupied primary folds contain one compound. "
                    f"The fold-mean-only predictor scores ρ {fmt(fold_mean['spearman'])}: on singleton folds, "
                    "holding out a higher pKd makes the training mean lower by construction. This arithmetic "
                    "property is expected and does not diagnose model failure or anti-prediction.")
        correlations = [fmt(c["model_a_diagnostics"]["agreement_with_fold_mean"]["spearman"])
                        for c in (lora_controls, frozen_controls)]
        sham_tallies = []
        for controls in (lora_controls, frozen_controls):
            counts = [controls["shams"][block]["n_at_or_above_real_delta"] for block in SHAMS]
            rates = ["undefined" if count is None else f"{count}/10" for count in counts]
            sham_tallies.append(f"{rates[0]} shuffled and {rates[1]} Gaussian")
        power = [
            f"Model A’s rank correlation with the training-fold mean is {correlations[0]} for LoRA and "
            f"{correlations[1]} for frozen. Sham draws at or above the real Δρ: {sham_tallies[0]} for LoRA; "
            f"{sham_tallies[1]} for frozen."
            + (" The positive LoRA gain is not attributable to cooperativity."
               if lora["decision"].startswith("Sham reproduced") else ""),
            "",
            "In the protocol owner’s power analysis, the unchanged decision rule fired 38% of the time at n=26 "
            "with a planted true Δρ of exactly 0.3. Regrouping folds changed this by less than one percentage point. "
            "The CI-excludes-zero condition alone had 31% power at n=26 and needed about 104 compounds for 80% power. "
            "At n=15 the sealed result is unresolved; grouped k-fold is a comparability analysis, not a rescue.",
            "",
        ]
    else:
        baseline = "A singleton fold has no within-fold rank correlation."
    lines = [
        "# Aurora A ChemBERTa A/B",
        "",
        headline,
        "",
        baseline,
        "",
        *power,
        "**Sealed primary:** The n=15 panel is underpowered for the 0.3 Δrho bar. A null is uninformative, not evidence of no effect.",
        "",
        "Model A is a ligand-only baseline, matching PPARγ: no geometry columns survive. PXD062561 has 0/24 ligand-resolved HDX coverage. Model B therefore adds only training-fold-standardised log10(cooperativity α); no substitutes or predicted complexes were used. SNS-314 remains included despite its poorly constrained α (only a lower bound >10 could be established).",
        "",
        "The primary endpoint is pKd from Klaeger et al. (2017), human cell-lysate Kinobeads binding assay CHEMBL3991630, not Lake data. Aurora A is scored on a different endpoint from the other three targets. The Lake no-TPX2 pKd sensitivity is **circular**: it shares the global fit that defined α, and is not headline or independent evidence.",
        "",
        f"Prespec SHA256: `{results['prespec_sha256']}`.",
        "",
        "The sealed prespec names fold-first execution. The requester asked to finish the primary endpoint before starting the Lake-coupled sensitivity. Only dispatch order across independent endpoint fits changed; the sealed prespec and its SHA256 above are unchanged.",
        "",
        "The first primary check stopped with an `IndexError`: a local loop variable shadowed the feature lookup in `check.py`. The checker-only correction did not change any fit or computed metric.",
        "",
        "Metrics use seed-mean out-of-fold predictions (pooled across models trained on different subsets). The RMSE scale is pKd. A fold with one compound has undefined rho; a two-compound rho is an ordering sign, not a stable estimate.",
        "",
        "| Endpoint | Variant | n / folds | ρ A | ρ B | RMSE A | RMSE B | Δρ | 95% paired CI | Undefined bootstrap | Shuffle / Gaussian at or above Δρ | Decision |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|",
    ]
    names = {"primary": "Independent (exploratory grouped)" if exploratory else "Independent",
             "coupled": "Lake coupled (exploratory grouped)" if exploratory else "Lake coupled (sensitivity)"}
    for label, entry in results["labels"].items():
        for variant, cell in entry["variants"].items():
            a, b = (cell["arms"][arm]["pooled"] for arm in ARMS)
            lo, hi = cell["bootstrap"]["ci95"]
            counts = [block["n_at_or_above_real_delta"] for block in cell["controls"]["shams"].values()]
            tally = " / ".join("undefined" if count is None else f"{count}/10" for count in counts)
            lines.append(f"| {names[label]} | {variant} | {entry['n']} / {entry['fold_count']} | "
                         f"{fmt(a['spearman'])} | {fmt(b['spearman'])} | {fmt(a['rmse_pkd'])} | {fmt(b['rmse_pkd'])} | "
                         f"{fmt(cell['delta_spearman'])} | [{fmt(lo)}, {fmt(hi)}] | "
                         f"{cell['bootstrap']['undefined_fraction']:.2%} | {tally} | {cell['decision']} |")
    lines += ["", "## Controls", "",
              "Sham Δρ is relative to the same real Model A. Near-constant means prediction SD ≤10% of label SD on an evaluable fold. Singleton or constant-label folds cannot diagnose near-constancy. Rank tracking means ρ ≥0.8 against the fold-mean predictor.",
              "",
              "| Endpoint | Variant | Fold-mean ρ / RMSE | Cooperativity-only ρ / RMSE | A vs fold mean ρ / RMSE | A tracks fold means by rank | Near-constant / evaluable folds |",
              "|---|---|---|---|---|---|---|"]
    for label, entry in results["labels"].items():
        for variant, cell in entry["variants"].items():
            controls = cell["controls"]
            diag = controls["model_a_diagnostics"]
            compare = diag["agreement_with_fold_mean"]
            tracking = {True: "Yes", False: "No", None: "Undefined"}[compare["tracks_rank"]]
            near = [fold["near_constant"] for fold in diag["per_fold"] if fold["near_constant"] is not None]
            pairs = " | ".join(f"{fmt(p['spearman'])} / {fmt(p['rmse_pkd'])}"
                               for p in (controls["fold_mean"], controls["features_only"], compare))
            lines.append(f"| {names[label]} | {variant} | {pairs} | {tracking} | {sum(near)}/{len(near)} |")
    lines += ["", "| Endpoint | Variant | Sham block | Draw | ρ | RMSE | Δρ vs A |",
              "|---|---|---|---:|---:|---:|---:|"]
    for label, entry in results["labels"].items():
        for variant, cell in entry["variants"].items():
            for block, summary in cell["controls"]["shams"].items():
                for point in summary["draws"]:
                    lines.append(f"| {names[label]} | {variant} | {block} | {point['draw']} | "
                                 f"{fmt(point['spearman'])} | {fmt(point['rmse_pkd'])} | {fmt(point['delta_spearman'])} |")
    lines += ["", "## Supporting pooled and individual-seed metrics", "",
              "| Endpoint | Variant | Arm | Pairwise accuracy | Top-5 enrichment | Seed pooled ρ (20261123–27) | Mean | SD | Min | Max |",
              "|---|---|---|---:|---:|---|---:|---:|---:|---:|"]
    for label, entry in results["labels"].items():
        for variant, cell in entry["variants"].items():
            for arm, result in cell["arms"].items():
                pooled, summary = result["pooled"], result["seed_spearman_summary"]
                seeds = ", ".join(fmt(point["spearman"]) for point in result["individual_seeds"])
                stats = " | ".join(fmt(summary[name]) for name in ("mean", "sd", "min", "max"))
                lines.append(f"| {names[label]} | {variant} | {arm} | {fmt(pooled['pairwise_accuracy'])} | "
                             f"{fmt(pooled['top5_enrichment'])} | {seeds} | {stats} |")
    lines += ["", f"## Scaffold membership ({results['scaffold_count']} total typed scaffolds)", "",
              "| Fold | Primary members | Coupled members |", "|---|---|---|"]
    for fold in results["scaffolds"]:
        primary, coupled = ("; ".join(fold["members"].get(label, [])) or "—" for label in LABELS)
        lines.append(f"| {fold['id']} | {primary} | {coupled} |")
    if exploratory:
        lines += ["", "## Grouped fold membership", "",
                  "| Fold | Typed scaffold IDs | Held-out compounds |", "|---|---|---|"]
        for fold in results["grouped_folds"]:
            groups = ", ".join(fold["scaffold_ids"])
            members = "; ".join(next(iter(fold["members"].values())))
            lines.append(f"| {fold['id']} | {groups} | {members} |")
    for label, entry in results["labels"].items():
        for variant, cell in entry["variants"].items():
            lines += ["", f"## Per-fold: {names[label]}, {variant}", "",
                      "| Fold | n | ρ A | ρ B | RMSE A | RMSE B | A prediction SD | pKd label SD | SD ratio |",
                      "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
            diagnostics = cell["controls"]["model_a_diagnostics"]["per_fold"]
            for a, b, diag in zip(*(cell["arms"][arm]["per_fold"] for arm in ARMS), diagnostics):
                assert a["fold"] == b["fold"] == diag["fold"]
                lines.append(f"| {a['fold']} | {a['n']} | {fmt(a['spearman'])} | {fmt(b['spearman'])} | "
                             f"{fmt(a['rmse_pkd'])} | {fmt(b['rmse_pkd'])} | {fmt(diag['prediction_sd_pkd'])} | "
                             f"{fmt(diag['label_sd_pkd'])} | {fmt(diag['sd_ratio'])} |")
    if "exploratory_grouped" in results:
        for label, sensitivity in results["exploratory_grouped"].items():
            title = "independent primary" if label == "primary" else "circular Lake-coupled"
            lines += ["", f"## Post-hoc exploratory grouped-scaffold sensitivity: {title}", ""]
            lines += [f"#{line}" if line.startswith("## ") else line
                      for line in results_markdown(sensitivity).splitlines()[2:]]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    labels, suffix = run_scope()
    result = evaluate(read_prepared(), read_json(RUN / f"training{suffix}.json"), labels, suffix)
    write_json(RESULTS / f"metrics{suffix}.json", result)
    (RESULTS / f"results{suffix}.md").write_text(results_markdown(result))
    for label, entry in result["labels"].items():
        for variant, cell in entry["variants"].items():
            print(f"{label}/{variant}: delta={fmt(cell['delta_spearman'])}, "
                  f"CI={[fmt(x) for x in cell['bootstrap']['ci95']]}, {cell['decision']}")
