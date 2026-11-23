import hashlib
import math
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import torch
from scipy.stats import spearmanr

from common import (
    ARMS, BATCH, DATA, EPOCHS, FIT_ARMS, LABELS, MODEL_MANIFEST, RESULTS, REVISION, RUN, SEEDS, SHAMS,
    SNAPSHOT, VARIANTS,
    prespec_sha, read_json, read_prepared, run_scope, sha256, verify_manifest, write_json,
)
from evaluate import (
    evaluate, interpretation, paired_bootstrap, pairwise_accuracy, results_markdown, spearman, top5_enrichment,
)
from prepare import prepare
from train import configure_cpu, feature_blocks, load_encoder, partition, tensor_digest, train_cell


def check_sources(prepared, training):
    assert prespec_sha() == prepared["prespec_sha256"] == training["prespec_sha256"]
    assert training["prepared_sha256"] == sha256(RUN / "prepared.json")
    assert prepared == prepare()
    expected = {path.relative_to(DATA).as_posix() for path in DATA.iterdir()
                if path.is_file() and path.name != "SHA256SUMS"}
    assert expected == {name.removeprefix("./") for name in prepared["source_sha256"]}
    print(f"PASS: sealed prespec and all {len(expected)} source files verified; prepared artifacts reproducible")


def check_model(prepared, training):
    report = training["model"]
    assert report["revision"] == REVISION and report["dtype"] == "torch.float32"
    assert report["hidden_size"] == 768 and report["num_hidden_layers"] == 6
    assert report["file_sha256"] == verify_manifest(MODEL_MANIFEST, SNAPSHOT)
    assert not any(report["loading"][k] for k in ("missing_keys", "mismatched_keys", "error_msgs"))
    workspace, loaded = load_encoder(prepared["compounds"])
    assert loaded == report
    assert all(run["base_sha256"] == report["base_sha256"] for run in training["runs"])
    print(f"PASS: pinned checkpoint and reloaded base weights match all fits; missing={loaded['loading']['missing_keys']}")
    return workspace


def keys(prepared, labels):
    return [(f["id"], label, variant, arm, seed)
            for f in prepared["folds"] for label in labels if f["members"][label]
            for variant in VARIANTS for arm in FIT_ARMS for seed in SEEDS]


def scaler(values):
    mean = math.fsum(values) / len(values)
    sd = math.sqrt(math.fsum((v - mean) ** 2 for v in values) / len(values))
    return {"mean": mean, "sd": sd or 1.0}


def assert_scalers(actual, expected):
    for block, stats in expected.items():
        for name, value in stats.items():
            assert math.isclose(actual[block][name], value, rel_tol=1e-14, abs_tol=1e-14), (actual, expected)


def input_sha(fold_data):
    return tensor_digest([
        ("targets", fold_data.train_pytorch_y), ("dynamics", fold_data.train_pytorch_feature),
    ])


def check_exact_shortcuts(prepared, workspace, label, cells):
    fold = next(f for f in prepared["folds"] if f["members"][label])
    truth = prepared["labels"][label]
    held = fold["members"][label]
    features = feature_blocks(prepared["compounds"], truth, label)
    uncached = replace(workspace, frozen={})
    max_error = 0.0
    for arm in ARMS:
        data = partition(features[arm], truth, held, arm)
        for seed in SEEDS:
            saved = cells[("frozen", arm, seed)]
            replay = train_cell(uncached, data, "frozen", seed, arm)
            for field in ("init_sha256", "batch_sha256", "steps"):
                assert saved[field] == replay[field]
            max_error = max(max_error, *(abs(saved["oof"][name] - replay["oof"][name]) for name in held))
    if max_error >= 1e-6:
        raise AssertionError(f"{label} frozen embedding cache: max prediction error {max_error:.9g} pKd >= 1e-6")
    data = partition(features["features_only"], truth, held, "features_only")
    for variant in VARIANTS:
        saved = cells[(variant, "features_only", SEEDS[0])]
        with patch.object(uncached, "embed", side_effect=lambda _, ids: torch.zeros((len(ids), 768))):
            replay = train_cell(uncached, data, variant, SEEDS[0], "B")
        for field in ("init_sha256", "final_sha256", "batch_sha256", "steps", "oof"):
            assert saved[field] == replay[field], (label, variant, field)
    print(f"PASS: {label} frozen A/B: two completed cells, all five seeds, max cached-vs-uncached error "
          f"{max_error:.9g} pKd; encoder-free cooperativity-only fits matched exactly in both variants", flush=True)


def check_training(prepared, training, workspace, labels_selected):
    names = ("fold", "label", "variant", "arm", "seed")
    runs = {tuple(record[name] for name in names): record for record in training["runs"]}
    expected_keys = keys(prepared, labels_selected)
    assert list(runs) == expected_keys and len(runs) == len(training["runs"]) == len(expected_keys)
    compounds = prepared["compounds"]
    blocks = {label: feature_blocks(compounds, prepared["labels"][label], label) for label in labels_selected}
    for label, features in blocks.items():
        names = sorted(prepared["labels"][label])
        values = [compounds[name]["log10_alpha"] for name in names]
        for arm in SHAMS["shuffle"]:
            assert sorted(features[arm].values()) == sorted(values)
        for group in SHAMS.values():
            assert len({tuple(features[arm].values()) for arm in group}) == 10
    replayed = 0
    replayed_labels = set()
    for fold in prepared["folds"]:
        for label in labels_selected:
            held = fold["members"][label]
            if not held:
                continue
            labels = prepared["labels"][label]
            train_ids = sorted(set(labels) - set(held))
            assert not set(held).intersection(train_ids)
            assert set(held).union(train_ids) == set(labels)
            for scaffold in prepared.get("scaffolds", prepared["folds"]):
                members = set(scaffold["members"][label])
                assert not (members.intersection(train_ids) and members.intersection(held))
            poisoned_labels = labels | {name: 100000.0 + i for i, name in enumerate(held)}
            expected_steps = EPOCHS * math.ceil(len(train_ids) / BATCH)
            for arm in FIT_ARMS:
                features = blocks[label][arm]
                expected = {"label": scaler([labels[name] for name in train_ids]),
                            "dynamics": scaler([features[name] for name in train_ids])}
                fold_data = partition(features, labels, held, arm)
                poisoned = partition(features, poisoned_labels, held, arm)
                assert_scalers(fold_data.scalers, expected)
                assert_scalers(poisoned.scalers, expected)
                assert fold_data.train_ids == poisoned.train_ids == tuple(train_ids)
                assert fold_data.test_ids == poisoned.test_ids == tuple(held)
                assert input_sha(fold_data) == input_sha(poisoned)
                assert torch.equal(fold_data.test_pytorch_feature, poisoned.test_pytorch_feature)
                for ids, tensor in ((train_ids, fold_data.train_pytorch_feature),
                                    (held, fold_data.test_pytorch_feature)):
                    values = np.array([features[name] for name in ids])
                    stats = expected["dynamics"]
                    expected_values = np.zeros_like(values) if arm == "A" else (values - stats["mean"]) / stats["sd"]
                    assert np.allclose(tensor.numpy().ravel(), expected_values, rtol=1e-6, atol=1e-7)
                    if arm == "A":
                        assert torch.count_nonzero(tensor).item() == 0

                for variant in VARIANTS:
                    for seed in SEEDS:
                        r = runs[fold["id"], label, variant, arm, seed]
                        assert r["train_ids"] == train_ids and r["test_ids"] == held
                        assert_scalers(r["scalers"], expected)
                        assert r["steps"] == expected_steps
                        assert r["training_inputs_sha256"] == input_sha(fold_data)
                        assert list(r["oof"]) == held and all(math.isfinite(p) for p in r["oof"].values())
                        expected_head_params = (769 * 32 + 32) + (32 + 1)
                        expected_adapters = 6 * 2 * 2 * 768 * 8 if variant == "lora" else 0
                        assert r["n_trainable_parameters"] == expected_head_params + expected_adapters
                        assert r["n_trainable_tensors"] == 4 + (24 if variant == "lora" else 0)
                        order = torch.Generator().manual_seed(seed)
                        batch_hash = hashlib.sha256()
                        for _ in range(EPOCHS):
                            batch_hash.update(torch.randperm(len(train_ids), generator=order).numpy().tobytes())
                        assert r["batch_sha256"] == batch_hash.hexdigest()
                    if label not in replayed_labels and arm in (*ARMS, "features_only", *(v[0] for v in SHAMS.values())):
                        seed = SEEDS[0]
                        r = runs[fold["id"], label, variant, arm, seed]
                        replay = train_cell(workspace, poisoned, variant, seed, arm)
                        for field in ("init_sha256", "final_sha256", "batch_sha256", "steps", "oof"):
                            assert r[field] == replay[field], (fold["id"], label, variant, arm, field)
                        replayed += 1
            replayed_labels.add(label)
            for variant in VARIANTS:
                for seed in SEEDS:
                    a, b = (runs[fold["id"], label, variant, arm, seed] for arm in ARMS)
                    assert a["scalers"] == b["scalers"]
                    for arm in FIT_ARMS[1:]:
                        other = runs[fold["id"], label, variant, arm, seed]
                        for field in ("init_sha256", "batch_sha256", "steps", "n_trainable_parameters",
                                      "n_trainable_tensors"):
                            assert a[field] == other[field], (fold["id"], label, variant, seed, arm, field)
    assert replayed == 10 * len(labels_selected)
    for label in labels_selected:
        fold = next(f["id"] for f in prepared["folds"] if f["members"][label])
        cells = {(variant, arm, seed): runs[fold, label, variant, arm, seed]
                 for variant in VARIANTS for arm in (*ARMS, "features_only") for seed in SEEDS}
        check_exact_shortcuts(prepared, workspace, label, cells)
    print(f"PASS: verified fold isolation, train-fold-only scalers, exact-zero A, and OOF coverage for all {len(runs)} fits")
    print("PASS: identical initialisation, parameter counts, batch order and optimiser steps across matched arms")
    count = sum(bool(fold["members"][label]) for fold in prepared["folds"] for label in labels_selected)
    print(f"PASS: poisoned every held-out label in all {count} fold/endpoint cells; {replayed} full deterministic replays unchanged")
    return replayed


def check_metrics(prepared, training, labels=LABELS, suffix=""):
    y = np.array([1.0, 2.0, 2.0, 4.0, 5.0, 6.0])
    p = np.array([1.0, 4.0, 2.0, 2.0, 6.0, 5.0])
    c = np.ones(len(y))
    assert spearman(y[:1], p[:1]) is None and spearman(y, c) is None and spearman(c, y) is None
    assert np.isclose(spearman(y, p), spearmanr(y, p).statistic)
    assert pairwise_accuracy(y, y) == 1 and pairwise_accuracy(y, c) == 0.5
    assert pairwise_accuracy(c, y) is None
    assert np.isclose(top5_enrichment(y, y), 6 / 5) and np.isclose(top5_enrichment(y, c), 1.0)
    assert paired_bootstrap(y, c, y)["undefined_fraction"] == 1
    assert paired_bootstrap(y, c, y)["ci95"] == [None, None]
    assert paired_bootstrap(y, p, p)["ci95"] == [0, 0]
    assert paired_bootstrap(y, -y, y)["ci95"] == [2, 2]
    shams = {"shams": {"shuffle": {"n_at_or_above_real_delta": 1}, "gaussian": {"n_at_or_above_real_delta": 0}}}
    assert interpretation("primary", 0.5, True, shams) == "Sham reproduced the gain; not attributable to cooperativity"
    assert interpretation("primary", 0, False, shams) == "Null/sub-threshold; uninformative"
    result = evaluate(prepared, training, labels, suffix)
    stored = read_json(RESULTS / f"metrics{suffix}.json")
    result["seed_order"] = list(result["seed_order"])
    assert stored == result
    assert (RESULTS / f"results{suffix}.md").read_text() == results_markdown(result)
    for label in labels:
        for variant in VARIANTS:
            for arm in ARMS:
                for fold in result["labels"][label]["variants"][variant]["arms"][arm]["per_fold"]:
                    if fold["n"] == 1:
                        assert fold["spearman"] is None and fold["pairwise_accuracy"] is None
    print("PASS: average ties, constant vectors, singleton folds, paired bootstrap and saved metric reproducibility")


if __name__ == "__main__":
    labels, suffix = run_scope()
    prepared = read_prepared()
    training = read_json(RUN / f"training{suffix}.json")
    configure_cpu(training["threads"])
    check_sources(prepared, training)
    workspace = check_model(prepared, training)
    replayed = check_training(prepared, training, workspace, labels)
    check_metrics(prepared, training, labels, suffix)
    write_json(RESULTS / f"checks{suffix}.json", {
        "status": "passed", "prespec_sha256": prepared["prespec_sha256"],
        "training_fits": len(training["runs"]), "replayed_fits": replayed,
        "metric_results_sha256": sha256(RESULTS / f"metrics{suffix}.json"),
    })
    print("All protocol checks passed.")
