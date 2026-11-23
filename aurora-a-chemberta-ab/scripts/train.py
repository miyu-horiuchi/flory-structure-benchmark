import hashlib
import math
import multiprocessing as mp
import os
import random
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import torch
from peft import LoraConfig, TaskType, get_peft_model
from torch import nn
from transformers import AutoModel, AutoTokenizer

from common import (
    BATCH, EPOCHS, FIT_ARMS, LABELS, MODEL_MANIFEST, REVISION, RUN, SEEDS, SHAMS, SNAPSHOT,
    VARIANTS, read_json, read_prepared, run_scope, sha256, verify_manifest, write_json,
)


@dataclass(frozen=True)
class FoldData:
    train_ids: tuple
    test_ids: tuple
    train_pytorch_y: torch.Tensor
    train_pytorch_feature: torch.Tensor
    test_pytorch_feature: torch.Tensor
    scalers: dict


def population_scale(values):
    mean, std = float(np.mean(values)), float(np.std(values, ddof=0))
    return {"mean": mean, "sd": std or 1.0}


def scale(values, stats):
    return torch.as_tensor((np.asarray(values) - stats["mean"]) / stats["sd"], dtype=torch.float32).reshape(-1, 1)


def feature_blocks(compounds, labels, label):
    names = sorted(labels)
    values = np.array([compounds[name]["log10_alpha"] for name in names])
    real = dict(zip(names, values))
    features = {arm: real for arm in ("A", "B", "features_only")}
    for block_index, (block, arms) in enumerate(SHAMS.items()):
        for draw, arm in enumerate(arms, 1):
            seed = np.random.SeedSequence([SEEDS[0], LABELS.index(label), block_index, draw])
            rng = np.random.default_rng(seed)
            sample = rng.permutation(values) if block == "shuffle" else rng.standard_normal(values.shape)
            features[arm] = dict(zip(names, sample))
    return features


def partition(features, labels, test_ids, arm):
    test_ids = tuple(test_ids)
    train_ids = tuple(sorted(set(labels).difference(test_ids)))
    train_y = [labels[name] for name in train_ids]
    train_dyn = [features[name] for name in train_ids]
    test_dyn = [features[name] for name in test_ids]
    scalers = {"label": population_scale(train_y), "dynamics": population_scale(train_dyn)}
    train_feature = scale(train_dyn, scalers["dynamics"])
    test_feature = scale(test_dyn, scalers["dynamics"])
    if arm == "A":
        train_feature = torch.zeros_like(train_feature)
        test_feature = torch.zeros_like(test_feature)
    return FoldData(train_ids, test_ids, scale(train_y, scalers["label"]), train_feature, test_feature, scalers)


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def tensor_digest(named_tensors):
    h = hashlib.sha256()
    for name, tensor in named_tensors:
        h.update(name.encode() + b"\0")
        h.update(str(tuple(tensor.shape)).encode() + b"\0")
        h.update(tensor.detach().contiguous().numpy().tobytes())
    return h.hexdigest()


def configure_cpu(threads):
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)


def load_encoder(compounds, *, precompute=True):
    paths = verify_manifest(MODEL_MANIFEST, SNAPSHOT)
    encoder, loading = AutoModel.from_pretrained(
        SNAPSHOT, local_files_only=True, add_pooling_layer=False, output_loading_info=True,
        hidden_dropout_prob=0.0, attention_probs_dropout_prob=0.0,
    )
    assert not any(loading[key] for key in ("missing_keys", "mismatched_keys", "error_msgs")), loading
    assert all(key.startswith(("lm_head.", "pooler.")) for key in loading["unexpected_keys"]), loading
    assert encoder.config.hidden_size == 768
    assert encoder.config.hidden_dropout_prob == encoder.config.attention_probs_dropout_prob == 0.0
    assert all(module.p == 0.0 for module in encoder.modules() if isinstance(module, nn.Dropout))
    encoder.requires_grad_(False).eval()
    report = {
        "revision": REVISION, "file_sha256": paths, "loading": loading, "dtype": str(encoder.dtype),
        "hidden_size": encoder.config.hidden_size, "num_hidden_layers": encoder.config.num_hidden_layers,
        "base_sha256": tensor_digest(encoder.named_parameters()),
    }
    tokenizer = AutoTokenizer.from_pretrained(SNAPSHOT, local_files_only=True)
    ids = sorted(compounds)
    smiles = {name: compounds[name]["smiles"] for name in ids}
    workspace = Workspace(encoder, tokenizer, smiles, {})
    if precompute:
        with torch.no_grad():
            for start in range(0, len(ids), BATCH):
                names = ids[start:start + BATCH]
                workspace.frozen.update(zip(names, workspace.embed(encoder, names)))
    return workspace, report


@dataclass
class Workspace:
    base: nn.Module
    tokenizer: object
    smiles: dict
    frozen: dict

    def embed(self, encoder, ids):
        batch = self.tokenizer(
            [self.smiles[name] for name in ids], padding=True, truncation=True, max_length=256, return_tensors="pt",
        )
        hidden = encoder(**batch).last_hidden_state
        mask = batch["attention_mask"].unsqueeze(-1)
        return (hidden * mask).sum(dim=1) / mask.sum(dim=1)


def train_cell(workspace, fold, variant, seed, arm):
    seed_all(seed)
    base = workspace.base
    encoder = base
    if variant == "lora":
        encoder = get_peft_model(base, LoraConfig(
            task_type=TaskType.FEATURE_EXTRACTION, r=8, lora_alpha=16, lora_dropout=0,
            bias="none", target_modules=["query", "value"],
        ))
    encoder.eval()
    try:
        return fit(workspace, encoder, fold, variant, seed, arm)
    finally:
        if variant == "lora":
            assert encoder.unload() is base
            del base.peft_config


def fit(workspace, encoder, fold, variant, seed, arm):
    head = nn.Sequential(nn.Linear(769, 32), nn.GELU(), nn.Linear(32, 1))
    named_adapters = [(name, p) for name, p in encoder.named_parameters() if p.requires_grad]
    if variant == "lora":
        assert len(named_adapters) == 4 * workspace.base.config.num_hidden_layers
        assert all((".query." in name or ".value." in name) and ".lora_" in name for name, _ in named_adapters)
    else:
        assert not named_adapters
    trainable = named_adapters + [(f"head.{name}", p) for name, p in head.named_parameters()]
    params = [p for _, p in trainable]
    init = tensor_digest(trainable)
    optimizer = torch.optim.AdamW(params, lr=1e-4, betas=(0.9, 0.999), eps=1e-8, weight_decay=0.01)
    order = torch.Generator().manual_seed(seed)
    batch_sha = hashlib.sha256()
    n = len(fold.train_ids)
    steps = 0

    def predict(names, features):
        if arm == "features_only":
            embedding = features.new_zeros((len(names), 768))
        elif variant == "frozen" and workspace.frozen:
            embedding = torch.stack([workspace.frozen[name] for name in names])
        else:
            embedding = workspace.embed(encoder, names)
        return head(torch.cat((embedding, features), dim=1))

    for _ in range(EPOCHS):
        permutation = torch.randperm(n, generator=order)
        batch_sha.update(permutation.numpy().tobytes())
        for indices in permutation.split(BATCH):
            names = [fold.train_ids[i] for i in indices.tolist()]
            optimizer.zero_grad(set_to_none=True)
            prediction = predict(names, fold.train_pytorch_feature[indices])
            nn.functional.mse_loss(prediction, fold.train_pytorch_y[indices]).backward()
            nn.utils.clip_grad_norm_(params, 1.0)
            optimizer.step()
            steps += 1

    with torch.no_grad():
        test = predict(fold.test_ids, fold.test_pytorch_feature).flatten().double().numpy()
    label = fold.scalers["label"]
    pred = (test * label["sd"] + label["mean"]).tolist()
    assert all(math.isfinite(value) for value in pred)
    record = {
        "train_ids": fold.train_ids, "test_ids": fold.test_ids, "scalers": fold.scalers,
        "steps": steps, "init_sha256": init, "final_sha256": tensor_digest(trainable),
        "n_trainable_parameters": sum(p.numel() for p in params), "n_trainable_tensors": len(params),
        "batch_sha256": batch_sha.hexdigest(),
        "training_inputs_sha256": tensor_digest([
            ("targets", fold.train_pytorch_y), ("dynamics", fold.train_pytorch_feature),
        ]),
        "oof": dict(zip(fold.test_ids, pred)),
    }
    return record


_worker = None
_worker_model = None


def init_worker(compounds, threads):
    global _worker, _worker_model
    configure_cpu(threads)
    _worker, _worker_model = load_encoder(compounds)


def run_seed(fold, variant, seed, arm):
    previous = {name: (id(param), param._version) for name, param in _worker.base.named_parameters()}
    record = train_cell(_worker, fold, variant, seed, arm)
    current = {name: (id(param), param._version) for name, param in _worker.base.named_parameters()}
    assert current == previous, "Frozen base parameters changed"
    record["base_sha256"] = _worker_model["base_sha256"]
    return record, _worker_model


def run_training(prepared, threads, workers, labels=LABELS):
    compounds = prepared["compounds"]
    blocks = {label: feature_blocks(compounds, prepared["labels"][label], label) for label in labels}
    prepared_sha = sha256(RUN / "prepared.json")
    checkpoint_root = RUN / f"cells-{prepared_sha[:12]}"
    checkpoint_root.mkdir(exist_ok=True)
    records = []
    model = None
    reused = 0
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=mp.get_context("spawn"), initializer=init_worker, initargs=(compounds, threads),
    ) as pool:
        for fold in prepared["folds"]:
            for label in labels:
                test_ids = fold["members"][label]
                if not test_ids:
                    continue
                for variant in VARIANTS:
                    for start in range(0, len(FIT_ARMS), max(1, workers // 2)):
                        jobs = []
                        for arm in FIT_ARMS[start:start + max(1, workers // 2)]:
                            data = partition(blocks[label][arm], prepared["labels"][label], test_ids, arm)
                            path = checkpoint_root / f"{fold['id']}-{label}-{variant}-{arm}.json"
                            if path.exists():
                                payload = read_json(path)
                                assert payload["threads"] == threads
                                assert len(payload["fits"]) == len(SEEDS)
                                reused += len(SEEDS)
                                jobs.append((arm, path, payload, ()))
                            else:
                                jobs.append((arm, path, None, [pool.submit(run_seed, data, variant, seed, arm)
                                                               for seed in SEEDS]))
                        for arm, path, payload, futures in jobs:
                            if payload is None:
                                fits = [future.result() for future in futures]
                                loaded = fits[0][1]
                                assert all(meta == loaded for _, meta in fits)
                                payload = {"threads": threads, "model": loaded, "fits": [fit for fit, _ in fits]}
                                write_json(path, payload)
                            loaded = payload["model"]
                            if model is None:
                                model = loaded
                                print(f"Loaded {model['revision']}; missing={model['loading']['missing_keys']}; "
                                      f"seed workers={workers}; threads per worker={threads}", flush=True)
                            assert loaded == model
                            records.extend({"fold": fold["id"], "label": label, "variant": variant,
                                            "arm": arm, "seed": seed, **record}
                                           for seed, record in zip(SEEDS, payload["fits"]))
                            if variant == "lora":
                                print(f"Finished {fold['id']} {label}/{variant}/{arm}; completed={len(records)}",
                                      flush=True)
                print(f"Finished {fold['id']} {label}; total completed fits={len(records)}", flush=True)
    print(f"Reused checkpointed fits: {reused}", flush=True)
    return {
        "prespec_sha256": prepared["prespec_sha256"], "prepared_sha256": prepared_sha,
        "model": model, "threads": threads, "seed_workers": workers, "runs": records,
    }


if __name__ == "__main__":
    labels, suffix = run_scope()
    threads = int(os.environ.get("AURORA_TORCH_THREADS", "1"))
    workers = int(os.environ.get("AURORA_SEED_WORKERS", str(len(SEEDS))))
    training = run_training(read_prepared(), threads, workers, labels)
    write_json(RUN / f"training{suffix}.json", training)
    print(f"Training fits: {len(training['runs'])}")
