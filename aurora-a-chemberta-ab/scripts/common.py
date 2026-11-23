import argparse
import hashlib
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA, RUN, RESULTS = (ROOT / name for name in ("data", "run", "results"))
GROUPED_LABEL = {"1": "primary", "coupled": "coupled"}.get(os.environ.get("AURORA_GROUPED"))
GROUPED_OUTPUTS = {"primary": "grouped", "coupled": "grouped-coupled"}
GROUPED_PLANS = {"primary": ROOT / "EXPLORATORY.md", "coupled": ROOT / "EXPLORATORY-LAKE.md"}
if GROUPED_LABEL:
    RUN, RESULTS = RUN / GROUPED_OUTPUTS[GROUPED_LABEL], RESULTS / GROUPED_OUTPUTS[GROUPED_LABEL]
REVISION = "761d6a18cf99db371e0b43baf3e2d21b3e865a20"
SNAPSHOT = ROOT / ".hf/hub/models--seyonec--ChemBERTa-zinc-base-v1/snapshots" / REVISION
MODEL_MANIFEST = ROOT / "build/seyonec-761d6a18.sha256"
LABELS = ("primary", "coupled")
VARIANTS = ("lora", "frozen")
ARMS = ("A", "B")
SHAMS = {block: tuple(f"{block}_{i}" for i in range(1, 11)) for block in ("shuffle", "gaussian")}
FIT_ARMS = (*ARMS, "features_only", *(arm for draws in SHAMS.values() for arm in draws))
SEEDS = tuple(range(20261123, 20261128))
BOOTSTRAP_N = 1500
BATCH = 8
EPOCHS = 20


def run_scope():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", choices=LABELS, help="Run only the selected endpoint")
    label = parser.parse_args().label
    if GROUPED_LABEL:
        if label and label != GROUPED_LABEL:
            parser.error(f"grouped sensitivity is configured for {GROUPED_LABEL}")
        return (GROUPED_LABEL,), ""
    return ((label,), f"-{label}") if label else (LABELS, "")


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else _digest(handle)


def _digest(handle):
    digest = hashlib.sha256()
    for block in iter(lambda: handle.read(1 << 20), b""):
        digest.update(block)
    return digest.hexdigest()


def verify_manifest(manifest, directory):
    entries = {}
    for line in manifest.read_text().splitlines():
        digest, filename = line.split(maxsplit=1)
        path = directory / filename
        if sha256(path) != digest:
            raise ValueError(f"SHA256 mismatch: {path}")
        entries[filename] = digest
    if not entries:
        raise ValueError(f"Empty manifest: {manifest}")
    return entries


def prespec_sha():
    return verify_manifest(ROOT / "PRESPEC.sha256", ROOT)["PRESPEC.md"]


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def read_prepared():
    prepared = read_json(RUN / "prepared.json")
    if prepared["prespec_sha256"] != prespec_sha():
        raise ValueError("Prepared inputs do not match the sealed prespec")
    if prepared["source_manifest_sha256"] != sha256(DATA / "SHA256SUMS"):
        raise ValueError("Prepared inputs do not match the source manifest")
    if GROUPED_LABEL and prepared["exploratory"]["plan_sha256"] != sha256(GROUPED_PLANS[GROUPED_LABEL]):
        raise ValueError("Grouped inputs do not match the exploratory plan")
    return prepared
