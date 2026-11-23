# Aurora A ChemBERTa A/B

This is a held-out functional comparison, separate from the parent repository’s structure benchmark. The independent primary endpoint is pKd in the human Aurora A Kinobeads binding assay [CHEMBL3991630](https://www.ebi.ac.uk/chembl/assay_report_card/CHEMBL3991630/) (15 compounds, 14 scaffold folds). No comparable independent numeric endpoint was found for all 24 compounds. The Lake no-TPX2 pKd sensitivity covers 24 compounds and 23 folds, but is circular: its label shares the global fit for the cooperativity feature. Do not treat it as independent evidence. A null on the small primary panel is uninformative.

Model A is ligand-only, matching PPARγ. Model B adds only training-fold-standardised log10(cooperativity); there were no eligible geometry or HDX columns. The LoRA comparison is primary; frozen ChemBERTa supports it. Controls are the training-fold mean, cooperativity alone, ten row shuffles and ten Gaussian shams. Reproduction by a sham rules out attributing a gain to cooperativity. See the sealed [prespec](PRESPEC.md) and [source notes](data/README.md).

## Rerun

From this directory on Linux, use Python 3.10 and a CPU PyTorch installation. The original run used eight workers with one Torch thread each; the LoRA fits take hours. The commands below write restartable fit checkpoints to `run/`. Do not run two training invocations against the same `run/` concurrently.

```sh
python3.10 -m venv .venv
. .venv/bin/activate
pip install torch==2.9.1 --index-url https://download.pytorch.org/whl/cpu
pip install numpy==2.2.6 scipy==1.15.3 rdkit==2026.9.2 transformers==4.57.6 peft==0.18.1 \
  accelerate==1.15.0 huggingface-hub==0.36.2 tokenizers==0.22.2
mkdir -p build run results
cp CHECKPOINT.sha256 build/seyonec-761d6a18.sha256
hf download seyonec/ChemBERTa-zinc-base-v1 \
  --revision 761d6a18cf99db371e0b43baf3e2d21b3e865a20 --cache-dir .hf/hub
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
export AURORA_TORCH_THREADS=1 AURORA_SEED_WORKERS=8 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python scripts/prepare.py
for endpoint in primary coupled; do
  python scripts/train.py --label "$endpoint" || exit 1
  python scripts/evaluate.py --label "$endpoint" || exit 1
  python scripts/check.py --label "$endpoint" || exit 1
done
python scripts/train.py       # Assembles existing fits after both endpoint runs.
python scripts/evaluate.py   # Writes results/results.md and results/metrics.json.
python scripts/check.py
```

The production run used endpoint-first dispatch at the requester’s direction; the sealed prespec remains unedited. The tracked tables suffice to rerun the experiment. The published run’s artifact hashes refer to its original 39-file audit; the public subset has its own manifest, so the source-manifest and dependent artifact hashes differ on a rerun. No omitted audit download is a model input.
