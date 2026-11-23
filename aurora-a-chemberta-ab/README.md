# Aurora A ChemBERTa A/B

This is a held-out functional comparison, separate from the parent repository’s structure benchmark. The sealed primary endpoint is independent pKd in the human Aurora A Kinobeads binding assay [CHEMBL3991630](https://www.ebi.ac.uk/chembl/assay_report_card/CHEMBL3991630/) (15 compounds, 14 scaffold folds). The decision rule did not resolve a Model B benefit at this size; that is not evidence against cooperativity or dynamics. No comparable independent numeric endpoint was found for all 24 compounds. The Lake no-TPX2 pKd sensitivity covers 24 compounds and 23 folds, but is circular: its label shares the global fit for the cooperativity feature.

Model A is ligand-only, matching PPARγ. Model B adds only training-fold-standardised log10(cooperativity); there were no eligible geometry or HDX columns. The LoRA comparison is primary; frozen ChemBERTa supports it. Controls are the training-fold mean, cooperativity alone, ten row shuffles and ten Gaussian shams. Reproduction by a sham rules out attributing a gain to cooperativity. See the sealed [prespec](PRESPEC.md) and [source notes](data/README.md).

On singleton folds the training mean decreases as the held-out label rises, so the primary fold-mean-only ρ of −0.9991 is expected arithmetic, not a model diagnostic. Model A’s correlation with that mean was 0.5170 (LoRA) and 0.5796 (frozen). At or above the real Δρ were 0/10 shuffled and 1/10 Gaussian draws for LoRA, and 8/10 shuffled and 9/10 Gaussian draws for frozen.

The protocol owner’s power analysis found that the unchanged rule fires 38% of the time at n=26 when the true Δρ is exactly 0.3; regrouping moved that by under one percentage point. CI-excludes-zero alone had 31% power at n=26 and required about 104 compounds for 80% power. Grouped k-fold is a comparability analysis, not a rescue, and does not replace the sealed headline.

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
for grouped in 1 coupled; do
  (
    export AURORA_GROUPED="$grouped"
    python scripts/prepare.py && python scripts/train.py &&
      python scripts/evaluate.py && python scripts/check.py
  ) || exit 1
done
python scripts/train.py       # Assembles existing fits after both endpoint runs.
python scripts/evaluate.py   # Writes results/results.md and results/metrics.json.
python scripts/check.py
```

The production run used endpoint-first dispatch at the requester’s direction; the sealed prespec remains unedited. The tracked tables suffice to rerun the experiment. The published run’s artifact hashes refer to its original 39-file audit; the public subset has its own manifest, so the source-manifest and dependent artifact hashes differ on a rerun. No omitted audit download is a model input.

The post-hoc grouped sensitivities use [three five-compound folds for the independent endpoint](EXPLORATORY.md) and [four six-compound folds for the circular Lake endpoint](EXPLORATORY-LAKE.md), with typed scaffold groups kept intact. Grouped checkpoints and reports are kept separately; the combined report includes each grouped result only after its checks pass.
