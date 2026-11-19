# Flory protein–ligand structure benchmark

Interactive results for seven models across PPARγ, GCK, AURKA and ERK2.

**Takeaway:** Chai-1 and ESMFold2 were most consistent for ligand pose and chemistry here; Protenix was close behind. This does not establish accurate partner-chain placement, dynamics or which model will fine-tune best.

Website: https://miyu-horiuchi.github.io/flory-structure-benchmark/

## Contents

- `index.html`: compiled interactive report, including the benchmark data.
- `src/report.body.html`: editable JSX source.
- `src/benchmark.data.json`: the original report data, unchanged for publication.

The site is published from the root of `main` using GitHub Pages, with no build step. The report loads the public design-system runtime from `static.devin.ai`; internet access is required. To preview locally, run `python3 -m http.server 8000` in this directory.

Source edits are compiled using Devin's `html_artifact_compile`, with `src/report.body.html` as `body_path`, `src/benchmark.data.json` as `data_path`, and `index.html` as `output_path`.

This retrospective structural/chemical screen used public complexes and is not a held-out functional or dynamics benchmark. No fine-tuning was performed. The report contains scoring provenance and the separate Aurora WT diagnostics.
