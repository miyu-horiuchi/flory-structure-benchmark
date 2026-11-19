# Flory protein–ligand structure benchmark

Interactive results for seven models across PPARγ, GCK, AURKA and ERK2.

**Takeaway:** Chai-1 and ESMFold2 were most consistent for ligand pose and chemistry here; Protenix was close behind. This does not establish accurate partner-chain placement, dynamics or which model will fine-tune best.

Website: https://miyu-horiuchi.github.io/flory-structure-benchmark/

## Contents

- `index.html`: compiled interactive report, including the benchmark data.
- `src/report.body.html`: editable JSX source.
- `src/benchmark.data.json`: the original report data, unchanged for publication.
- `src/preview.head.html`: static title and social-sharing metadata.
- `src/publish.py`: combines the compiled report with the hosting metadata.

The site is published from the root of `main` using GitHub Pages, with no build step. The report loads the public design-system runtime from `static.devin.ai`; internet access is required. To preview locally, run `python3 -m http.server 8000` in this directory.

Compile source edits using Devin's `html_artifact_compile`, with `src/report.body.html` as `body_path`, `src/benchmark.data.json` as `data_path`, and an output path outside the repository. Then run `python3 src/publish.py /path/to/compiled-report.html` to produce `index.html`. This preserves the compiled artifact while adding static metadata that WhatsApp and other non-JavaScript crawlers can read. Run publishing tests with `python3 -m unittest discover -s src -p 'test_*.py'`.

This retrospective structural/chemical screen used public complexes and is not a held-out functional or dynamics benchmark. No fine-tuning was performed. The report contains scoring provenance and the separate Aurora WT diagnostics.
