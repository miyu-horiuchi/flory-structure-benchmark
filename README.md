# Flory protein–ligand structure benchmark

Interactive results for seven models across PPARγ, GCK, AURKA and ERK2, plus a separately labelled personal, non-commercial AlphaFold 3 Server comparison.

In the original seven-model screen, Chai-1 and ESMFold2 were most consistent for ligand pose and chemistry; Protenix was close behind. This does not establish accurate partner-chain placement, dynamics or which model will fine-tune best. AF3 is not incorporated into that controlled comparison: its GCK exports are incomplete and its Aurora ligand charge differs.

Website: https://miyu-horiuchi.github.io/flory-structure-benchmark/

## Contents

- `index.html`: compiled interactive report, including the benchmark data.
- `src/report.body.html`: editable JSX source.
- `src/benchmark.data.json`: the original report data, unchanged for publication.
- `personal-af3/results.json`: separate personal AF3 metrics, native confidence, hashes, manifests and output-terms notice. Aligned coordinate arrays are omitted from this display export; raw structures are not redistributed here.
- `personal-af3/output-terms.txt`: the preserved AlphaFold Server Output Terms.
- `src/prepare_data.py`: combines both datasets without changing any original benchmark fields.
- `src/preview.head.html`: static title and social-sharing metadata.
- `src/publish.py`: combines the compiled report with the hosting metadata.

The site is published from the root of `main` using GitHub Pages, with no build step. The report loads the public design-system runtime from `static.devin.ai`; internet access is required. To preview locally, run `python3 -m http.server 8000` in this directory.

First run `python3 src/prepare_data.py /path/to/combined-data.json`, with an output path outside the repository. Compile source edits using Devin's `html_artifact_compile`, with `src/report.body.html` as `body_path`, the combined JSON as `data_path`, and an output path outside the repository. Then run `python3 src/publish.py /path/to/compiled-report.html` to produce `index.html`. This preserves the compiled artifact while adding static metadata that WhatsApp and other non-JavaScript crawlers can read. Run publishing/data tests with `python3 -m unittest discover -s src -p 'test_*.py'`.

This retrospective structural/chemical screen used public complexes and is not a held-out functional or dynamics benchmark. No fine-tuning was performed. The report contains scoring provenance and the separate Aurora WT diagnostics.

## Personal AlphaFold 3 comparison

The AF3 section is personal, non-commercial exploratory work, separate from the original Flory benchmark and not for company research. GitHub Pages is publicly accessible; a personal GitHub account does not make the page private or grant commercial rights. See [the comparison notes](personal-af3/README.md) for the GCK failure, AURKA chemistry difference, partner-placement caveat and citation.

**This information is subject to AlphaFold Server Output Terms of Use found at https://alphafoldserver.com/output-terms**

This notice applies to the AF3 output and derivatives, including the displayed metrics, JSON and AF3-related page content. The [preserved full terms](personal-af3/output-terms.txt) apply separately from rights in the original seven-model data. Do not use AF3 output for company/commercial work, automated docking or interaction-screening systems, or training similar biomolecular structure-prediction models. No Google endorsement is implied.
