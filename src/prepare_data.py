"""Combine the original benchmark and separately labelled personal AF3 results."""

import argparse
import json
from pathlib import Path


def combine_data(original: dict, personal: dict) -> dict:
    if "personal_af3" in original:
        raise ValueError("Original benchmark already contains personal AF3 data")
    if not personal.get("notice", "").endswith("https://alphafoldserver.com/output-terms"):
        raise ValueError("Personal AF3 results must retain the output terms notice")
    targets = set(original["target_names"])
    cells = personal["cells"]
    samples = personal["samples"]
    if len(cells) != len(targets) or {c["target"] for c in cells} != targets:
        raise ValueError("Personal AF3 must contain one cell for each original target")
    if len(samples) != 5 * len(cells):
        raise ValueError("Personal AF3 must retain all five samples per target")
    for cell in cells:
        rows = [r for r in samples if r["case"] == cell["case"]]
        if len(rows) != 5 or {r["sample"] for r in rows} != set(range(5)):
            raise ValueError("Personal AF3 samples must be unique and complete")
        if any(r["target"] != cell["target"] for r in rows):
            raise ValueError("Personal AF3 target and case must agree")
        if sum(r["evaluation_status"] == "scored" for r in rows) != cell["scored_samples"]:
            raise ValueError("Personal AF3 scorable count disagrees with samples")
        if sum(bool(r.get("clean_screen_hit")) for r in rows) != cell["clean_hits"]:
            raise ValueError("Personal AF3 pass count disagrees with samples")
        selected = max(rows, key=lambda r: (r["ranking_value"], -r["sample"]))
        if selected["sample"] != cell["selected_sample"]:
            raise ValueError("Personal AF3 selection must follow native confidence")
    return {**original, "personal_af3": personal}


def main() -> None:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    original = json.loads((root / "benchmark.data.json").read_text())
    personal = json.loads((root.parent / "personal-af3/results.json").read_text())
    args.output.write_text(json.dumps(combine_data(original, personal)))


if __name__ == "__main__":
    main()
