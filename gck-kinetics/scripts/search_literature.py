"""Repeat the targeted literature and patent searches; retain result metadata."""

from concurrent.futures import ThreadPoolExecutor
import json
from urllib.parse import urlencode

import requests

from fetch_databases import DATA_DIR


EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
PUBMED = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PATENTS = "https://patents.google.com/xhr/query"
TERMS = '(SPR OR "surface plasmon resonance" OR "binding kinetics" OR "dissociation rate" OR "residence time" OR sensorgram OR "stopped-flow")'
NAMED = [
    "piragliatin OR RO4389620", "dorzagliatin OR HMS5552", "AZD1656 OR AZD6370",
    "MK-0941", "PF-04937319 OR PF-04991532", "GKA50", "LY2608204", "TTP399",
    "AMG-151 OR ARRY-403", "GKM-001", "ZYGK1", '"glucokinase activator"',
]
LIT_QUERIES = [
    '(glucokinase OR hexokinase-4) AND ("surface plasmon resonance" OR SPR OR sensorgram)',
    'glucokinase AND ("residence time" OR "binding kinetics" OR "dissociation rate" OR koff OR kon OR "stopped flow" OR "stopped-flow") AND (activator OR GKA)',
    f"TITLE_ABS:glucokinase AND (TITLE_ABS:activator OR TITLE_ABS:allosteric) AND {TERMS}",
    '"Optimizing glucokinase activator binding kinetics"',
    "DOI:10.1039/C4MD00027G", "DOI:10.1021/bi900374c", "DOI:10.1021/acs.jmedchem.6b00632",
] + [f"glucokinase AND {TERMS} AND ({name})" for name in NAMED]
PUBMED_QUERY = (
    'glucokinase[tiab] AND (activator*[tiab] OR allosteric[tiab]) AND '
    '("surface plasmon resonance"[tiab] OR SPR[tiab] OR "residence time"[tiab] OR '
    '"binding kinetics"[tiab] OR "dissociation rate"[tiab] OR koff[tiab] OR '
    'kon[tiab] OR stopped-flow[tiab])'
)
PATENT_QUERIES = [
    '"glucokinase activator"', '"glucokinase activator" "residence time"',
    '"glucokinase activator" "surface plasmon resonance"', '"glucokinase activator" "dissociation rate"',
    '"glucokinase activator" koff', '"glucokinase" "TAFMT"',
]


def fetch(job):
    service, url, params = job
    response = requests.get(url, params=params, timeout=60)
    try:
        content = response.json()
    except requests.exceptions.JSONDecodeError:
        content = response.text[:1000]
    return {"service": service, "url": response.url, "http_status": response.status_code,
            "response": content}


def main():
    jobs = [("Europe PMC", EPMC + "/search", {"query": q, "pageSize": 100, "format": "json"})
            for q in LIT_QUERIES]
    jobs += [("PubMed", PUBMED, {"db": "pubmed", "term": PUBMED_QUERY, "retmode": "json", "retmax": 100})]
    jobs += [("Google Patents", PATENTS, {"url": urlencode({"q": q}), "exp": ""})
             for q in PATENT_QUERIES]
    with ThreadPoolExecutor(max_workers=6) as pool:
        responses = list(pool.map(fetch, jobs))
    (DATA_DIR / "search_results.json").write_text(
        json.dumps(responses, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for r in responses:
        payload = r["response"]
        count = ""
        if isinstance(payload, dict):
            count = payload.get("hitCount", payload.get("esearchresult", {}).get(
                "count", payload.get("results", {}).get("total_num_results", "")))
        print(r["http_status"], count, r["url"])


if __name__ == "__main__":
    main()
