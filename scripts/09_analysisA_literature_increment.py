#!/usr/bin/env python3
"""Analysis A: literature-increment audit.

Purpose (answers Reviewer #1 point 1 / Reviewer #2 point 1):
quantify what a literature review can and cannot deliver for the two claims
the manuscript is built on, so the rebuttal rests on counts rather than assertion.

Outputs (written next to this script):
  analysisA_query_counts.tsv      record counts per query
  analysisA_tlr5_sentences.tsv    every abstract sentence mentioning TLR5

The sentence file is deliberately not auto-classified: direction statements are
read by a human/agent from the extracted sentences, so no claim about the
literature is made by keyword matching alone.
"""

from __future__ import annotations

import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HERE = os.path.dirname(os.path.abspath(__file__))

QUERIES = {
    "q1_tlr5_and_pa_tiab": '(TLR5[tiab] OR "toll-like receptor 5"[tiab]) AND '
                          '("Pseudomonas aeruginosa"[tiab] OR "P. aeruginosa"[tiab])',
    "q2_tlr5_pa_expression": '((TLR5[tiab] OR "toll-like receptor 5"[tiab]) AND '
                             '("Pseudomonas aeruginosa"[tiab])) AND '
                             '(expression[tiab] OR transcript*[tiab] OR downregulat*[tiab] '
                             'OR upregulat*[tiab] OR attenuat*[tiab])',
    "q3_flagellin_tlr5_pa": '(flagellin[tiab]) AND (TLR5[tiab]) AND '
                            '("Pseudomonas aeruginosa"[tiab])',
    "q4_pa_interactome_network": '("Pseudomonas aeruginosa"[tiab]) AND '
                                 '(interactome[tiab] OR "host-pathogen interaction*"[tiab]) '
                                 'AND (network*[tiab] OR hub[tiab] OR centrality[tiab])',
    "q5_pa_interactome_bias": '("Pseudomonas aeruginosa"[tiab]) AND (bias[tiab] OR '
                              'ascertainment[tiab]) AND ("protein-protein interaction*"[tiab] '
                              'OR interactome[tiab] OR "interaction network"[tiab])',
    "q6_pa_interactome_cross_species": '("Pseudomonas aeruginosa"[tiab]) AND '
                                       '(cross-species[tiab] OR "cross-pathogen"[tiab] OR '
                                       'convergen*[tiab]) AND (interactome[tiab] OR '
                                       '"protein-protein interaction*"[tiab])',
}


def _get(endpoint: str, params: dict) -> str:
    params = dict(params)
    params.setdefault("retmode", "json" if endpoint.endswith("fcgi") else "text")
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8", "replace")


def esearch_count(term: str) -> int:
    import json

    payload = json.loads(_get("esearch.fcgi", {"db": "pubmed", "term": term, "retmax": "0"}))
    return int(payload["esearchresult"]["count"])


def esearch_ids(term: str, retmax: int = 200) -> list[str]:
    import json

    payload = json.loads(
        _get("esearch.fcgi", {"db": "pubmed", "term": term, "retmax": str(retmax)})
    )
    return payload["esearchresult"]["idlist"]


def efetch_abstracts(pmids: list[str]) -> list[dict]:
    records = []
    for start in range(0, len(pmids), 50):
        chunk = pmids[start:start + 50]
        xml_text = _get(
            "efetch.fcgi",
            {"db": "pubmed", "id": ",".join(chunk), "retmode": "xml"},
        )
        root = ET.fromstring(xml_text)
        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//PMID") or ""
            title = "".join(article.find(".//ArticleTitle").itertext()) if article.find(".//ArticleTitle") is not None else ""
            pieces = []
            for node in article.findall(".//Abstract/AbstractText"):
                label = node.get("Label")
                text = "".join(node.itertext())
                pieces.append(("%s: " % label if label else "") + text)
            year = article.findtext(".//PubDate/Year") or ""
            journal = article.findtext(".//Journal/ISOAbbreviation") or ""
            records.append(
                {"pmid": pmid, "year": year, "journal": journal, "title": title,
                 "abstract": " ".join(pieces)}
            )
        time.sleep(0.4)
    return records


def split_sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]


def main() -> None:
    counts = {}
    for name, term in QUERIES.items():
        counts[name] = esearch_count(term)
        print("%-34s %s" % (name, counts[name]))
        time.sleep(0.4)

    with open(os.path.join(HERE, "analysisA_query_counts.tsv"), "w", encoding="utf-8") as handle:
        handle.write("query\tn_pubmed_records\n")
        for name, count in counts.items():
            handle.write("%s\t%d\n" % (name, count))

    pmids = esearch_ids(QUERIES["q1_tlr5_and_pa_tiab"], retmax=200)
    print("fetching %d abstracts for the TLR5 x PA query" % len(pmids))
    records = efetch_abstracts(pmids)

    hits = []
    for record in records:
        for sentence in split_sentences(record["abstract"]):
            if re.search(r"\bTLR5\b|toll-like receptor 5", sentence, re.IGNORECASE):
                hits.append((record["pmid"], record["year"], record["journal"],
                             record["title"], sentence))

    with open(os.path.join(HERE, "analysisA_tlr5_sentences.tsv"), "w", encoding="utf-8") as handle:
        handle.write("pmid\tyear\tjournal\ttitle\tsentence\n")
        for row in hits:
            handle.write("\t".join(field.replace("\t", " ").replace("\n", " ") for field in row) + "\n")

    print("abstracts retrieved: %d" % len(records))
    print("sentences mentioning TLR5: %d" % len(hits))
    print("wrote analysisA_query_counts.tsv and analysisA_tlr5_sentences.tsv")


if __name__ == "__main__":
    main()
