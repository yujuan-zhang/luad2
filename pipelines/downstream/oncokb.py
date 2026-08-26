"""Targeted-drug matching without an OncoKB token.

Two sources, combined:
  1. `DRUG_KB` — curated FDA-approved gene/mutation -> drug knowledge base
     (NCCN NSCLC guideline pattern). Always available offline.
  2. CIViC (civicdb.org) GraphQL API — free, no token. Queried live by
     molecular profile name (e.g. "EGFR L858R"). Best-effort: network
     failures return [] rather than breaking the pipeline.

`is_actionable()` decides what counts as a real targeted-therapy match for
routing purposes: KB hits (curated as FDA-approved only) and CIViC evidence
of level A/B with SENSITIVITYRESPONSE significance. Weaker CIViC evidence
(level C/D, case studies) is still returned for visibility but does not
count as actionable -- those mutations still flow to neoantigen ranking.

A real OncoKB integration (once a license token is available) can be added
as a third source later without changing `match_drugs`'s interface.
"""
import re

import requests

CIVIC_GRAPHQL_URL = "https://civicdb.org/api/graphql"

DRUG_KB = {
    "EGFR": [
        {"drug": "Osimertinib", "drug_class": "3rd-gen EGFR TKI", "line": "1L",
         "pattern": r"L858R|exon19|G719|S768|L861"},
        {"drug": "Erlotinib", "drug_class": "1st-gen EGFR TKI", "line": "1L",
         "pattern": r"L858R|exon19|G719|S768|L861"},
        {"drug": "Afatinib", "drug_class": "2nd-gen EGFR TKI", "line": "1L",
         "pattern": r"L858R|exon19|G719|S768|L861"},
        {"drug": "Amivantamab", "drug_class": "EGFR+MET bispecific Ab", "line": "1L",
         "pattern": r"exon20ins"},
    ],
    "KRAS": [
        {"drug": "Sotorasib", "drug_class": "KRAS G12C inhibitor", "line": "2L", "pattern": r"G12C"},
        {"drug": "Adagrasib", "drug_class": "KRAS G12C inhibitor", "line": "2L", "pattern": r"G12C"},
    ],
    "BRAF": [
        {"drug": "Dabrafenib + Trametinib", "drug_class": "BRAF+MEK inhibitor", "line": "1L",
         "pattern": r"V600E"},
    ],
    "MET": [
        {"drug": "Tepotinib", "drug_class": "MET inhibitor", "line": "1L", "pattern": r"exon14|splice"},
        {"drug": "Capmatinib", "drug_class": "MET inhibitor", "line": "1L", "pattern": r"exon14|splice"},
    ],
    "ALK": [
        {"drug": "Alectinib", "drug_class": "2nd-gen ALK TKI", "line": "1L", "pattern": r"fusion|.*"},
        {"drug": "Lorlatinib", "drug_class": "3rd-gen ALK TKI", "line": "1L", "pattern": r"fusion|.*"},
    ],
    "ROS1": [
        {"drug": "Entrectinib", "drug_class": "ROS1/NTRK inhibitor", "line": "1L", "pattern": r"fusion|.*"},
        {"drug": "Crizotinib", "drug_class": "ALK/ROS1/MET inhibitor", "line": "1L", "pattern": r"fusion|.*"},
    ],
    "RET": [
        {"drug": "Selpercatinib", "drug_class": "RET-selective inhibitor", "line": "1L", "pattern": r"fusion|.*"},
        {"drug": "Pralsetinib", "drug_class": "RET-selective inhibitor", "line": "1L", "pattern": r"fusion|.*"},
    ],
    "ERBB2": [
        {"drug": "Trastuzumab deruxtecan", "drug_class": "HER2 ADC", "line": "2L", "pattern": r".*"},
    ],
}

_CIVIC_QUERY = """
query($profileName: String!) {
  evidenceItems(molecularProfileName: $profileName, evidenceType: PREDICTIVE, first: 15) {
    edges {
      node {
        evidenceLevel
        significance
        disease { name }
        therapies { name }
      }
    }
  }
}
"""


def _kb_matches(gene, protein_change, consequence):
    matches = []
    target = f"{protein_change} {consequence}"
    for d in DRUG_KB.get(gene, []):
        if re.search(d["pattern"], target, re.IGNORECASE):
            matches.append({
                "source": "internal_kb",
                "gene": gene,
                "protein_change": protein_change,
                "drug": d["drug"],
                "drug_class": d["drug_class"],
                "line": d["line"],
                "level": "FDA-Approved",
                "significance": "SENSITIVITYRESPONSE",
            })
    return matches


def _civic_matches(gene, protein_change, timeout=8):
    short = re.sub(r"^p\.", "", str(protein_change))
    if not short or short == "NA":
        return []
    profile_name = f"{gene} {short}"
    try:
        resp = requests.post(
            CIVIC_GRAPHQL_URL,
            json={"query": _CIVIC_QUERY, "variables": {"profileName": profile_name}},
            timeout=timeout,
        )
        resp.raise_for_status()
        edges = resp.json()["data"]["evidenceItems"]["edges"]
    except (requests.RequestException, KeyError, ValueError, TypeError):
        return []

    matches = []
    for edge in edges:
        node = edge["node"]
        for therapy in node.get("therapies", []):
            matches.append({
                "source": "civic",
                "gene": gene,
                "protein_change": protein_change,
                "drug": therapy["name"],
                "drug_class": "CIViC evidence",
                "line": "—",
                "level": node.get("evidenceLevel", ""),
                "significance": node.get("significance", ""),
                "disease": node.get("disease", {}).get("name", ""),
            })
    return matches


def is_actionable(match):
    """Whether a match is strong enough to route the mutation to the
    targeted-therapy branch instead of neoantigen ranking."""
    if match["source"] == "internal_kb":
        return True
    return match.get("level") in ("A", "B") and match.get("significance") == "SENSITIVITYRESPONSE"


def match_drugs(variants, use_civic=True):
    matches = []
    for v in variants:
        gene, protein_change, consequence = v["gene"], v["protein_change"], v["consequence"]
        matches.extend(_kb_matches(gene, protein_change, consequence))
        if use_civic:
            matches.extend(_civic_matches(gene, protein_change))
    return matches
