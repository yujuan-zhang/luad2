"""Mock OncoKB-style targeted-drug matching.

MOCK: real integration will call the OncoKB API (requires a license
token). The interface (`match_drugs`) stays the same.
"""

_knowledge_base = {
    ("EGFR", "p.L858R"): [
        {"drug": "Osimertinib", "level": "1"},
        {"drug": "Erlotinib", "level": "1"},
    ],
    ("KRAS", "p.G12C"): [{"drug": "Sotorasib", "level": "1"}],
    ("BRAF", "p.V600E"): [{"drug": "Dabrafenib + Trametinib", "level": "1"}],
}


def match_drugs(variants):
    matches = []
    for v in variants:
        key = (v["gene"], v["protein_change"])
        for hit in _knowledge_base.get(key, []):
            matches.append({"gene": v["gene"], "protein_change": v["protein_change"], **hit})
    return matches
