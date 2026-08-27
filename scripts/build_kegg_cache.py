"""One-time (re-runnable) build of the local KEGG pathway cache.

Downloads a PNG + gene-box-coordinate JSON for every KEGG pathway in five
of KEGG's own official BRITE categories (Signal transduction, Cancer:
overview, Cancer: specific types, Cell growth and death, Immune system) --
an objective, KEGG-defined scope rather than a hand-picked list, and
comprehensive enough that most driver/passenger genes in a real uploaded
VCF will land in a rendered pathway instead of falling outside an
arbitrarily small set.

Everything is written to pipelines/downstream/kegg_cache/, committed to
git, and read entirely offline at request time (see pathway.py) -- no
network call happens when the app itself runs, in prod or locally. Re-run
this script only when you want to refresh/expand the cache; it's not part
of the request path.

Usage:
    python -m scripts.build_kegg_cache
"""
import json
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import gseapy as gp
import requests

CACHE_DIR = Path(__file__).resolve().parents[1] / "pipelines" / "downstream" / "kegg_cache" / "pathways"
INDEX_PATH = CACHE_DIR.parent / "pathway_index.json"

TARGET_CATEGORIES = {
    "09132 Signal transduction",
    "09161 Cancer: overview",
    "09162 Cancer: specific types",
    "09143 Cell growth and death",
    "09151 Immune system",
}

# KEGG's official pathway name doesn't always match gseapy's KEGG_2021_Human
# key verbatim -- KEGG has renamed/split a few pathways since that library
# snapshot (2021). Manual overrides for the ones with an unambiguous older
# equivalent; anything else with no match is skipped rather than guessed.
GSEAPY_KEY_OVERRIDES = {
    "hsa04392": "Hippo signaling pathway",  # KEGG: "... - multiple species"
    "hsa04215": "Apoptosis",                # KEGG: "... - multiple species"
}

REQUEST_DELAY_SECONDS = 0.35  # keep well under KEGG's informal ~3 req/s guidance


def select_pathways():
    """{pathway_id: kegg_name} for every pathway in TARGET_CATEGORIES."""
    text = requests.get("https://rest.kegg.jp/get/br:hsa00001", timeout=30).text
    cur_b = None
    selected = {}
    for line in text.splitlines():
        if line.startswith("B"):
            cur_b = line[2:].strip() if len(line) > 2 else None
        elif line.startswith("C"):
            m = re.search(r"\d{5} (.+?)\s*\[PATH:(hsa\d+)\]", line)
            if m and cur_b in TARGET_CATEGORIES:
                name, pid = m.groups()
                selected[pid] = name.strip()
    return selected


def resolve_gseapy_keys(pathways):
    """(pid -> (kegg_name, gseapy_key)) for pathways with a usable gene set;
    prints and drops anything with no match."""
    kegg_sets = gp.get_library("KEGG_2021_Human")
    lower_to_key = {k.lower(): k for k in kegg_sets.keys()}
    resolved, dropped = {}, []
    for pid, name in pathways.items():
        gsea_key = GSEAPY_KEY_OVERRIDES.get(pid) or lower_to_key.get(name.lower())
        if gsea_key:
            resolved[pid] = (name, gsea_key)
        else:
            dropped.append((pid, name))
    if dropped:
        print(f"Skipping {len(dropped)} pathway(s) with no gseapy KEGG_2021_Human gene set:")
        for pid, name in dropped:
            print(f"  {pid}  {name}")
    return resolved


def build_id_to_symbol():
    """{'hsa:1956': 'EGFR', ...} for every human gene KEGG knows about."""
    text = requests.get("https://rest.kegg.jp/list/hsa", timeout=60).text
    mapping = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        kid, names_desc = parts[0], parts[3]
        primary = names_desc.split(";")[0].split(",")[0].strip()
        if primary:
            mapping[kid] = primary
    return mapping


def fetch_nodes(pathway_id, id_to_symbol):
    kgml = requests.get(f"https://rest.kegg.jp/get/{pathway_id}/kgml", timeout=30).text
    root = ET.fromstring(kgml)
    nodes = {}
    for entry in root.findall("entry"):
        if entry.get("type") != "gene":
            continue
        graphics = entry.find("graphics")
        if graphics is None:
            continue
        x, y = int(float(graphics.get("x"))), int(float(graphics.get("y")))
        w, h = int(float(graphics.get("width"))), int(float(graphics.get("height")))
        for kid in entry.get("name", "").split():
            symbol = id_to_symbol.get(kid)
            if symbol:
                nodes[symbol.upper()] = [x, y, w, h]
    return nodes


def main():
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    print("Selecting pathways from KEGG's BRITE hierarchy...")
    pathways = select_pathways()
    print(f"  {len(pathways)} pathways in target categories")

    resolved = resolve_gseapy_keys(pathways)
    print(f"  {len(resolved)} have a usable gseapy gene set")

    print("Building KEGG gene-ID -> symbol map (~24k genes, one request)...")
    id_to_symbol = build_id_to_symbol()

    index = {}
    failures = []
    for i, (pid, (kegg_name, gsea_key)) in enumerate(sorted(resolved.items()), start=1):
        print(f"[{i}/{len(resolved)}] {pid}  {kegg_name}")
        try:
            nodes = fetch_nodes(pid, id_to_symbol)
            time.sleep(REQUEST_DELAY_SECONDS)
            png_bytes = requests.get(f"https://rest.kegg.jp/get/{pid}/image", timeout=30).content
            time.sleep(REQUEST_DELAY_SECONDS)
            if not nodes or not png_bytes:
                raise ValueError("empty nodes or image")
        except Exception as exc:  # noqa: BLE001 -- best-effort batch build, log and continue
            print(f"    FAILED: {exc}")
            failures.append((pid, kegg_name, str(exc)))
            continue

        (CACHE_DIR / f"{pid}.png").write_bytes(png_bytes)
        (CACHE_DIR / f"{pid}_nodes.json").write_text(json.dumps(nodes))
        index[kegg_name] = [pid, gsea_key]

    INDEX_PATH.write_text(json.dumps(index, indent=2, sort_keys=True))
    print(f"\nWrote {len(index)} pathways to {CACHE_DIR}")
    print(f"Wrote index to {INDEX_PATH}")
    if failures:
        print(f"\n{len(failures)} failure(s):")
        for pid, name, err in failures:
            print(f"  {pid}  {name}  -- {err}")


if __name__ == "__main__":
    main()
