"""KEGG pathway viewer: overlay mutation/expression status on cached KEGG diagrams.

Not a call to pathview/cytoscape -- this renders the pathway image itself,
from a locally-cached KEGG PNG + pre-parsed gene-box coordinates
(`kegg_cache/pathways/*_nodes.json`, extracted from KEGG's KGML via the
KEGG REST API and not re-fetched here). Gene-to-pathway membership comes
from gseapy's KEGG_2021_Human gene sets, also resolved locally -- no
network calls at request time, in prod or locally.

`LUAD_PATHWAYS` covers 79 pathways from five of KEGG's own BRITE
categories (Signal transduction, Cancer: overview, Cancer: specific types,
Cell growth and death, Immune system) -- an objective, KEGG-defined scope
rather than a hand-picked list, so most driver/passenger genes in a real
uploaded VCF land in a rendered pathway rather than falling outside an
arbitrarily small set. Built (and re-buildable) by
`scripts/build_kegg_cache.py`; the mapping itself is data, not code --
loaded from `kegg_cache/pathway_index.json`.

Overlay colors (semi-transparent boxes; a gene with more than one status
gets its box split into vertical strips, one per status):
  green  = mutated gene (from vep.py output)
  red    = highly expressed (TPM >= HIGH_EXPR_TPM)
  yellow = expressed, not high (TPM >= EXPRESSED_TPM)

Unlike luad_workflow's version, there's no differential-expression fold
change here -- just raw tumor TPM against the two thresholds above.
"""
import base64
import io
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

import gseapy as gp
from PIL import Image, ImageDraw

CACHE_DIR = Path(__file__).parent / "kegg_cache" / "pathways"
_INDEX_PATH = Path(__file__).parent / "kegg_cache" / "pathway_index.json"

EXPRESSED_TPM = 1.0
HIGH_EXPR_TPM = 5.0

# {display_name: (kegg_pathway_id, gseapy_KEGG_2021_Human_key)}
LUAD_PATHWAYS = {
    name: tuple(v) for name, v in json.loads(_INDEX_PATH.read_text()).items()
}

NODE_COLORS = {
    "mut": (44, 160, 44, 200),  # green
    "high": (214, 39, 40, 190),  # red
    "expr": (255, 221, 70, 170),  # yellow
}
COLOR_ORDER = ["mut", "high", "expr"]


def load_all_pathway_genes():
    """{pathway_id: [symbols]} for the 8 LUAD pathways, via gseapy (offline, cached)."""
    kegg_sets = gp.get_library("KEGG_2021_Human")
    return {
        pid: [g.upper() for g in kegg_sets.get(gsea_key, [])]
        for _, (pid, gsea_key) in LUAD_PATHWAYS.items()
    }


def mutated_genes_from_variants(variants):
    return sorted({v["gene"] for v in variants if v["gene"] != "NA"})


def expression_tiers(expression_df, expressed_tpm=EXPRESSED_TPM, high_tpm=HIGH_EXPR_TPM):
    """Gene symbols clearing each TPM threshold. `high` implies `expressed`
    -- callers treat them as an ordinal ladder, not independent flags."""
    expr_lookup = dict(zip(expression_df["SYMBOL"].astype(str).str.upper(), expression_df["TPM_GENE"]))
    expressed = {g for g, tpm in expr_lookup.items() if tpm >= expressed_tpm}
    high = {g for g, tpm in expr_lookup.items() if tpm >= high_tpm}
    return expressed, high


def get_hit_pathways(mutated_genes, pathway_genes):
    """[(name, pathway_id, hit_genes)] for pathways with >=1 mutated gene,
    sorted by hit count descending."""
    mut_set = {g.upper() for g in mutated_genes}
    hits = []
    for name, (pid, _) in LUAD_PATHWAYS.items():
        hit = sorted(mut_set & set(pathway_genes.get(pid, [])))
        if hit:
            hits.append((name, pid, hit))
    hits.sort(key=lambda x: len(x[2]), reverse=True)
    return hits


def _node_positions(pathway_id):
    path = CACHE_DIR / f"{pathway_id}_nodes.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def render_pathway(pathway_id, mutated_genes, expressed_genes, high_expr_genes):
    """Overlay status on the cached KEGG PNG. Returns PNG bytes, or b""
    if this pathway isn't cached."""
    png_path = CACHE_DIR / f"{pathway_id}.png"
    nodes = _node_positions(pathway_id)
    if not png_path.exists() or not nodes:
        return b""

    mut_set = {g.upper() for g in mutated_genes}
    high_set = {g.upper() for g in high_expr_genes}
    expr_set = {g.upper() for g in expressed_genes} - high_set

    box_layers = defaultdict(set)
    for symbol, coords in nodes.items():
        key = tuple(coords)
        if symbol in mut_set:
            box_layers[key].add("mut")
        if symbol in high_set:
            box_layers[key].add("high")
        elif symbol in expr_set:
            box_layers[key].add("expr")

    img = Image.open(png_path).convert("RGBA")
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    for (cx, cy, w, h), keys in box_layers.items():
        layers = [NODE_COLORS[k] for k in COLOR_ORDER if k in keys]
        if not layers:
            continue
        x0, y0 = cx - w // 2, cy - h // 2
        x1, y1 = cx + w // 2, cy + h // 2
        strip_w = (x1 - x0) / len(layers)
        for i, color in enumerate(layers):
            sx0 = int(x0 + i * strip_w)
            sx1 = int(x0 + (i + 1) * strip_w)
            draw.rectangle([sx0, y0, sx1, y1], fill=color)

    result = Image.alpha_composite(img, overlay).convert("RGB")
    buf = io.BytesIO()
    result.save(buf, format="PNG")
    return buf.getvalue()


def build_kegg_url(pathway_id, mutated_genes, expressed_genes, high_expr_genes):
    """KEGG's own online colored-pathway link, kept as a fallback/reference view."""
    mut_set = {g.upper() for g in mutated_genes}
    high_set = {g.upper() for g in high_expr_genes} - mut_set
    expr_set = ({g.upper() for g in expressed_genes} - mut_set) - high_set

    lines = [f"{g}\tgreen" for g in sorted(mut_set)]
    lines += [f"{g}\tred" for g in sorted(high_set)]
    lines += [f"{g}\tyellow" for g in sorted(expr_set)]

    if not lines:
        return f"https://www.kegg.jp/pathway/{pathway_id}"
    multi_query = quote("\n".join(lines))
    return f"https://www.kegg.jp/kegg-bin/show_pathway?{pathway_id}&multi_query={multi_query}"


def build_gene_table(pathway_genes, variants, expressed_genes, high_expr_genes, tpm_lookup):
    """One row per mutated variant site in this pathway's genes -- a gene
    with two mutations (e.g. EGFR L62R and L858R) gets two rows -- with
    genomic position, DNA change, protein change, and the gene's own
    expression status attached (expression confirmation). Most highly
    expressed genes in any pathway are unmutated housekeeping/signaling
    genes -- listing every one of those would flood the table with noise,
    so only mutated genes appear here; expression on non-mutated genes is
    still shown in the image."""
    members = set(pathway_genes)
    high_set = {g.upper() for g in high_expr_genes}
    expr_set = {g.upper() for g in expressed_genes} - high_set

    rows = []
    for v in variants:
        gene = v["gene"].upper()
        if gene not in members:
            continue
        rows.append({
            "gene": v["gene"],
            "position": f"{v['chrom']}:{v['pos']}",
            "dna_change": f"{v['ref']}>{v['alt']}",
            "protein_change": v["protein_change"],
            "tumor_tpm": tpm_lookup.get(gene),
            "status": "High expression" if gene in high_set else ("Expressed" if gene in expr_set else "Not expressed"),
        })
    rows.sort(key=lambda r: (r["gene"], -(r["tumor_tpm"] or 0)))
    return rows


def analyze_pathways(variants, expression_df):
    """For each of LUAD_PATHWAYS' pathways with >=1 mutated gene, render an
    annotated PNG + gene status table. Returns [] if none were hit."""
    pathway_genes = load_all_pathway_genes()
    mutated = mutated_genes_from_variants(variants)
    expressed, high = expression_tiers(expression_df)
    tpm_lookup = dict(zip(expression_df["SYMBOL"].astype(str).str.upper(), expression_df["TPM_GENE"]))

    results = []
    for name, pid, hit_genes in get_hit_pathways(mutated, pathway_genes):
        # Restrict to this pathway's own members -- mutated/expressed/high
        # are otherwise genome-wide sets, and coloring every highly
        # expressed gene in the genome would make the KEGG URL absurdly long.
        members = set(pathway_genes[pid])
        pw_mutated = [g for g in mutated if g in members]
        pw_expressed = [g for g in expressed if g in members]
        pw_high = [g for g in high if g in members]

        png_bytes = render_pathway(pid, pw_mutated, pw_expressed, pw_high)
        results.append({
            "name": name,
            "pathway_id": pid,
            "mutated_hit_genes": hit_genes,
            "image_png_base64": base64.b64encode(png_bytes).decode() if png_bytes else None,
            "kegg_url": build_kegg_url(pid, pw_mutated, pw_expressed, pw_high),
            "gene_table": build_gene_table(pathway_genes[pid], variants, pw_expressed, pw_high, tpm_lookup),
        })
    return results
