import base64
import html
import inspect
import json
from pathlib import Path

import pandas as pd
import requests
import streamlit as st

API_URL = "http://localhost:8000"
PRECOMPUTED_PATH = Path(__file__).resolve().parent.parent / "data" / "demo" / "precomputed_result.json"

# st.image()'s width kwarg was renamed use_column_width -> use_container_width
# across streamlit versions; local dev and Streamlit Cloud can end up on
# different ones (Cloud resolves frontend/requirements.txt independently),
# so pick whichever this installed version actually supports.
_IMAGE_WIDTH_KWARG = (
    "use_container_width"
    if "use_container_width" in inspect.signature(st.image).parameters
    else "use_column_width"
)

_LEVEL_ORDER = {"FDA-Approved": 0, "A": 1, "B": 2, "C": 3, "D": 4}

# Generic-name systemic therapies that sometimes show up in CIViC evidence
# alongside real targeted agents (e.g. an EGFR L858R evidence item citing
# "Carboplatin" as part of a combo regimen) -- these aren't targeted
# therapies and shouldn't be presented as if they were.
_NON_TARGETED_DRUGS = {
    "chemotherapy", "carboplatin", "cisplatin", "pemetrexed", "paclitaxel",
    "docetaxel", "gemcitabine", "vinorelbine", "etoposide", "doxorubicin",
    "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "ipilimumab",
}


def _is_targeted_drug(name):
    return name.strip().lower() not in _NON_TARGETED_DRUGS


def _ic50_badge(ic50):
    if ic50 < 50:
        return "luad-badge--strong", "Strong binder"
    if ic50 < 150:
        return "luad-badge--info", "Binder"
    return "luad-badge--muted", "Weak binder"


def _level_badge_class(level):
    if level == "FDA-Approved":
        return "luad-badge--fda"
    if level in ("A", "B"):
        return "luad-badge--info"
    return "luad-badge--muted"


def _response_label(significance):
    sig = (significance or "").upper()
    if "SENSITIVITY" in sig or "RESPONSE" in sig:
        return "Sensitive", "luad-badge--fda"
    if "RESISTANCE" in sig:
        return "Resistant", "luad-badge--danger"
    return (significance.title() if significance else "—"), "luad-badge--muted"


def _drug_table_html(matches):
    """Small evidence table rendered as plain HTML (with real badges) instead
    of st.dataframe -- a native Streamlit grid can't put a colored pill in a
    cell, and the raw SENSITIVITYRESPONSE/drug_class strings read as
    database internals rather than a clinical summary."""
    rows = []
    for m in matches:
        level = m.get("level") or "—"
        line = m.get("line") or "—"
        resp_label, resp_class = _response_label(m.get("significance"))
        rows.append(
            "<tr style='border-bottom:1px solid #F1F5F9;'>"
            f"<td style='padding:0.5rem 0.7rem;'>{html.escape(m['gene'])}</td>"
            f"<td style='padding:0.5rem 0.7rem;'>{html.escape(m['protein_change'])}</td>"
            f"<td style='padding:0.5rem 0.7rem; font-weight:600;'>{html.escape(m['drug'])}</td>"
            f"<td style='padding:0.5rem 0.7rem;'><span class='luad-badge {_level_badge_class(level)}'>{html.escape(level)}</span></td>"
            f"<td style='padding:0.5rem 0.7rem; opacity:0.75;'>{html.escape(line)}</td>"
            f"<td style='padding:0.5rem 0.7rem;'><span class='luad-badge {resp_class}'>{resp_label}</span></td>"
            "</tr>"
        )
    header = (
        "<tr style='border-bottom:1px solid #E5E7EB; text-align:left; opacity:0.55; font-size:0.78rem;'>"
        "<th style='padding:0.4rem 0.7rem;'>Gene</th><th style='padding:0.4rem 0.7rem;'>Mutation</th>"
        "<th style='padding:0.4rem 0.7rem;'>Drug</th><th style='padding:0.4rem 0.7rem;'>Level</th>"
        "<th style='padding:0.4rem 0.7rem;'>Line</th><th style='padding:0.4rem 0.7rem;'>Response</th></tr>"
    )
    return (
        "<div class='luad-card' style='padding:0.4rem 0.2rem;'>"
        "<table style='width:100%; border-collapse:collapse; font-size:0.85rem;'>"
        f"<thead>{header}</thead><tbody>{''.join(rows)}</tbody></table>"
        "</div>"
    )


st.set_page_config(page_title="LUADtx", layout="wide")

st.markdown(
    """
    <style>
    .stApp { background-color: #F8FAFC; }
    /* Room for the fixed navbar/footer bars (see below) so page content
       and the sidebar's own content never sit underneath them. */
    .main .block-container { padding-top: 4.5rem; padding-bottom: 5rem; }
    [data-testid="stSidebar"] {
        min-width: 230px; max-width: 230px;
        background-color: #FFFFFF; border-right: 1px solid #E5E7EB;
    }
    [data-testid="stSidebarUserContent"] { padding-top: 4.5rem; }
    /* Streamlit's own toolbar (Deploy button, hamburger) sits at the very
       top and would otherwise visually collide with our fixed navbar --
       it's replaced by the navbar, not just covered by it. */
    [data-testid="stHeader"] { display: none; }
    [data-testid="stMetricValue"] { color: #2563EB; }
    .stButton > button[kind="primary"] { background-color: #2563EB; border-color: #2563EB; }

    .stApp h1 { color: #1E3A8A; }
    .stApp h2, .stApp h3 {
        color: #1E3A8A;
        border-left: 3px solid #2563EB;
        padding-left: 0.8rem;
        margin-top: 1.8rem;
    }

    /* Fixed to the actual viewport edges -- unlike percentage/vw-based
       "breakout" margins, position:fixed isn't affected by the sidebar
       pushing the main content column off-center, so this reliably spans
       the true full browser width including over the sidebar. */
    .luad-navbar, .luad-footer {
        position: fixed;
        left: 0;
        width: 100vw;
        z-index: 999999;
        box-sizing: border-box;
    }

    .luad-navbar {
        top: 0;
        background: #1E3A8A;
        padding: 0.9rem 2rem;
        display: flex;
        justify-content: space-between;
        align-items: center;
        flex-wrap: wrap;
        gap: 0.6rem;
    }
    .luad-navbar .luad-brand { color: #FFFFFF; font-weight: 800; font-size: 1.1rem; }
    .luad-navbar .luad-brand-divider { color: rgba(255,255,255,0.35); margin: 0 0.65rem; font-weight: 300; }
    .luad-navbar .luad-brand-sub { color: #93C5FD; font-weight: 400; font-size: 0.86rem; }
    .luad-navbar nav a {
        color: #BFDBFE; text-decoration: none; font-size: 0.86rem;
        padding: 0.3rem 0.7rem; border-radius: 999px; transition: background 0.15s, color 0.15s;
    }
    .luad-navbar nav a:hover { color: #FFFFFF; background: rgba(255,255,255,0.12); }
    .luad-navbar nav a.luad-nav-external { color: #93C5FD; }

    .luad-footer {
        bottom: 0;
        background: #FFFFFF;
        border-top: 1px solid #E2E8F0;
        padding: 0.55rem 2rem;
        color: #64748B;
        font-size: 0.78rem;
        line-height: 1.4;
    }
    .luad-footer a { color: #2563EB; text-decoration: none; font-weight: 600; }
    .luad-footer a:hover { text-decoration: underline; }

    /* So the nav links' anchor-jump doesn't land the target section
       right underneath the fixed navbar. */
    #home, #analysis, #results, #about { scroll-margin-top: 4.5rem; }

    .stTabs [data-baseweb="tab-list"] { gap: 1.6rem; }
    .stTabs [data-baseweb="tab"] { color: #6B7280; font-weight: 500; }
    .stTabs [data-baseweb="tab-list"] button:nth-child(1)[aria-selected="true"] { color: #2563EB !important; }
    .stTabs [data-baseweb="tab-list"] button:nth-child(2)[aria-selected="true"] { color: #0F766E !important; }
    .stTabs [data-baseweb="tab-highlight"] { background-color: #2563EB; }

    [data-testid="stExpander"] {
        border: 1px solid #E5E7EB !important;
        border-radius: 12px !important;
        box-shadow: 0 1px 3px rgba(15, 23, 42, 0.04);
    }

    [data-testid="stDataFrame"] {
        border: 1px solid #E5E7EB;
        border-radius: 10px;
        overflow: hidden;
        box-shadow: 0 1px 3px rgba(15, 23, 42, 0.04);
    }

    [data-testid="stFileUploader"] {
        border: 1px solid #E5E7EB;
        border-radius: 10px;
        padding: 0.5rem 0.6rem;
        background: #F8FAFC;
        margin-bottom: 0.5rem;
    }
    [data-testid="stFileUploader"] label p { font-size: 0.85rem; font-weight: 600; }
    [data-testid="stFileUploaderDropzone"] {
        background: #FFFFFF !important;
        border: 1px dashed #CBD5E1 !important;
        border-radius: 8px !important;
        min-height: 2.6rem !important;
    }

    .luad-card {
        border: 1px solid #E5E7EB;
        border-radius: 12px;
        padding: 1rem 1.25rem;
        box-shadow: 0 1px 3px rgba(15, 23, 42, 0.05);
        background: #FFFFFF;
        margin-bottom: 0.9rem;
    }
    .luad-card--blue { border-left: 3px solid #2563EB; }
    .luad-card--flow { background: linear-gradient(180deg, #EFF6FF 0%, #FFFFFF 65%); border-color: #DBEAFE; }
    .luad-card--teal { border-left: 3px solid #0F766E; }
    .luad-card--rank1 { border: 1.5px solid #2563EB; box-shadow: 0 4px 10px rgba(37, 99, 235, 0.12); }

    .luad-badge {
        display: inline-block; padding: 0.16rem 0.6rem; border-radius: 999px;
        font-size: 0.72rem; font-weight: 600; letter-spacing: 0.02em;
    }
    .luad-badge--fda { background: #DCFCE7; color: #15803D; }
    .luad-badge--strong { background: #CCFBF1; color: #0F766E; }
    .luad-badge--info { background: #DBEAFE; color: #2563EB; }
    .luad-badge--muted { background: #F1F5F9; color: #475569; }
    .luad-badge--warning { background: #FEF3C7; color: #B45309; }
    .luad-badge--danger { background: #FEE2E2; color: #B91C1C; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div class="luad-navbar">
      <div><span class="luad-brand">LUADtx</span><span class="luad-brand-divider">|</span><span class="luad-brand-sub">End-to-End Personalized Therapy Prioritization</span></div>
      <nav>
        <a href="#home">Home</a>
        <a href="#analysis">Analysis</a>
        <a href="#results">Results</a>
        <a href="#about">About</a>
        <a class="luad-nav-external" href="https://github.com/yujuan-zhang/luad2" target="_blank">GitHub ↗</a>
      </nav>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("<div id='home'></div>", unsafe_allow_html=True)
st.markdown("### Patient / Case")
st.markdown(
    """
    <div class="luad-card luad-card--blue">
      <span class="luad-badge luad-badge--info">TCGA-38-4627</span>
      <span class="luad-badge luad-badge--info">LUAD</span>
      <span class="luad-badge luad-badge--muted">Demo Case</span>
      <div style="margin-top:0.6rem; font-size:0.85rem;">
        <span style="color:#15803D;">✓ Somatic variants</span>&nbsp;&nbsp;
        <span style="color:#15803D;">✓ Tumor RNA</span>&nbsp;&nbsp;
        <span class="luad-badge luad-badge--warning">HLA Synthetic</span>
      </div>
    </div>
    """,
    unsafe_allow_html=True,
)


@st.cache_data
def load_precomputed_demo():
    return json.loads(PRECOMPUTED_PATH.read_text())


st.sidebar.markdown("#### Input Data")
st.sidebar.caption("Leave blank to use the built-in TCGA-38-4627 case.")
vcf_file = st.sidebar.file_uploader("① Somatic VCF", type=["gz", "vcf"])
expression_file = st.sidebar.file_uploader("② Tumor Expression", type=["gz", "tsv"])
hla_file = st.sidebar.file_uploader("③ HLA Typing", type=["tsv"])
has_upload = bool(vcf_file or expression_file or hla_file)

if st.sidebar.button("Run Analysis →", type="primary"):
    if not has_upload:
        # Demo case's answer doesn't change between runs -- load the
        # precomputed result instead of re-running the ~1 minute real MHC
        # binding prediction.
        st.session_state["result"] = load_precomputed_demo()
    else:
        files = {}
        if vcf_file:
            files["vcf"] = (vcf_file.name, vcf_file.getvalue())
        if expression_file:
            files["expression"] = (expression_file.name, expression_file.getvalue())
        if hla_file:
            files["hla"] = (hla_file.name, hla_file.getvalue())
        try:
            with st.spinner("Running analysis... (real MHC binding prediction takes about a minute)"):
                resp = requests.post(f"{API_URL}/analyze", files=files, timeout=300)
            resp.raise_for_status()
            st.session_state["result"] = resp.json()
        except requests.RequestException:
            st.error(
                "Couldn't reach the analysis backend. Live analysis of uploaded files needs "
                "the FastAPI backend running locally (`uvicorn backend.main:app`) -- it isn't "
                "available on this hosted demo, which only shows the precomputed TCGA-38-4627 case."
            )

# Show the demo result on first load, without requiring a click.
if "result" not in st.session_state:
    st.session_state["result"] = load_precomputed_demo()

result = st.session_state["result"]

# ── Analysis Flow: horizontal KPI strip ──────────────────────────────────────
st.markdown("<div id='analysis'></div>", unsafe_allow_html=True)
st.subheader("Analysis Flow")
stages = list(result["funnel"].items())
n = len(stages)
flow_html = "<div class='luad-card luad-card--flow' style='display:flex; align-items:flex-start;'>"
for i, (label, value) in enumerate(stages):
    flow_html += (
        "<div style='flex:1; text-align:center;'>"
        f"<div style='font-size:0.8rem; opacity:0.65; line-height:1.25; min-height:2.1em;'>{html.escape(label)}</div>"
        f"<div style='font-size:1.9rem; font-weight:700; line-height:1.1; color:#2563EB;'>{value:,}</div>"
        "</div>"
    )
    if i < n - 1:
        flow_html += "<div style='font-size:1.6rem; padding-top:2.1rem; color:#93C5FD;'>→</div>"
flow_html += "</div>"
st.markdown(flow_html, unsafe_allow_html=True)
with st.expander("ⓘ Filtering criteria"):
    st.markdown(
        "- **Expressed:** TPM ≥ 1\n"
        "- **Variant-derived Peptides:** real UniProt sequence + missense substitution\n"
        "- **Peptide-HLA Evaluations:** candidate 8-11mer windows × HLA alleles evaluated\n"
        "- **Presented Candidates:** IC50 ≤ 500nM"
    )

# ── Key Clinical Finding ──────────────────────────────────────────────────────
st.markdown("<div id='results'></div>", unsafe_allow_html=True)
st.subheader("Key Clinical Finding")
drug_matches = result["drug_matches"]
if not drug_matches:
    st.info("No actionable targeted-therapy driver in this case's variants.")
else:
    drug_df = pd.DataFrame(drug_matches)
    drug_df["_rank"] = drug_df["level"].map(_LEVEL_ORDER).fillna(9)
    for (gene, change), group in drug_df.groupby(["gene", "protein_change"], sort=False):
        group = group.sort_values("_rank")
        targeted = group[group["drug"].apply(_is_targeted_drug)]
        other = group[~group["drug"].apply(_is_targeted_drug)]

        targeted_drugs = list(dict.fromkeys(targeted["drug"]))  # unique, evidence-ranked order
        preferred, alternatives = (targeted_drugs[0], targeted_drugs[1:3]) if targeted_drugs else (None, [])
        more_drugs = targeted_drugs[3:]
        top_level = targeted["level"].iloc[0] if not targeted.empty else group["level"].iloc[0]
        badge_class = _level_badge_class(top_level)

        alt_html = (
            f"<div style='margin-top:0.6rem; font-size:0.85rem; opacity:0.75;'>"
            f"<b>Alternative evidence</b> &nbsp; {html.escape(' · '.join(alternatives))}</div>"
            if alternatives else ""
        )
        st.markdown(
            "<div class='luad-card luad-card--blue'>"
            f"<div style='display:flex; justify-content:space-between; align-items:center;'>"
            f"<div style='font-size:1.1rem; font-weight:700;'>{html.escape(gene)} {html.escape(change)}</div>"
            f"<span class='luad-badge {badge_class}'>{html.escape(top_level)}</span>"
            "</div>"
            "<div style='font-size:0.8rem; opacity:0.6; margin-top:0.1rem;'>Actionable Driver</div>"
            "<div style='margin-top:0.7rem; font-size:0.78rem; opacity:0.6; text-transform:uppercase; letter-spacing:0.03em;'>Preferred therapy</div>"
            f"<div style='font-size:1.3rem; font-weight:700; color:#1E3A8A;'>{html.escape(preferred) if preferred else '—'}</div>"
            f"{alt_html}"
            "</div>",
            unsafe_allow_html=True,
        )

        with st.expander("View all evidence →"):
            if more_drugs:
                st.markdown(f"**Additional targeted-therapy evidence:** {', '.join(more_drugs)}")
            if not other.empty:
                other_drugs = list(dict.fromkeys(other["drug"]))
                st.markdown(f"**Other systemic therapy evidence (not targeted):** {', '.join(other_drugs)}")
            st.markdown(_drug_table_html(group.drop(columns=["_rank"]).to_dict("records")), unsafe_allow_html=True)

# ── Tabs: the two main analysis branches ─────────────────────────────────────
tab_therapy, tab_neo = st.tabs(["Targeted Therapy", "Neoantigen / INT"])

with tab_therapy:
    if not drug_matches:
        st.info("No targeted drug matches found.")
    else:
        st.markdown(_drug_table_html(drug_matches), unsafe_allow_html=True)

with tab_neo:
    neoantigens = result["neoantigens"]
    if not neoantigens:
        st.info("No HLA-presented neoantigens passed the filtering thresholds.")
    else:
        st.markdown("##### Top 5 Neoantigens")
        st.caption("Ranked by IC50 (lower = stronger predicted HLA binding). No ML ranking model yet -- rank here is by IC50 only.")
        top5 = neoantigens[:5]
        cols = st.columns(len(top5))
        for rank, (col, epitope) in enumerate(zip(cols, top5), start=1):
            badge_class, badge_text = _ic50_badge(epitope["ic50_nm"])
            vaf = epitope.get("dna_vaf")
            vaf_text = f"{vaf:.2f}" if vaf is not None else "—"
            rank_class = " luad-card--rank1" if rank == 1 else ""
            card_html = (
                f"<div class='luad-card luad-card--teal{rank_class}'>"
                f"<div style='font-size:0.72rem; opacity:0.6; font-weight:700; letter-spacing:0.03em;'>RANK #{rank}</div>"
                f"<div style='font-size:1rem; font-weight:700; margin:0.15rem 0 0.5rem;'>"
                f"{html.escape(epitope['gene'])} {html.escape(epitope['protein_change'])}</div>"
                f"<div style='font-family:monospace; font-size:0.82rem; background:#F1F5F9; "
                f"border-radius:6px; padding:0.3rem 0.45rem; margin-bottom:0.5rem; word-break:break-all;'>"
                f"{html.escape(epitope['peptide'])}</div>"
                f"<span class='luad-badge {badge_class}'>{badge_text} · {epitope['ic50_nm']} nM</span>"
                f"<div style='font-size:0.8rem; opacity:0.75; margin-top:0.5rem;'>HLA {html.escape(epitope['hla_allele'])}</div>"
                f"<div style='font-size:0.8rem; opacity:0.75;'>TPM {epitope['tumor_tpm']:.1f} · VAF {vaf_text}</div>"
                "</div>"
            )
            with col:
                st.markdown(card_html, unsafe_allow_html=True)

        st.markdown(f"##### Full Ranking ({len(neoantigens)})")
        st.dataframe(pd.DataFrame(neoantigens), use_container_width=True, hide_index=True)

# ── Pathways: summary cards first, KEGG diagram only on expand ──────────────
st.subheader("Affected Pathways")
pathways = result.get("pathways", [])
if not pathways:
    st.info("No core LUAD pathway was hit (no variant falls in these 8 pathways' genes).")
else:
    st.caption("Green = mutated gene | Red = highly expressed (TPM ≥ 5) | Yellow = expressed (TPM ≥ 1)")

    summary_cols = st.columns(len(pathways))
    for col, pw in zip(summary_cols, pathways):
        n_genes = len(pw["mutated_hit_genes"])
        with col:
            st.markdown(
                "<div class='luad-card luad-card--blue' style='text-align:center; padding:0.8rem 0.5rem;'>"
                f"<div style='font-weight:700; font-size:0.9rem;'>{html.escape(pw['name'])}</div>"
                f"<span class='luad-badge luad-badge--info' style='margin-top:0.4rem;'>{n_genes} gene{'s' if n_genes != 1 else ''}</span>"
                "</div>",
                unsafe_allow_html=True,
            )

    for pw in pathways:
        with st.expander(f"{pw['name']} — {', '.join(pw['mutated_hit_genes'])}", expanded=False):
            if pw["image_png_base64"]:
                st.image(base64.b64decode(pw["image_png_base64"]), **{_IMAGE_WIDTH_KWARG: True})
            st.dataframe(pd.DataFrame(pw["gene_table"]), use_container_width=True, hide_index=True)
            st.caption(f"[KEGG's own colored pathway link (fallback/reference)]({pw['kegg_url']})")

# ── Supporting detail ─────────────────────────────────────────────────────────
with st.expander(f"All Somatic Variants ({len(result['variants'])})"):
    st.dataframe(pd.DataFrame(result["variants"]), use_container_width=True, hide_index=True)

# ── About ──────────────────────────────────────────────────────────────────
st.markdown("<div id='about'></div>", unsafe_allow_html=True)
st.subheader("About")
st.markdown(
    "<div class='luad-card'>"
    "LUADtx is an end-to-end precision oncology platform for lung adenocarcinoma, integrating "
    "somatic variant analysis, tumor expression, HLA typing, targeted-therapy evidence, pathway "
    "interpretation, and neoantigen prioritization. The platform combines Nextflow/nf-core, "
    "<a href='https://civicdb.org' target='_blank'>CIViC</a>, KEGG, UniProt, pVACtools/MHCflurry, "
    "FastAPI, Streamlit, Docker, and AWS into a reproducible research workflow."
    "</div>",
    unsafe_allow_html=True,
)

# ── Footer ─────────────────────────────────────────────────────────────────
st.markdown(
    """
    <div class="luad-footer">
      Yujuan Zhang, PhD | Professor<br>
      <a href="https://github.com/yujuan-zhang" target="_blank">GitHub</a> · <a href="https://github.com/yujuan-zhang/luad2" target="_blank">⭐ Star this project</a><br>
      © 2026 Yujuan Zhang. All rights reserved.
    </div>
    """,
    unsafe_allow_html=True,
)
