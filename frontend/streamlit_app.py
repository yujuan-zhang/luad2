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
CASE_METADATA_PATH = Path(__file__).resolve().parent.parent / "data" / "demo" / "case_metadata.json"

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
    database internals rather than a clinical summary. No vertical grid
    lines, light header background, hover row highlight -- reads as a
    product table rather than a spreadsheet dump."""
    rows = []
    for m in matches:
        level = m.get("level") or "—"
        line = m.get("line") or "—"
        resp_label, resp_class = _response_label(m.get("significance"))
        rows.append(
            "<tr>"
            f"<td style='padding:0.5rem 0.7rem;'>{html.escape(m['gene'])}</td>"
            f"<td style='padding:0.5rem 0.7rem;'>{html.escape(m['protein_change'])}</td>"
            f"<td style='padding:0.5rem 0.7rem; font-weight:600;'>{html.escape(m['drug'])}</td>"
            f"<td style='padding:0.5rem 0.7rem;'><span class='luad-badge {_level_badge_class(level)}'>{html.escape(level)}</span></td>"
            f"<td style='padding:0.5rem 0.7rem; opacity:0.75;'>{html.escape(line)}</td>"
            f"<td style='padding:0.5rem 0.7rem;'><span class='luad-badge {resp_class}'>{resp_label}</span></td>"
            "</tr>"
        )
    header = (
        "<tr>"
        "<th style='padding:0.45rem 0.7rem;'>Gene</th><th style='padding:0.45rem 0.7rem;'>Mutation</th>"
        "<th style='padding:0.45rem 0.7rem;'>Drug</th><th style='padding:0.45rem 0.7rem;'>Level</th>"
        "<th style='padding:0.45rem 0.7rem;'>Line</th><th style='padding:0.45rem 0.7rem;'>Response</th></tr>"
    )
    return (
        "<div class='luad-card luad-evidence-table' style='padding:0.4rem 0.2rem;'>"
        "<table style='width:100%; border-collapse:collapse; font-size:0.85rem;'>"
        f"<thead>{header}</thead><tbody>{''.join(rows)}</tbody></table>"
        "</div>"
    )


st.set_page_config(page_title="LUADtx", layout="wide")

st.markdown(
    """
    <style>
    .stApp { background-color: #F8FAFC; }
    /* Room for the fixed navbar (see below) so page content and the
       sidebar's own content don't sit underneath it. The footer is no
       longer fixed (see .luad-footer), so no bottom padding needed here. */
    .main .block-container {
        padding-top: 6.2rem; padding-bottom: 1.5rem;
        padding-left: 2.5rem; padding-right: 2.5rem;
    }
    [data-testid="stSidebar"] {
        min-width: 172px; max-width: 172px;
        background-color: #FFFFFF; border-right: 1px solid #E5E7EB;
    }
    [data-testid="stSidebarUserContent"] { padding-top: 6.2rem; }
    /* Thinner scrollbar in the sidebar instead of the browser's default
       thick one. */
    [data-testid="stSidebar"] { scrollbar-width: thin; scrollbar-color: #CBD5E1 transparent; }
    [data-testid="stSidebar"] ::-webkit-scrollbar { width: 6px; }
    [data-testid="stSidebar"] ::-webkit-scrollbar-thumb { background: #CBD5E1; border-radius: 3px; }
    [data-testid="stSidebar"] ::-webkit-scrollbar-thumb:hover { background: #94A3B8; }
    /* Streamlit's own toolbar (Deploy button, hamburger, sidebar-collapse
       arrow) sits at the very top, in the same space our fixed navbar
       occupies -- keeping it visible would mean either it's painted over
       and unclickable (our navbar is on top) or the navbar has to start
       lower to leave it room, which reopens the whole fixed-positioning
       math this page already had two rounds of trouble with. Hiding it
       entirely was the deliberate trade-off; restoring just the sidebar
       toggle needs that restructure, not a one-line fix. */
    [data-testid="stHeader"] { display: none; }
    [data-testid="stMetricValue"] { color: #2563EB; }
    .stButton > button[kind="primary"] {
        background-color: #EFF6FF;
        border: 1.5px solid #2563EB;
        color: #1E3A8A;
        font-weight: 700;
        position: relative;
        padding-top: 0.35rem; padding-bottom: 0.35rem;
    }
    .stButton > button[kind="primary"]::before {
        content: ""; position: absolute; top: 6px; left: 10px;
        width: 6px; height: 6px; border-radius: 50%; background: #2563EB;
    }
    .stButton > button[kind="secondary"] { padding-top: 0.35rem; padding-bottom: 0.35rem; }

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
       the true full browser width including over the sidebar. Only the
       navbar is fixed/sticky; the footer is a normal trailing element (see
       .luad-footer) so it doesn't sit on top of content while scrolling. */
    .luad-navbar {
        position: fixed;
        top: 0;
        left: 0;
        width: 100vw;
        z-index: 999999;
        box-sizing: border-box;
        background: #1E3A8A;
        padding: 1.3rem 2.6rem;
        display: flex;
        justify-content: space-between;
        align-items: center;
        flex-wrap: wrap;
        gap: 0.6rem;
    }
    .luad-navbar .luad-brand { color: #FFFFFF; font-weight: 900; font-size: 2rem; }
    .luad-navbar .luad-brand-divider { color: rgba(255,255,255,0.35); margin: 0 0.85rem; font-weight: 300; }
    .luad-navbar .luad-brand-sub { color: #93C5FD; font-weight: 600; font-size: 1.25rem; opacity: 0.9; }
    .luad-navbar nav a {
        color: #BFDBFE; text-decoration: none; font-size: 1.2rem; font-weight: 600;
        padding: 0.5rem 1.1rem; border-radius: 999px; transition: background 0.15s, color 0.15s;
    }
    .luad-navbar nav a:hover { color: #FFFFFF; background: rgba(255,255,255,0.12); }
    .luad-navbar nav a.luad-nav-external { color: #93C5FD; }

    .luad-footer {
        margin-top: 2rem;
        border-top: 1px solid #E2E8F0;
        padding: 0.6rem 0;
        color: #B0BAC7;
        font-size: 0.71rem;
        line-height: 1.3;
        text-align: center;
    }
    .luad-footer a { color: #2563EB; text-decoration: none; font-weight: 600; }
    .luad-footer a:hover { text-decoration: underline; }

    /* So the nav links' anchor-jump doesn't land the target section
       right underneath the fixed navbar. */
    #home, #analysis, #results, #about { scroll-margin-top: 6.2rem; }

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

    /* One border (the dropzone's own), not a wrapper box around it too --
       double-boxing read as more "raw widget" clutter, not less. */
    [data-testid="stFileUploader"] { margin-bottom: 0.3rem; }
    [data-testid="stFileUploader"] label p { font-size: 0.8rem; font-weight: 600; }
    /* "Drag and drop file here / Limit 200MB per file..." -- redundant
       now that the label itself says the format and size limit. */
    [data-testid="stFileUploaderDropzoneInstructions"] { display: none; }
    [data-testid="stFileUploaderDropzone"] {
        background: #F8FAFC !important;
        border: 1px dashed #CBD5E1 !important;
        border-radius: 8px !important;
        min-height: 2.4rem !important;
        padding: 0.4rem !important;
    }

    .luad-card {
        border: 1px solid #E5E7EB;
        border-radius: 12px;
        padding: 0.8rem 1.25rem;
        box-shadow: 0 1px 3px rgba(15, 23, 42, 0.05);
        background: #FFFFFF;
        margin-bottom: 0.7rem;
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
    .luad-badge--success { background: #DCFCE7; color: #15803D; }

    /* Result cards (Patient/Case, Key Clinical Finding) sit right under a
       section header that already has a blue left bar -- repeating the
       same bar on the card read as redundant, so these get a thin top
       accent instead. */
    .luad-card--accent-top { border-top: 3px solid #2563EB; }

    .luad-evidence-table thead tr { background: #F8FAFC; text-align: left; }
    .luad-evidence-table th { opacity: 0.55; font-size: 0.78rem; font-weight: 600; }
    .luad-evidence-table tbody tr { border-top: 1px solid #F1F5F9; }
    .luad-evidence-table tbody tr:hover { background: #EFF6FF; }

    .luad-hero { padding: 0 0 0.5rem; }
    .luad-hero h2 { margin: 0 0 0.2rem; font-size: 1.5rem; color: #1E3A8A; border: none; padding: 0; }
    .luad-hero p { margin: 0; font-size: 0.92rem; color: #64748B; }

    /* Best-effort "active" nav link: CSS :target only updates when a nav
       link is actually clicked (no JS scroll-spy), so this reflects the
       last section jumped to, not necessarily what's on screen. */
    /* Home reads as active by default (nothing else is a real page you can
       "leave" here) and steps aside once another section is targeted. */
    .luad-navbar nav a[href="#home"] { color: #FFFFFF; box-shadow: inset 0 -2px 0 #FFFFFF; }
    body:has(#analysis:target) .luad-navbar nav a[href="#home"],
    body:has(#results:target) .luad-navbar nav a[href="#home"],
    body:has(#about:target) .luad-navbar nav a[href="#home"] {
        color: #BFDBFE; box-shadow: none;
    }
    body:has(#analysis:target) .luad-navbar nav a[href="#analysis"],
    body:has(#results:target) .luad-navbar nav a[href="#results"],
    body:has(#about:target) .luad-navbar nav a[href="#about"] {
        color: #FFFFFF; box-shadow: inset 0 -2px 0 #FFFFFF;
    }

    @media print {
        .luad-navbar { position: static; width: 100%; }
        .main .block-container { padding-top: 0.5rem; }
    }
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
        <a class="luad-nav-external" href="https://github.com/yujuan-zhang/luadtx" target="_blank">GitHub ↗</a>
      </nav>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("<div id='home'></div>", unsafe_allow_html=True)
st.markdown(
    """
    <div class="luad-hero">
      <h2>Personalized Therapy Prioritization for LUAD</h2>
      <p>From genomic and transcriptomic profiles to targeted therapy and neoantigen prioritization.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("### Patient / Case")
_case_meta = json.loads(CASE_METADATA_PATH.read_text())
_clinical = _case_meta.get("clinical", {})
_clinical_line = " &nbsp;·&nbsp; ".join(
    html.escape(str(v)) for v in [
        _clinical.get("sex"),
        f"Age {_clinical.get('age_at_diagnosis')}" if _clinical.get("age_at_diagnosis") is not None else None,
        _clinical.get("stage"),
        f"Vital status: {_clinical.get('vital_status')}" if _clinical.get("vital_status") else None,
        f"Smoking history: {_clinical.get('smoking_history')}" if _clinical.get("smoking_history") else None,
    ] if v
)
st.markdown(
    f"""
    <div class="luad-card luad-card--accent-top">
      <span class="luad-badge luad-badge--info">TCGA-38-4627</span>
      <span class="luad-badge luad-badge--info">LUAD</span>
      <span class="luad-badge luad-badge--muted">Demo Case</span>
      <div style="margin-top:0.5rem; font-size:0.85rem; opacity:0.8;">{_clinical_line}</div>
      <div style="margin-top:0.6rem;">
        <span class="luad-badge luad-badge--success">✓ Somatic VCF</span>
        <span class="luad-badge luad-badge--success">✓ Tumor RNA</span>
        <span class="luad-badge luad-badge--success">✓ Real Clinical (GDC)</span>
        <span class="luad-badge luad-badge--warning">⚠ Synthetic HLA</span>
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
vcf_file = st.sidebar.file_uploader("① Somatic VCF — GZ/VCF, max 200MB", type=["gz", "vcf"])
expression_file = st.sidebar.file_uploader("② Tumor Expression — GZ/TSV, max 200MB", type=["gz", "tsv"])
hla_file = st.sidebar.file_uploader("③ HLA Typing — TSV, max 200MB", type=["tsv"])
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
        "<div style='height:2.3em; display:flex; align-items:flex-end; justify-content:center; "
        f"font-size:0.8rem; opacity:0.65; line-height:1.15;'>{html.escape(label)}</div>"
        f"<div style='font-size:1.9rem; font-weight:700; line-height:1.1; margin-top:0.3rem; color:#2563EB;'>{value:,}</div>"
        "</div>"
    )
    if i < n - 1:
        flow_html += "<div style='font-size:1.6rem; padding-top:2.6rem; color:#BFDBFE;'>→</div>"
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
            "<div class='luad-card luad-card--accent-top'>"
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

# ── Pathways: click a summary card to select it, one detail view below ──────
st.subheader("Affected Pathways")
pathways = result.get("pathways", [])
if not pathways:
    st.info("No pathway was hit (no variant falls in any of the 79 cancer-related KEGG pathways tracked).")
else:
    st.markdown(
        "<div style='font-size:0.82rem; color:#64748B; margin-bottom:0.6rem;'>"
        "<span style='color:#22C55E;'>&#9632;</span> Mutated &nbsp;&nbsp;"
        "<span style='color:#EF4444;'>&#9632;</span> High expression (TPM ≥ 5) &nbsp;&nbsp;"
        "<span style='color:#EAB308;'>&#9632;</span> Expressed (TPM ≥ 1)"
        "</div>",
        unsafe_allow_html=True,
    )

    if "selected_pathway" not in st.session_state:
        st.session_state["selected_pathway"] = pathways[0]["pathway_id"]

    # pathways is already sorted by hit-gene count descending (get_hit_pathways) --
    # a quick-select row of the top few, plus a dropdown for the full list, scales
    # to any hit count instead of squeezing one column per pathway (up to 79).
    TOP_N_BUTTONS = 6
    top_pathways = pathways[:TOP_N_BUTTONS]
    quick_cols = st.columns(len(top_pathways))
    for col, pw in zip(quick_cols, top_pathways):
        n_genes = len(pw["mutated_hit_genes"])
        is_active = st.session_state["selected_pathway"] == pw["pathway_id"]
        with col:
            if st.button(
                f"{pw['name']} · {n_genes} gene{'s' if n_genes != 1 else ''}",
                key=f"pwbtn_{pw['pathway_id']}",
                use_container_width=True,
                type="primary" if is_active else "secondary",
            ):
                st.session_state["selected_pathway"] = pw["pathway_id"]

    if len(pathways) > TOP_N_BUTTONS:
        all_ids = [p["pathway_id"] for p in pathways]
        labels = {p["pathway_id"]: f"{p['name']} · {len(p['mutated_hit_genes'])} gene(s)" for p in pathways}
        current_idx = all_ids.index(st.session_state["selected_pathway"]) if st.session_state["selected_pathway"] in all_ids else 0
        chosen = st.selectbox(
            f"Or choose from all {len(pathways)} affected pathways",
            options=all_ids,
            index=current_idx,
            format_func=lambda pid: labels[pid],
        )
        st.session_state["selected_pathway"] = chosen

    selected = next((p for p in pathways if p["pathway_id"] == st.session_state["selected_pathway"]), pathways[0])
    with st.container(border=True):
        st.markdown(f"**{selected['name']} — {', '.join(selected['mutated_hit_genes'])}**")
        if selected["image_png_base64"]:
            st.image(base64.b64decode(selected["image_png_base64"]), **{_IMAGE_WIDTH_KWARG: True})
        st.dataframe(pd.DataFrame(selected["gene_table"]), use_container_width=True, hide_index=True)
        st.markdown(f"[KEGG's own colored pathway link (fallback/reference)]({selected['kegg_url']})")

# ── Supporting detail ─────────────────────────────────────────────────────────
with st.expander(f"All Somatic Variants ({len(result['variants'])})"):
    st.dataframe(pd.DataFrame(result["variants"]), use_container_width=True, hide_index=True)

# ── About ──────────────────────────────────────────────────────────────────
st.markdown("<div id='about'></div>", unsafe_allow_html=True)
st.subheader("About")
_ABOUT_STACK = [
    ("Workflow", ["Nextflow", "nf-core"]),
    ("Annotation", ["Ensembl VEP"]),
    ("Evidence", ["CIViC", "KEGG"]),
    ("Prediction", ["UniProt", "pVACtools", "MHCflurry"]),
    ("Deployment", ["FastAPI", "Streamlit", "Docker", "AWS"]),
]
about_col1, about_col2 = st.columns([3, 2])
with about_col1:
    st.markdown(
        "<div class='luad-card' style='height:100%;'>"
        "LUADtx is an end-to-end precision oncology platform for lung adenocarcinoma, spanning raw "
        "sequencing reads, somatic variant analysis, tumor expression, HLA typing, targeted-therapy "
        "evidence, pathway interpretation, and neoantigen prioritization. Raw-read processing with "
        "Nextflow/nf-core runs on local HPC to reduce AWS compute costs, while downstream analysis "
        "and reporting are deployed on AWS in a hybrid HPC–cloud architecture. Developed as a "
        "reproducible research and portfolio platform."
        "</div>",
        unsafe_allow_html=True,
    )
with about_col2:
    stack_html = "<div class='luad-card' style='height:100%;'>"
    for group, tools in _ABOUT_STACK:
        tags = " ".join(
            f"<a href='https://civicdb.org' target='_blank' style='text-decoration:none;'>"
            f"<span class='luad-badge luad-badge--muted'>{html.escape(t)}</span></a>"
            if t == "CIViC" else f"<span class='luad-badge luad-badge--muted'>{html.escape(t)}</span>"
            for t in tools
        )
        stack_html += (
            f"<div style='margin-bottom:0.5rem;'>"
            f"<div style='font-size:0.72rem; opacity:0.6; text-transform:uppercase; letter-spacing:0.03em; margin-bottom:0.25rem;'>{group}</div>"
            f"{tags}</div>"
        )
    stack_html += "</div>"
    st.markdown(stack_html, unsafe_allow_html=True)

# ── Footer ─────────────────────────────────────────────────────────────────
st.markdown(
    """
    <div class="luad-footer">
      Yujuan Zhang, PhD | Professor &nbsp;·&nbsp;
      <a href="https://github.com/yujuan-zhang" target="_blank">GitHub</a> &nbsp;·&nbsp;
      <a href="https://github.com/yujuan-zhang/luadtx" target="_blank">⭐ Star this project</a> &nbsp;·&nbsp;
      © 2026 Yujuan Zhang. All rights reserved.
    </div>
    """,
    unsafe_allow_html=True,
)
