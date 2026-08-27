import base64
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

st.set_page_config(page_title="LUAD Neoantigen & Targeted Therapy", layout="wide")

# Streamlit has no set_page_config knob for sidebar width, so narrow it with CSS.
st.markdown(
    """
    <style>
    [data-testid="stSidebar"] { min-width: 230px; max-width: 230px; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("LUAD Neoantigen & Targeted Therapy")
st.caption("Case: TCGA-38-4627 — real somatic variants + real tumor expression, synthetic HLA typing")


@st.cache_data
def load_precomputed_demo():
    return json.loads(PRECOMPUTED_PATH.read_text())


st.sidebar.header("Upload data (leave blank to use the built-in TCGA-38-4627 case)")
vcf_file = st.sidebar.file_uploader("Somatic VCF (variants.vcf.gz)", type=["gz", "vcf"])
expression_file = st.sidebar.file_uploader("Tumor expression (expression.tsv.gz)", type=["gz", "tsv"])
hla_file = st.sidebar.file_uploader("HLA typing (hla.tsv)", type=["tsv"])
has_upload = bool(vcf_file or expression_file or hla_file)

if st.sidebar.button("Run analysis", type="primary"):
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

# ── Pipeline Funnel: horizontal KPI strip ────────────────────────────────────
st.subheader("Pipeline Funnel")
stages = list(result["funnel"].items())
n = len(stages)
cols = st.columns([4, 1] * (n - 1) + [4])
for i, (label, value) in enumerate(stages):
    with cols[i * 2]:
        # st.metric's built-in label truncates with an ellipsis instead of
        # wrapping when the column is this narrow -- a plain markdown label
        # above the number wraps onto a second line instead.
        st.markdown(
            f"<div style='text-align:center;'>"
            f"<div style='font-size:0.8rem; opacity:0.65; line-height:1.25; min-height:2.1em;'>{label}</div>"
            f"<div style='font-size:1.9rem; font-weight:600; line-height:1.1;'>{value}</div>"
            f"</div>",
            unsafe_allow_html=True,
        )
    if i < n - 1:
        with cols[i * 2 + 1]:
            st.markdown(
                "<div style='text-align:center; font-size:1.6rem; padding-top:2.1rem; opacity:0.5;'>→</div>",
                unsafe_allow_html=True,
            )
st.caption("Expressed: TPM ≥ 1  |  HLA-presented: IC50 ≤ 500nM")

# ── Key Clinical Finding ──────────────────────────────────────────────────────
st.subheader("Key Clinical Finding")
drug_matches = result["drug_matches"]
if not drug_matches:
    st.info("No actionable targeted-therapy driver in this case's variants.")
else:
    drug_df = pd.DataFrame(drug_matches)
    drug_df["_rank"] = drug_df["level"].map(_LEVEL_ORDER).fillna(9)
    for (gene, change), group in drug_df.groupby(["gene", "protein_change"], sort=False):
        group = group.sort_values("_rank")
        drugs = list(dict.fromkeys(group["drug"]))  # unique, evidence-ranked order
        levels = ", ".join(sorted(group["level"].unique(), key=lambda lv: _LEVEL_ORDER.get(lv, 9)))
        sources = ", ".join(sorted(group["source"].unique()))
        with st.container(border=True):
            st.markdown(f"#### {gene} {change} — Actionable Driver")
            st.markdown(f"**Recommended targeted therapies:** {' / '.join(drugs)}")
            st.caption(f"Evidence: {levels}  ·  Source: {sources}")

# ── Tabs: the two main analysis branches ─────────────────────────────────────
tab_therapy, tab_neo = st.tabs(["Targeted Therapy", "Neoantigen / INT"])

with tab_therapy:
    if not drug_matches:
        st.info("No targeted drug matches found.")
    else:
        st.dataframe(pd.DataFrame(drug_matches), use_container_width=True)

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
            with col:
                with st.container(border=True):
                    st.markdown(f"**#{rank} {epitope['gene']} {epitope['protein_change']}**")
                    st.code(epitope["peptide"], language=None)
                    st.caption(f"HLA {epitope['hla_allele']}")
                    st.metric("IC50 (nM)", epitope["ic50_nm"])
                    vaf = epitope.get("dna_vaf")
                    st.caption(f"TPM {epitope['tumor_tpm']:.1f}  ·  VAF {vaf:.2f}" if vaf is not None else f"TPM {epitope['tumor_tpm']:.1f}")

        st.markdown(f"##### Full Ranking ({len(neoantigens)})")
        st.dataframe(pd.DataFrame(neoantigens), use_container_width=True)

# ── Pathways: collapsed by default, expand for the KEGG diagram ─────────────
st.subheader("Pathways")
pathways = result.get("pathways", [])
if not pathways:
    st.info("No core LUAD pathway was hit (no variant falls in these 8 pathways' genes).")
else:
    st.caption("Green = mutated gene | Red = highly expressed (TPM ≥ 5) | Yellow = expressed (TPM ≥ 1)")
    for pw in pathways:
        with st.expander(f"{pw['name']} — {', '.join(pw['mutated_hit_genes'])}", expanded=False):
            if pw["image_png_base64"]:
                st.image(base64.b64decode(pw["image_png_base64"]), **{_IMAGE_WIDTH_KWARG: True})
            st.dataframe(pd.DataFrame(pw["gene_table"]), use_container_width=True)
            st.caption(f"[KEGG's own colored pathway link (fallback/reference)]({pw['kegg_url']})")

# ── Supporting detail ─────────────────────────────────────────────────────────
with st.expander(f"All Somatic Variants ({len(result['variants'])})"):
    st.dataframe(pd.DataFrame(result["variants"]), use_container_width=True)
