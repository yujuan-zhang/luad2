import base64
import json
from pathlib import Path

import pandas as pd
import requests
import streamlit as st

API_URL = "http://localhost:8000"
PRECOMPUTED_PATH = Path(__file__).resolve().parent.parent / "data" / "demo" / "precomputed_result.json"

st.set_page_config(page_title="LUAD Neoantigen & Targeted Therapy", layout="wide")
st.title("LUAD Neoantigen & Targeted Therapy")
st.caption("Default case: TCGA-38-4627 — real somatic variants + real tumor expression, synthetic HLA typing")


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

st.subheader("Pipeline Funnel")
funnel_df = pd.DataFrame(result["funnel"].items(), columns=["Stage", "Count"])
st.dataframe(funnel_df, use_container_width=False, hide_index=True)
st.caption("Expressed: TPM ≥ 1  |  HLA-presented: IC50 ≤ 500nM")

st.subheader("Somatic Variants")
st.dataframe(pd.DataFrame(result["variants"]), use_container_width=True)

st.subheader("Targeted Drug Matches")
if result["drug_matches"]:
    st.dataframe(pd.DataFrame(result["drug_matches"]), use_container_width=True)
else:
    st.info("No targeted drug matches found.")

st.subheader("Neoantigen Ranking")
st.dataframe(pd.DataFrame(result["neoantigens"]), use_container_width=True)

st.subheader("Pathways")
pathways = result.get("pathways", [])
if not pathways:
    st.info("No core LUAD pathway was hit (no variant falls in these 8 pathways' genes).")
else:
    st.caption("Green = mutated gene | Red = highly expressed (TPM ≥ 5) | Yellow = expressed (TPM ≥ 1)")
    for pw in pathways:
        with st.expander(f"{pw['name']} — {', '.join(pw['mutated_hit_genes'])}", expanded=True):
            if pw["image_png_base64"]:
                st.image(base64.b64decode(pw["image_png_base64"]), use_column_width=True)
            st.dataframe(pd.DataFrame(pw["gene_table"]), use_container_width=True)
            st.caption(f"[KEGG's own colored pathway link (fallback/reference)]({pw['kegg_url']})")
