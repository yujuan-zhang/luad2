import pandas as pd
import requests
import streamlit as st

API_URL = "http://localhost:8000"

st.set_page_config(page_title="LUAD Neoantigen & Targeted Therapy", layout="wide")
st.title("LUAD Neoantigen & Targeted Therapy")
st.caption("Default case: TCGA-38-4627 — real somatic variants + real tumor expression, synthetic HLA typing")

st.sidebar.header("上传数据（留空则用内置 TCGA-38-4627 病例）")
vcf_file = st.sidebar.file_uploader("Somatic VCF (variants.vcf.gz)", type=["gz", "vcf"])
expression_file = st.sidebar.file_uploader("Tumor expression (expression.tsv.gz)", type=["gz", "tsv"])
hla_file = st.sidebar.file_uploader("HLA typing (hla.tsv)", type=["tsv"])

if st.sidebar.button("运行分析", type="primary"):
    files = {}
    if vcf_file:
        files["vcf"] = (vcf_file.name, vcf_file.getvalue())
    if expression_file:
        files["expression"] = (expression_file.name, expression_file.getvalue())
    if hla_file:
        files["hla"] = (hla_file.name, hla_file.getvalue())

    with st.spinner("分析中..."):
        resp = requests.post(f"{API_URL}/analyze", files=files)
    resp.raise_for_status()
    st.session_state["result"] = resp.json()

if "result" in st.session_state:
    result = st.session_state["result"]

    st.subheader("Pipeline Funnel")
    funnel = result["funnel"]
    cols = st.columns(len(funnel))
    for col, (label, value) in zip(cols, funnel.items()):
        col.metric(label, value)

    st.subheader("Somatic Variants")
    st.dataframe(pd.DataFrame(result["variants"]), use_container_width=True)

    st.subheader("Targeted Drug Matches")
    if result["drug_matches"]:
        st.dataframe(pd.DataFrame(result["drug_matches"]), use_container_width=True)
    else:
        st.info("没有匹配到已知靶点药物。")

    st.subheader("Neoantigen Ranking")
    st.dataframe(pd.DataFrame(result["neoantigens"]), use_container_width=True)
else:
    st.info("点击左侧「运行分析」查看结果（默认用内置 TCGA-38-4627 病例）。")
