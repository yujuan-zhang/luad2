# LUAD Precision Platform

LUAD neoantigen / targeted-therapy prioritization pipeline. See project plan for full architecture.

## 当前阶段：Phase 2 — 真实病例数据 + mock peptide/binding

默认 demo 病例是 **TCGA-38-4627**（来自 `luad_workflow` 项目已经跑过的真实结果）：

- `data/demo/variants.tsv.gz`：真实的 25 个 protein-altering somatic variants，VEP 注释（gene/consequence/HGVSp/DNA VAF/hotspot/functional impact 都是真的）
- `data/demo/expression.tsv.gz`：真实的全基因组 tumor expression（TPM + GTEx-lung z-score/percentile）
- `data/demo/hla.tsv`：**synthetic HLA**（这个病例没有真实 HLA typing，按项目计划明确标注为 synthetic）

`vep.py` 不调用真实 VEP，而是加载已经注释好的变异表（VEP 标注这一步在 `luad_workflow` 里已经做过了）。
`civic.py` 是真实实现：curated 药物知识库 + CIViC（civicdb.org）GraphQL API 实时查询，不需要 OncoKB token。
`pvactools.py` 里的 mutant peptide 序列和 MHC binding IC50 还是 **mock**（`MOCKPEP-` 开头、hash 出来的 IC50）——
真实实现的下一步是接 Ensembl REST（取真实蛋白序列）+ IEDB REST（真实 binding 预测），两个都免费不需要 token。

## Pipeline Funnel

`main.py` 现在会输出每一步筛掉多少变异，而不只是最终三张表：

```
Somatic variants (protein-altering)         25
Actionable variants (targeted therapy)       1   -> 进 drug_matches 分支
Neoantigen candidate variants               24   -> 进 pvactools 分支
Expression-supported variants (TPM >= 1)    17
Peptide-HLA pairs evaluated                102   (17 candidates x 6 HLA alleles)
HLA-presented neoantigens (IC50 <= 500nM)   ~29
```

## 快速开始

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 1. 命令行测试 downstream pipeline（不需要开服务）
python -m pipelines.downstream.main

# 2. 启动后端 API（新开一个终端）
uvicorn backend.main:app --reload --port 8000

# 3. 启动前端网页（再开一个终端）
streamlit run frontend/streamlit_app.py
```

浏览器打开 Streamlit 给出的地址（默认 http://localhost:8501），点击「运行分析」即可看到内置 demo 数据的结果。

## 目录结构

```
data/demo/            默认病例 TCGA-38-4627：真实 variants + 真实 expression + synthetic HLA
pipelines/downstream/ 核心分析逻辑：vep.py / civic.py / pvactools.py / main.py（串联）
backend/               FastAPI，把 pipeline 包成 /analyze 接口
frontend/               Streamlit 网页，调用 FastAPI 展示结果（含 pipeline funnel）
```
