# LUAD Precision Platform

LUAD neoantigen / targeted-therapy prioritization pipeline. See project plan for full architecture.

## 当前阶段：Phase 2 — 真实病例数据 + mock peptide/binding

默认 demo 病例是 **TCGA-38-4627**（来自 `luad_workflow` 项目已经跑过的真实结果）：

- `data/demo/variants.vcf.gz`：真实的 25 个 protein-altering somatic variants，未注释的原始 VCF —— 符合项目设计（用户上传 VCF，AWS VEP API 负责注释）
- `data/demo/expression.tsv.gz`：真实的全基因组 tumor expression（TPM + GTEx-lung z-score/percentile）
- `data/demo/hla.tsv`：**synthetic HLA**（这个病例没有真实 HLA typing，按项目计划明确标注为 synthetic）

`vep.py` 还不是真的调 VEP API（AWS 那一步还没搭），而是内部查表注释 —— 但表里的内容是真实的 VEP 注释结果
（来自 `luad_workflow`，gene/consequence/HGVSp/DNA VAF/hotspot/functional impact 都是真的），只是注释这个
*动作* 现在是本地查表完成，接口（`annotate_variants(vcf_path)`）跟真的调 VEP API 完全一样，以后接上真实
VEP API 时只需要换内部实现。
`civic.py` 是真实实现：curated 药物知识库 + CIViC（civicdb.org）GraphQL API 实时查询，不需要 OncoKB token。
`pvactools.py` 里的 mutant peptide 序列和 MHC binding IC50 还是 **mock**（`MOCKPEP-` 开头、hash 出来的 IC50）——
真实实现的下一步是接 Ensembl REST（取真实蛋白序列）+ IEDB REST（真实 binding 预测），两个都免费不需要 token。

## Pipeline Funnel

`main.py` 现在会输出每一步筛掉多少变异，而不只是最终三张表：

```
Protein-altering variants    25
Actionable variants           1   -> 进 drug_matches 分支
Neoantigen candidates        24   -> 进 pvactools 分支
Expressed variants           17   (TPM >= 1)
Peptide-HLA pairs           102   (17 candidates x 6 HLA alleles)
HLA-presented                29   (IC50 <= 500nM)
```

## Pathway 可视化

`pathway.py`（同一套设计思路来自 `luad_workflow/modules/06_pathway/kegg_viewer.py`）：自己在
本地缓存的 KEGG 官方 PNG 上用 Pillow 叠色块，不调 pathview/cytoscape。通路成员基因用 gseapy 的
KEGG_2021_Human gene set（本地缓存，不联网）；底图 + 基因框坐标是提前从 KEGG REST API/KGML 下载好
缓存在 `pipelines/downstream/kegg_cache/pathways/` 的（这次是直接从 `luad_workflow` 拷贝过来的，
没有重新下载）。

固定检查 8 条 LUAD 核心通路（MAPK / PI3K-AKT / ErbB-EGFR / p53 / Cell Cycle / TGF-β / Wnt / VEGF），
只渲染有命中突变基因的通路。颜色（跟 `luad_workflow` 版本不同，没有差异表达倍数，只看是否表达/高表达）：

- 绿色 = 突变基因
- 黄色 = 有表达（TPM ≥ 1）
- 红色 = 较高表达（TPM ≥ 5）
- 一个基因框同时符合多种状态时，切成竖条分别染色

`build_kegg_url()` 保留了生成 KEGG 官网彩色链接的功能，作为备用/对照。

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
data/demo/                        默认病例 TCGA-38-4627：真实 VCF + 真实 expression + synthetic HLA
pipelines/downstream/              核心分析逻辑：vep.py / civic.py / pvactools.py / pathway.py / main.py（串联）
pipelines/downstream/kegg_cache/   KEGG 通路底图 PNG + 基因框坐标缓存（不联网）
backend/                           FastAPI，把 pipeline 包成 /analyze 接口
frontend/                          Streamlit 网页，调用 FastAPI 展示结果（含 pipeline funnel + 通路图）
```
