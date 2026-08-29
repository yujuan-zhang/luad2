# LUAD Precision Platform

LUAD neoantigen / targeted-therapy prioritization pipeline. See project plan for full architecture.

## 当前阶段：Phase 3 — 全链路真实

默认 demo 病例是 **TCGA-38-4627**（来自 `luad_workflow` 项目已经跑过的真实结果）：

- `data/demo/variants.vcf.gz`：真实的 25 个 protein-altering somatic variants，未注释的原始 VCF，INFO 里带真实
  DNA VAF —— 符合项目设计（用户上传 VCF，VEP API 负责注释）
- `data/demo/expression.tsv.gz`：真实的全基因组 tumor expression（TPM + GTEx-lung z-score/percentile）
- `data/demo/hla.tsv`：**synthetic HLA**（这个病例没有真实 HLA typing，按项目计划明确标注为 synthetic）

`vep.py` 现在真的调 Ensembl VEP REST API（`rest.ensembl.org/vep/human/region`，免费不需要 token，GRCh38）
做注释，`canonical=1` 锁定 canonical transcript，`gene`/`consequence`/`functional_impact`（VEP 原生的
HIGH/MODERATE/LOW/MODIFIER 分级）都是接口直接返回的；`protein_change` 是从返回的 `amino_acids` +
`protein_start` 自己拼的（没用 `hgvs=1` 参数——这个参数在这个接口上会 500）。DNA VAF 不是 VEP 的概念，
是从输入 VCF 自己的 `INFO/VAF` 字段读出来的（跟真实变异检测流程一样，caller 出 VAF，VEP 只管注释）。
`hotspot` 也不是接口直接给的——COSMIC 共定位（`colocated_variants[].somatic`）几乎所有体细胞变异都命中，
不能当热点信号，所以换成一个小的、保守的 curated 真实 LUAD driver hotspot 表（`_KNOWN_LUAD_HOTSPOTS`），
跟 `civic.DRUG_KB` 是同一个套路。API 调用失败会直接抛错，不会静默退化成假注释。
`civic.py` 是真实实现：curated 药物知识库 + CIViC（civicdb.org）GraphQL API 实时查询，不需要 OncoKB token。

`pvactools.py` 现在也是真实实现了：
- **mutant peptide**：从 UniProt REST（真实蛋白序列，免费不需要 token）取该基因的 canonical 序列，在突变位置替换氨基酸，切出真实的 flanking peptide。只对 missense 变异做（stop_gained/frameshift/splice 会产生全新的下游序列，需要真正的 CDS 层建模才能算对，这里没做，会被跳过而不是编一个假的出来）。
- **MHC binding**：真的装了 `pvactools`（7.1.2），并用它依赖的 `mhcflurry` 真实预测模型算 IC50 —— 不是通过完整 `pvacseq run`（那条路需要真 VEP 标注 + Wildtype/Frameshift plugin，我们没有），而是直接批量调 `mhcflurry-predict`（跟 pVACtools 内部包装类调的是同一个模型）。用 CMV/流感的经典强结合表位验证过预测结果是对的。
- 装 mhcflurry 踩了个坑：它内部用的是老版 TF1 Keras API，新版 Keras 3 删掉了，需要装 `tf-keras` 兼容层 + 设 `TF_USE_LEGACY_KERAS=1`（`pvactools.py` 里已经处理了）。
- 首次使用前还需要手动跑一次模型下载（不在 pip 依赖里，是单独的模型文件，135MB）：
  ```bash
  mhcflurry-downloads fetch models_class1_presentation
  ```
- **vaccine construct（`design_vaccine_construct()`）**：把 top 5 个 neoantigen 拼成一条疫苗肽链，核心问题是
  两个肽拼接处可能意外产生一个新的、没设计过的强结合表位（junctional epitope）——这是 pVACtools 自带的
  `pvacvector` 工具要解决的事，但**没有直接调用它的 CLI**：实测它对每一个（HLA型别 × 表位长度 × spacer）
  组合都单独起一次进程重新加载 MHCflurry 模型，这个case（5个候选×6个HLA型别）跑一轮要 1.5 小时以上，是
  跟当初 pVACtools 自带 wrapper 类同一个性能问题（10+分钟 vs 批量调用几十秒），只是这次严重得多。用的是同一个
  解法：把所有候选拼接点（每对肽 × 每种spacer）产生的候选表位一次性批量丢给 `_run_mhcflurry`，再对5个候选的
  全排列（120种，直接暴力枚举，不用模拟退火）挑出"最弱那个拼接点的结合力最强"这个目标下最优的排列+spacer组合。
  真实MHCflurry模型、真实的"避免拼接处产生强结合表位"这个科学目标，只是没有照搬pVACvector自己的代码
  （它的模拟退火寻路 + 多算法取中位数的打分方式在只用一个算法（MHCflurry）时也用不上）。

## Pipeline Funnel

`main.py` 现在会输出每一步筛掉多少变异，而不只是最终三张表（数字是真实跑出来的，非固定）：

```
Protein-altering variants    25
Actionable variants           1   -> 进 drug_matches 分支
Neoantigen candidates        24   -> 进 pvactools 分支
Expressed variants           17   (TPM >= 1)
Real peptide generated        12  (missense only, 真实蛋白序列取到 + 位点对得上)
HLA-presented                 69  (IC50 <= 500nM，真实 mhcflurry 预测)
```

## Pathway 可视化

`pathway.py`（同一套设计思路来自 `luad_workflow/modules/06_pathway/kegg_viewer.py`）：自己在
本地缓存的 KEGG 官方 PNG 上用 Pillow 叠色块，不调 pathview/cytoscape。通路成员基因用 gseapy 的
KEGG_2021_Human gene set（本地缓存，不联网）；底图 + 基因框坐标是真的用 KEGG REST API/KGML 下载好
缓存在 `pipelines/downstream/kegg_cache/pathways/` 的。

覆盖范围：最早只挑了8条手选的"LUAD核心通路"，范围太窄——真实上传的VCF里突变基因很容易落在这8条之外。
现在改成 KEGG 自己的 BRITE 分类体系里 5 个跟肿瘤直接相关的官方类别（Signal transduction / Cancer:
overview / Cancer: specific types / Cell growth and death / Immune system），一共 **79 条通路**
（含 KEGG 自己的 Non-small cell lung cancer / Small cell lung cancer 通路图），不是随手挑的，是
KEGG 官方的分类。构建脚本是 `scripts/build_kegg_cache.py`（可重新跑，联网只发生在这一步，构建产物
79×(PNG + 坐标json) ≈ 9.4MB 全部提交进仓库，`pathway.py` 运行时只读本地文件，不联网）。只渲染实际
命中突变基因的通路。颜色（跟 `luad_workflow` 版本不同，没有差异表达倍数，只看是否表达/高表达）：

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

浏览器打开 Streamlit 给出的地址（默认 http://localhost:8501），默认打开就直接看到 TCGA-38-4627 的结果
（读的是预先算好的缓存，见下）。上传自己的 VCF/expression/HLA 再点「Run analysis」才会真的调后端跑一遍
（约 2 分钟，真实 MHC binding 预测 + vaccine construct 设计）。

## Demo 结果预计算（避免每次都重新跑 ~2 分钟的真实预测）

默认病例的结果不会变，没必要每次打开网页都重新跑一遍真实 pipeline。`scripts/precompute_demo.py` 把
`run_pipeline()` 的结果存成 `data/demo/precomputed_result.json`（~550KB，已提交进仓库），前端默认直接读
这个文件，秒开。改了 demo 数据或 pipeline 逻辑之后要记得重新生成：

```bash
python -m scripts.precompute_demo
```

## 部署（GitHub + Streamlit Community Cloud）

`frontend/streamlit_app.py` 本身只依赖 `streamlit`/`pandas`/`requests`（不 import 任何 pipeline 代码），
所以云端只需要 `frontend/requirements.txt` 这份轻量依赖，不需要装 `pvactools`/`tensorflow`/`mhcflurry`
这些重依赖 —— 云端版本只展示预先算好的 `precomputed_result.json`；上传自定义文件走真实分析这个功能，
需要本地起 FastAPI 后端（`uvicorn backend.main:app`）才能用，云端连不上后端会给出提示而不是报错崩溃。

部署到 Streamlit Community Cloud：仓库推到 GitHub 后，在 share.streamlit.io 里选这个仓库，
**Main file path 填 `frontend/streamlit_app.py`**（Cloud 会自动找同目录下的 `requirements.txt`）。

## 目录结构

```
data/demo/                        默认病例 TCGA-38-4627：真实 VCF + 真实 expression + synthetic HLA + 预计算结果
scripts/precompute_demo.py         重新生成 data/demo/precomputed_result.json
pipelines/downstream/              核心分析逻辑：vep.py / civic.py / pvactools.py / pathway.py / main.py（串联）
pipelines/downstream/kegg_cache/   KEGG 通路底图 PNG + 基因框坐标缓存（不联网）
backend/                           FastAPI，把 pipeline 包成 /analyze 接口（本地跑自定义上传时才需要）
frontend/                          Streamlit 网页；frontend/requirements.txt 是云端部署用的轻量依赖
```
