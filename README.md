# LUAD Precision Platform

LUAD neoantigen / targeted-therapy prioritization pipeline. See project plan for full architecture.

## 当前阶段：Phase 1 — 本地全链路打通

`vep.py` / `oncokb.py` / `pvactools.py` 现在是 **mock 实现**（用固定的查表/伪随机逻辑模拟真实工具的输出），
目的是先把 `demo 数据 → downstream 分析 → FastAPI → Streamlit` 这条链路跑通、看到网页结果。
后续会把每个 mock 换成真实的 VEP / OncoKB / pVACtools 调用，接口（函数签名）保持不变。

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

## 目录结构（Phase 1）

```
data/demo/            内置的假 demo 数据（VCF、expression、HLA）
pipelines/downstream/ 核心分析逻辑：vep.py / oncokb.py / pvactools.py / main.py（串联）
backend/               FastAPI，把 pipeline 包成 /analyze 接口
frontend/               Streamlit 网页，调用 FastAPI 展示结果
```
