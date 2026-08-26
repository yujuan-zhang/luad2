import shutil
import tempfile
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, UploadFile
from fastapi.middleware.cors import CORSMiddleware

from pipelines.downstream.main import run_pipeline

app = FastAPI(title="LUAD Precision Platform API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

DEMO_DIR = Path(__file__).resolve().parent.parent / "data" / "demo"


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/analyze")
async def analyze(
    vcf: Optional[UploadFile] = File(None),
    expression: Optional[UploadFile] = File(None),
    hla: Optional[UploadFile] = File(None),
):
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        vcf_path = await _save_or_default(vcf, tmp_path, "variants.vcf.gz", DEMO_DIR / "variants.vcf.gz")
        expression_path = await _save_or_default(expression, tmp_path, "expression.tsv.gz", DEMO_DIR / "expression.tsv.gz")
        hla_path = await _save_or_default(hla, tmp_path, "hla.tsv", DEMO_DIR / "hla.tsv")
        result = run_pipeline(vcf_path, expression_path, hla_path)
    return result


async def _save_or_default(upload, tmp_path, filename, default_path):
    if upload is None:
        return default_path
    dest = tmp_path / filename
    with dest.open("wb") as out:
        shutil.copyfileobj(upload.file, out)
    return dest
