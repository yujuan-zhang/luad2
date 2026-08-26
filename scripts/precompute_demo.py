#!/usr/bin/env python
"""Run the full pipeline once on the built-in TCGA-38-4627 case and cache
the result as JSON.

The default (no-upload) view in the Streamlit app loads this file directly
instead of calling the backend -- the demo case's answer doesn't change
between runs, so there's no reason to re-run the ~1 minute real MHC
binding prediction on every page load. Re-run this script whenever the
demo data or pipeline logic changes.

Usage: python scripts/precompute_demo.py
"""
import json
from pathlib import Path

from pipelines.downstream.main import run_pipeline

PROJECT_DIR = Path(__file__).resolve().parent.parent
DEMO_DIR = PROJECT_DIR / "data" / "demo"
OUT_PATH = DEMO_DIR / "precomputed_result.json"


def main():
    result = run_pipeline(
        DEMO_DIR / "variants.vcf.gz",
        DEMO_DIR / "expression.tsv.gz",
        DEMO_DIR / "hla.tsv",
    )
    OUT_PATH.write_text(json.dumps(result, indent=2))
    print(f"Wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
