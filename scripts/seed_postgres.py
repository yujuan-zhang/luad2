"""One-time (re-runnable) seed of Postgres with the demo case + its
precomputed result, so the DB isn't empty before the first real /analyze
call.

Usage:
    python -m scripts.seed_postgres
"""
import json
from pathlib import Path

from backend import db

DEMO_DIR = Path(__file__).resolve().parents[1] / "data" / "demo"


def main():
    db.init_db()

    case = json.loads((DEMO_DIR / "case_metadata.json").read_text())
    result = json.loads((DEMO_DIR / "precomputed_result.json").read_text())

    db.save_case(
        case_id=case["case_id"],
        cancer_type=case["cancer_type"],
        clinical=case["clinical"],
        note=case.get("note"),
    )
    db.save_result(case_id=case["case_id"], result=result)
    print(f"Seeded case {case['case_id']} + 1 analysis result into Postgres.")


if __name__ == "__main__":
    main()
