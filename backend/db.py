import os

import psycopg2
from psycopg2.extras import Json

DATABASE_URL = os.environ.get(
    "DATABASE_URL", "postgresql://luadtx:luadtx_dev@localhost:5432/luadtx"
)


def get_connection():
    return psycopg2.connect(DATABASE_URL)


def init_db():
    with get_connection() as conn, conn.cursor() as cur:
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS cases (
                case_id TEXT PRIMARY KEY,
                cancer_type TEXT,
                clinical JSONB,
                note TEXT
            )
            """
        )
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS analysis_results (
                id SERIAL PRIMARY KEY,
                case_id TEXT,
                result JSONB,
                created_at TIMESTAMPTZ DEFAULT now()
            )
            """
        )


def save_case(case_id, cancer_type, clinical, note):
    with get_connection() as conn, conn.cursor() as cur:
        cur.execute(
            """
            INSERT INTO cases (case_id, cancer_type, clinical, note)
            VALUES (%s, %s, %s, %s)
            ON CONFLICT (case_id) DO UPDATE
            SET cancer_type = EXCLUDED.cancer_type,
                clinical = EXCLUDED.clinical,
                note = EXCLUDED.note
            """,
            (case_id, cancer_type, Json(clinical), note),
        )


def save_result(case_id, result):
    with get_connection() as conn, conn.cursor() as cur:
        cur.execute(
            "INSERT INTO analysis_results (case_id, result) VALUES (%s, %s)",
            (case_id, Json(result)),
        )
