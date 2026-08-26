from pathlib import Path

import pandas as pd

from . import oncokb, pvactools, vep


def run_pipeline(vcf_path, expression_path, hla_path):
    variants = vep.annotate_variants(vcf_path)
    drug_matches = oncokb.match_drugs(variants)
    expression_df = pd.read_csv(expression_path, sep="\t")
    hla_alleles = pvactools.load_hla(hla_path)
    neoantigens = pvactools.predict_neoantigens(variants, expression_df, hla_alleles)
    return {
        "variants": variants,
        "drug_matches": drug_matches,
        "neoantigens": neoantigens,
    }


if __name__ == "__main__":
    demo_dir = Path(__file__).resolve().parents[2] / "data" / "demo"
    result = run_pipeline(
        demo_dir / "variants.vcf.gz",
        demo_dir / "expression.tsv",
        demo_dir / "hla.tsv",
    )
    for section, rows in result.items():
        print(f"\n== {section} ==")
        for row in rows:
            print(row)
