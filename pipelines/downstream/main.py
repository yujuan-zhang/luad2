from pathlib import Path

import pandas as pd

from . import civic, pathway, pvactools, vep


def run_pipeline(vcf_path, expression_path, hla_path):
    variants = vep.annotate_variants(vcf_path)
    protein_altering = pvactools.filter_protein_altering(variants)
    drug_matches = civic.match_drugs(protein_altering)

    # Actionable mutations go to the targeted-therapy branch; everything
    # else that clears the expression + HLA-presentation bar in pvactools
    # goes to neoantigen ranking instead. Mutually exclusive by design.
    actionable = {
        (m["gene"], m["protein_change"])
        for m in drug_matches if civic.is_actionable(m)
    }
    neoantigen_candidates = [
        v for v in protein_altering if (v["gene"], v["protein_change"]) not in actionable
    ]

    expression_df = pd.read_csv(expression_path, sep="\t")
    hla_alleles = pvactools.load_hla(hla_path)
    expressed = pvactools.filter_expressed(neoantigen_candidates, expression_df)
    neoantigens = pvactools.predict_neoantigens(expressed, hla_alleles)
    pathways = pathway.analyze_pathways(protein_altering, expression_df)

    # Short labels -- these are column headers in the UI's funnel table.
    # Threshold details (TPM >= 1, IC50 <= 500nM) are spelled out in a
    # caption next to the table instead of crammed into the header.
    funnel = {
        "Protein-altering variants": len(protein_altering),
        "Actionable variants": len(actionable),
        "Neoantigen candidates": len(neoantigen_candidates),
        "Expressed variants": len(expressed),
        "Peptide-HLA pairs": len(expressed) * len(hla_alleles),
        "HLA-presented": len(neoantigens),
    }

    return {
        "funnel": funnel,
        "variants": variants,
        "drug_matches": drug_matches,
        "neoantigens": neoantigens,
        "pathways": pathways,
    }


if __name__ == "__main__":
    demo_dir = Path(__file__).resolve().parents[2] / "data" / "demo"
    result = run_pipeline(
        demo_dir / "variants.vcf.gz",
        demo_dir / "expression.tsv.gz",
        demo_dir / "hla.tsv",
    )
    for section, rows in result.items():
        print(f"\n== {section} ==")
        if isinstance(rows, dict):
            for k, v in rows.items():
                print(f"  {k}: {v}")
        else:
            for row in rows:
                row = {k: (f"<{len(v)} b64 chars>" if k == "image_png_base64" and v else v)
                       for k, v in row.items()} if isinstance(row, dict) else row
                print(row)
