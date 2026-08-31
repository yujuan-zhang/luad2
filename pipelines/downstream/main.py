import time
from pathlib import Path

import pandas as pd

from . import civic, pathway, pvactools, vep


class _StageTimer:
    """Prints wall-clock time per pipeline stage -- these are network/model
    calls to external services (VEP, UniProt, CIViC, mhcflurry), so runtime
    is dominated by I/O latency, not CPU; this makes it visible which stage
    to optimize instead of guessing."""

    def __init__(self):
        self._t = time.perf_counter()

    def lap(self, label):
        now = time.perf_counter()
        print(f"[timing] {label}: {now - self._t:.1f}s")
        self._t = now


def run_pipeline(vcf_path, expression_path, hla_path):
    timer = _StageTimer()
    variants = vep.annotate_variants(vcf_path)
    timer.lap("vep.annotate_variants")
    protein_altering = pvactools.filter_protein_altering(variants)
    drug_matches = civic.match_drugs(protein_altering)
    timer.lap("civic.match_drugs")

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
    # generate_mutant_peptide only succeeds for missense variants with a
    # fetchable real sequence and a matching reference residue -- see
    # pvactools.py's module docstring for why nonsense/frameshift/splice
    # are out of scope. Counted separately here since predict_neoantigens
    # only returns the peptide-allele pairs that also clear the IC50 bar.
    with_peptide = pvactools.variants_with_peptide(expressed)
    neoantigens, peptide_hla_pairs = pvactools.predict_neoantigens(expressed, hla_alleles)
    timer.lap("pvactools.predict_neoantigens")
    vaccine_construct = pvactools.design_vaccine_construct(neoantigens, hla_alleles)
    timer.lap("pvactools.design_vaccine_construct")
    pathways = pathway.analyze_pathways(protein_altering, expression_df)
    timer.lap("pathway.analyze_pathways")

    # This is the neoantigen branch's flow specifically, so the actionable
    # variant that diverged to the targeted-therapy branch isn't a stage
    # here (it's surfaced in Key Clinical Finding instead) -- including it
    # would make the "funnel" visually grow partway through, which reads as
    # broken rather than as "one peptide can match several HLA alleles".
    # Short labels -- these are column headers in the UI. Threshold details
    # (TPM >= 1, IC50 <= 500nM) are spelled out in a caption instead of
    # crammed into the header.
    funnel = {
        "Variants": len(protein_altering),
        "Neoantigen Candidates": len(neoantigen_candidates),
        "Expressed": len(expressed),
        "Variant-derived Peptides": len(with_peptide),
        "Peptide-HLA Evaluations": peptide_hla_pairs,
        "Presented Candidates": len(neoantigens),
    }

    return {
        "funnel": funnel,
        "variants": variants,
        "drug_matches": drug_matches,
        "neoantigens": neoantigens,
        "vaccine_construct": vaccine_construct,
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
