"""Mock pVACtools-style neoantigen prediction.

MOCK: real integration will call pVACtools (mutant peptide generation +
MHC binding prediction) -- or, short of installing pVACtools itself, the
free Ensembl REST API (real mutant peptide sequence) + IEDB REST API (real
MHC binding prediction), both tokenless. The interface stays the same.
IC50 here is a deterministic hash, not a real binding prediction.

Everything upstream of this mock (which variants are candidates, whether
they're expressed, real DNA VAF/hotspot/functional-impact annotation) is
already real -- see vep.py.
"""
import hashlib

PROTEIN_ALTERING_CONSEQUENCES = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "splice_acceptor_variant", "splice_donor_variant",
}


def load_hla(hla_path):
    alleles = []
    with open(hla_path) as f:
        next(f)  # header
        for line in f:
            _locus, allele = line.rstrip("\n").split("\t")
            alleles.append(allele)
    return alleles


def _is_protein_altering(consequence):
    return any(c in PROTEIN_ALTERING_CONSEQUENCES for c in consequence.split(";"))


def filter_protein_altering(variants):
    return [v for v in variants if _is_protein_altering(v["consequence"])]


def filter_expressed(variants, expression_df, expression_threshold=1.0):
    """Keep variants whose gene clears the tumor-expression bar. Attaches
    the real TPM onto each variant for downstream display."""
    expr_lookup = dict(zip(expression_df["SYMBOL"], expression_df["TPM_GENE"]))
    kept = []
    for v in variants:
        tpm = expr_lookup.get(v["gene"], 0.0)
        if tpm >= expression_threshold:
            kept.append({**v, "tumor_tpm": tpm})
    return kept


def _mock_ic50(peptide, allele):
    # Range wide enough that the 500nM cutoff below actually filters some
    # pairs out, like a real binder/non-binder split would.
    digest = hashlib.md5(f"{peptide}{allele}".encode()).hexdigest()
    return 20 + (int(digest[:4], 16) % 1980)  # arbitrary 20-2000 nM range


def predict_neoantigens(expressed_variants, hla_alleles, ic50_threshold=500.0):
    """Rank candidate neoantigens among variants that already cleared the
    expression bar, by HLA presentation (IC50 below the standard
    weak-binder cutoff of 500 nM)."""
    results = []
    for v in expressed_variants:
        mutant_peptide = f"MOCKPEP-{v['protein_change']}"
        for allele in hla_alleles:
            ic50 = _mock_ic50(mutant_peptide, allele)
            if ic50 > ic50_threshold:
                continue
            results.append({
                "gene": v["gene"],
                "protein_change": v["protein_change"],
                "hla_allele": allele,
                "peptide": mutant_peptide,
                "ic50_nm": ic50,
                "tumor_tpm": v["tumor_tpm"],
                "dna_vaf": v.get("dna_vaf"),
                "hotspot": v.get("hotspot"),
                "functional_impact": v.get("functional_impact"),
            })
    results.sort(key=lambda r: r["ic50_nm"])
    return results
