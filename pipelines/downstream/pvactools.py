"""Mock pVACtools-style neoantigen prediction.

MOCK: real integration will call pVACtools (mutant peptide generation +
MHC binding prediction). The interface (`predict_neoantigens`) stays the
same. IC50 here is a deterministic hash, not a real binding prediction.
"""
import hashlib


def load_hla(hla_path):
    alleles = []
    with open(hla_path) as f:
        next(f)  # header
        for line in f:
            _locus, allele = line.rstrip("\n").split("\t")
            alleles.append(allele)
    return alleles


def _mock_ic50(peptide, allele):
    digest = hashlib.md5(f"{peptide}{allele}".encode()).hexdigest()
    return 20 + (int(digest[:4], 16) % 480)  # arbitrary 20-500 nM range


def predict_neoantigens(variants, expression_df, hla_alleles, expression_threshold=1.0):
    expr_lookup = dict(zip(expression_df["gene"], expression_df["tpm"]))
    results = []
    for v in variants:
        if v["consequence"] != "missense_variant":
            continue
        tpm = expr_lookup.get(v["gene"], 0.0)
        if tpm < expression_threshold:
            continue
        mutant_peptide = f"MOCKPEP-{v['protein_change']}"
        for allele in hla_alleles:
            results.append({
                "gene": v["gene"],
                "protein_change": v["protein_change"],
                "hla_allele": allele,
                "peptide": mutant_peptide,
                "ic50_nm": _mock_ic50(mutant_peptide, allele),
                "tumor_tpm": tpm,
            })
    results.sort(key=lambda r: r["ic50_nm"])
    return results
