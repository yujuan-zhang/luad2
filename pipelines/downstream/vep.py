"""Mock VEP-style variant annotation.

MOCK: real VEP integration (Ensembl VEP REST API, or local VEP Docker for
batch use) will replace `_lookup_table` with an actual annotation call
later. The interface (`annotate_variants`) stays the same so nothing else
needs to change when that happens.
"""
import gzip

_lookup_table = {
    ("chr7", 55191822, "T", "G"): {
        "gene": "EGFR", "consequence": "missense_variant",
        "protein_change": "p.L858R", "transcript_id": "ENST00000275493",
    },
    ("chr12", 25245350, "C", "A"): {
        "gene": "KRAS", "consequence": "missense_variant",
        "protein_change": "p.G12C", "transcript_id": "ENST00000256078",
    },
    ("chr7", 140753336, "A", "T"): {
        "gene": "BRAF", "consequence": "missense_variant",
        "protein_change": "p.V600E", "transcript_id": "ENST00000288602",
    },
    ("chr17", 7674220, "C", "T"): {
        "gene": "TP53", "consequence": "missense_variant",
        "protein_change": "p.R175H", "transcript_id": "ENST00000269305",
    },
}


def _open(vcf_path):
    vcf_path = str(vcf_path)
    return gzip.open(vcf_path, "rt") if vcf_path.endswith(".gz") else open(vcf_path)


def annotate_variants(vcf_path):
    variants = []
    with _open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, *_ = line.rstrip("\n").split("\t")
            key = (chrom, int(pos), ref, alt)
            anno = _lookup_table.get(key, {
                "gene": "NA", "consequence": "unknown",
                "protein_change": "NA", "transcript_id": "NA",
            })
            variants.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, **anno})
    return variants
