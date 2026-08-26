"""Load pre-annotated somatic variants.

This does not call VEP itself -- per the project architecture, VEP
annotation happens upstream (nf-core/sarek -> VEP), so this stage expects
an already-annotated variant table matching that output schema (SYMBOL,
CHROM, POS, REF, ALT, CONSEQUENCE, HGVSp_Short, VAF_TUMOR, hotspot,
FUNCTIONAL_IMPACT, ...). The demo file is real VEP-annotated output for
TCGA-38-4627, not mock data.
"""
import csv
import gzip


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def annotate_variants(path):
    variants = []
    with _open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variants.append({
                "chrom": row["CHROM"],
                "pos": int(row["POS"]),
                "ref": row["REF"],
                "alt": row["ALT"],
                "gene": row["SYMBOL"],
                "consequence": row["CONSEQUENCE"],
                "protein_change": row.get("HGVSp_Short") or "NA",
                "dna_vaf": float(row["VAF_TUMOR"]) if row.get("VAF_TUMOR") else None,
                "hotspot": row.get("hotspot") == "Y",
                "functional_impact": row.get("FUNCTIONAL_IMPACT", ""),
            })
    return variants
