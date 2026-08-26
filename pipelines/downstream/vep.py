"""Somatic variant annotation.

MOCK mechanism, real content: annotates a somatic VCF by looking up each
(chrom, pos, ref, alt) in `_ANNOTATION_TABLE`, which holds the real VEP
annotation results for TCGA-38-4627 (from `luad_workflow`'s real pipeline)
-- gene, consequence, HGVSp, DNA VAF, hotspot flag, functional impact.

Matches the intended architecture: users upload a VCF, and (once wired
up) an AWS-hosted VEP API annotates it. `annotate_variants`'s interface
won't change when the lookup table is replaced with that real call.
"""
import gzip

_ANNOTATION_TABLE = {
    ("chr1", 207899997, "C", "T"): {"gene": "CD34", "consequence": "missense_variant", "protein_change": "p.G29E", "dna_vaf": 0.1429, "hotspot": False, "functional_impact": "VUS"},
    ("chr2", 189877340, "C", "G"): {"gene": "PMS1", "consequence": "missense_variant", "protein_change": "p.I901M", "dna_vaf": 0.0889, "hotspot": False, "functional_impact": "VUS"},
    ("chr3", 141577093, "G", "A"): {"gene": "RASA2", "consequence": "missense_variant", "protein_change": "p.R526Q", "dna_vaf": 0.1102, "hotspot": False, "functional_impact": "Possibly_Functional"},
    ("chr3", 143832090, "C", "A"): {"gene": "SLC9A9", "consequence": "stop_gained", "protein_change": "p.E103*", "dna_vaf": 0.1, "hotspot": False, "functional_impact": "Likely_Functional"},
    ("chr5", 138192270, "G", "A"): {"gene": "CDC23", "consequence": "stop_gained;splice_region_variant", "protein_change": "p.R429*", "dna_vaf": 0.1471, "hotspot": False, "functional_impact": "Likely_Functional"},
    ("chr5", 141854106, "G", "A"): {"gene": "PCDH1", "consequence": "missense_variant", "protein_change": "p.S1217L", "dna_vaf": 0.25, "hotspot": False, "functional_impact": "VUS"},
    ("chr6", 20758610, "C", "T"): {"gene": "CDKAL1", "consequence": "missense_variant", "protein_change": "p.R162C", "dna_vaf": 0.1111, "hotspot": False, "functional_impact": "VUS"},
    ("chr6", 25813195, "A", "G"): {"gene": "SLC17A1", "consequence": "missense_variant", "protein_change": "p.V212A", "dna_vaf": 0.1167, "hotspot": False, "functional_impact": "VUS"},
    ("chr6", 33412776, "G", "A"): {"gene": "PHF1", "consequence": "missense_variant", "protein_change": "p.C107Y", "dna_vaf": 0.0667, "hotspot": False, "functional_impact": "Possibly_Functional"},
    ("chr7", 12230396, "A", "G"): {"gene": "TMEM106B", "consequence": "missense_variant", "protein_change": "p.Y197C", "dna_vaf": 0.053, "hotspot": False, "functional_impact": "Possibly_Functional"},
    ("chr7", 55142382, "T", "G"): {"gene": "EGFR", "consequence": "missense_variant", "protein_change": "p.L62R", "dna_vaf": 0.15, "hotspot": True, "functional_impact": "VUS"},
    ("chr7", 55191822, "T", "G"): {"gene": "EGFR", "consequence": "missense_variant", "protein_change": "p.L858R", "dna_vaf": 0.1343, "hotspot": True, "functional_impact": "Possibly_Functional"},
    ("chr7", 149103739, "G", "A"): {"gene": "ZNF425", "consequence": "missense_variant", "protein_change": "p.A711V", "dna_vaf": 0.186, "hotspot": False, "functional_impact": "VUS"},
    ("chr8", 87873503, "C", "T"): {"gene": "DCAF4L2", "consequence": "missense_variant", "protein_change": "p.V157M", "dna_vaf": 0.0645, "hotspot": False, "functional_impact": "VUS"},
    ("chr12", 85283653, "G", "A"): {"gene": "ALX1", "consequence": "missense_variant", "protein_change": "p.R103Q", "dna_vaf": 0.12, "hotspot": False, "functional_impact": "VUS"},
    ("chr13", 69740501, "G", "T"): {"gene": "KLHL1", "consequence": "missense_variant", "protein_change": "p.S565R", "dna_vaf": 0.08, "hotspot": False, "functional_impact": "Possibly_Functional"},
    ("chr14", 20117962, "C", "T"): {"gene": "OR4K17", "consequence": "missense_variant", "protein_change": "p.H155Y", "dna_vaf": 0.1067, "hotspot": False, "functional_impact": "Possibly_Functional"},
    ("chr16", 53292902, "G", "A"): {"gene": "CHD9", "consequence": "missense_variant", "protein_change": "p.R1787H", "dna_vaf": 0.1484, "hotspot": False, "functional_impact": "Possibly_Functional"},
    ("chr17", 3734821, "T", "A"): {"gene": "ITGAE", "consequence": "missense_variant", "protein_change": "p.Q884L", "dna_vaf": 0.098, "hotspot": False, "functional_impact": "VUS"},
    ("chr17", 7513369, "C", "T"): {"gene": "POLR2A", "consequence": "missense_variant", "protein_change": "p.T1702I", "dna_vaf": 0.1786, "hotspot": False, "functional_impact": "VUS"},
    ("chr17", 10639422, "A", "-"): {"gene": "MYH3", "consequence": "frameshift_variant", "protein_change": "p.L993*", "dna_vaf": 0.1064, "hotspot": False, "functional_impact": "Likely_Functional"},
    ("chr17", 60947218, "G", "T"): {"gene": "BCAS3", "consequence": "splice_acceptor_variant", "protein_change": "p.X363_splice", "dna_vaf": 0.0903, "hotspot": False, "functional_impact": "Likely_Functional"},
    ("chr22", 19410772, "G", "A"): {"gene": "HIRA", "consequence": "missense_variant", "protein_change": "p.P15L", "dna_vaf": 0.1051, "hotspot": False, "functional_impact": "Possibly_Functional"},
    ("chr22", 32437824, "A", "T"): {"gene": "BPIFC", "consequence": "missense_variant", "protein_change": "p.L228Q", "dna_vaf": 0.0957, "hotspot": False, "functional_impact": "VUS"},
    ("chrX", 31478264, "T", "A"): {"gene": "DMD", "consequence": "stop_gained", "protein_change": "p.R2927*", "dna_vaf": 0.0704, "hotspot": False, "functional_impact": "Likely_Functional"},
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
            anno = _ANNOTATION_TABLE.get(key, {
                "gene": "NA", "consequence": "unknown",
                "protein_change": "NA", "dna_vaf": None,
                "hotspot": False, "functional_impact": "",
            })
            variants.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, **anno})
    return variants
