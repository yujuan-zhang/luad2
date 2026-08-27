"""Somatic variant annotation via the real Ensembl VEP REST API.

Free, no token, GRCh38 (matches this project's coordinates). Variants are
converted from VCF notation to Ensembl's region/allele format and POSTed in
batches (the API accepts up to `_BATCH_SIZE` per call) to
`/vep/human/region`, with `canonical=1` so annotation is pinned to each
gene's canonical transcript. `gene`, `consequence` and `impact` come
straight from the response; `protein_change` is built from `amino_acids` +
`protein_start` (e.g. "L/R" at 858 -> "p.L858R") rather than requesting
HGVS strings, since Ensembl's `hgvs=1` option 500-errors on this endpoint
as of 2026-08.

`functional_impact` is VEP's own IMPACT tier (HIGH/MODERATE/LOW/MODIFIER)
-- not the PCGR/cancer-specific tiering an earlier mock version of this
module used, which isn't a real VEP output. `hotspot` is *not* derived from
the API response: COSMIC co-location (`colocated_variants[].somatic`) flags
almost any observed somatic variant, not just recurrent drivers, so it's
checked against `_KNOWN_LUAD_HOTSPOTS`, a small curated list of
well-established LUAD driver hotspots -- same pattern as `civic.DRUG_KB`.

DNA VAF isn't a VEP concept at all (VEP annotates genomic consequence, not
genotype/frequency); it's read directly from the input VCF's own `VAF=` INFO
field, which is where a real variant caller would put it.

No annotation is fabricated: a variant with no INFO/VAF gets `dna_vaf=None`,
and a failed API call raises rather than silently returning empty/wrong
annotations.
"""
import gzip

import requests

ENSEMBL_VEP_URL = "https://rest.ensembl.org/vep/human/region"
_BATCH_SIZE = 200

# Well-established LUAD driver hotspots (gene -> set of HGVSp short forms).
# Deliberately small and conservative -- textbook driver mutations only, not
# a full cancerhotspots.org-style statistical hotspot table.
_KNOWN_LUAD_HOTSPOTS = {
    "EGFR": {"p.L858R", "p.T790M", "p.G719A", "p.G719S", "p.G719C", "p.L861Q", "p.S768I"},
    "KRAS": {"p.G12C", "p.G12D", "p.G12V", "p.G12A", "p.G12S", "p.G13D", "p.Q61H", "p.Q61K", "p.Q61L"},
    "BRAF": {"p.V600E"},
    "ERBB2": {"p.L755S"},
}


def _is_hotspot(gene, protein_change):
    return protein_change in _KNOWN_LUAD_HOTSPOTS.get(gene, set())


def _open(vcf_path):
    vcf_path = str(vcf_path)
    return gzip.open(vcf_path, "rt") if vcf_path.endswith(".gz") else open(vcf_path)


def _parse_info_vaf(info):
    for field in info.split(";"):
        if field.startswith("VAF="):
            try:
                return float(field[len("VAF="):])
            except ValueError:
                return None
    return None


def _read_vcf(vcf_path):
    """[(chrom, pos, ref, alt, dna_vaf), ...] in file order."""
    records = []
    with _open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, _qual, _filt, info, *_ = line.rstrip("\n").split("\t")
            records.append((chrom, int(pos), ref, alt, _parse_info_vaf(info)))
    return records


def _vcf_to_region(chrom, pos, ref, alt):
    """VCF (chrom, 1-based pos, ref, alt) -> Ensembl region-endpoint
    (chrom, start, end, "ref/alt") triple. Handles plain SNVs, this
    project's MAF-style unanchored indels ("-" for the missing allele),
    and standard VCF-anchored indels (shared-prefix trimming)."""
    ref = "" if ref == "-" else ref
    alt = "" if alt == "-" else alt
    chrom = chrom.replace("chr", "")

    if len(ref) == 1 and len(alt) == 1 and ref and alt:
        return chrom, pos, pos, f"{ref}/{alt}"
    if ref == "":
        return chrom, pos, pos - 1, f"-/{alt}"
    if alt == "":
        return chrom, pos, pos + len(ref) - 1, f"{ref}/-"

    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    trimmed_ref, trimmed_alt = ref[i:], alt[i:]
    start = pos + i
    if trimmed_ref == "":
        return chrom, start, start - 1, f"-/{trimmed_alt}"
    if trimmed_alt == "":
        return chrom, start, start + len(trimmed_ref) - 1, f"{trimmed_ref}/-"
    return chrom, start, start + len(trimmed_ref) - 1, f"{trimmed_ref}/{trimmed_alt}"


def _query_vep(region_strings, timeout=60):
    """One batched POST to the Ensembl VEP REST API. Raises on failure --
    annotation is core pipeline output, not a best-effort enrichment, so a
    broken call should surface loudly rather than annotate silently wrong."""
    resp = requests.post(
        ENSEMBL_VEP_URL,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        json={"variants": region_strings, "canonical": 1},
        timeout=timeout,
    )
    resp.raise_for_status()
    results = resp.json()
    if len(results) != len(region_strings):
        raise RuntimeError(
            f"Ensembl VEP returned {len(results)} results for {len(region_strings)} "
            "variants -- can't safely match annotations back to input order."
        )
    return results


def _pick_transcript(vep_result):
    consequences = vep_result.get("transcript_consequences", [])
    canonical = [t for t in consequences if t.get("canonical") == 1]
    return (canonical or consequences or [{}])[0]


def _protein_change(transcript):
    amino_acids = transcript.get("amino_acids")
    pos = transcript.get("protein_start")
    if not amino_acids or "/" not in amino_acids or not pos:
        return "NA"
    wt, mut = amino_acids.split("/")
    return f"p.{wt}{pos}{mut}"


def annotate_variants(vcf_path):
    records = _read_vcf(vcf_path)
    if not records:
        return []

    variants = []
    for batch_start in range(0, len(records), _BATCH_SIZE):
        batch = records[batch_start:batch_start + _BATCH_SIZE]
        region_strings = [
            "{} {} {} {} 1".format(*_vcf_to_region(chrom, pos, ref, alt))
            for chrom, pos, ref, alt, _vaf in batch
        ]
        results = _query_vep(region_strings)

        for (chrom, pos, ref, alt, dna_vaf), vep_result in zip(batch, results):
            transcript = _pick_transcript(vep_result)
            gene = transcript.get("gene_symbol", "NA")
            consequence = ";".join(
                transcript.get("consequence_terms")
                or [vep_result.get("most_severe_consequence", "unknown")]
            )
            protein_change = _protein_change(transcript)
            variants.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "gene": gene,
                "consequence": consequence,
                "protein_change": protein_change,
                "dna_vaf": dna_vaf,
                "hotspot": _is_hotspot(gene, protein_change),
                "functional_impact": transcript.get("impact", ""),
            })
    return variants
