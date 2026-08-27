"""Neoantigen prediction: real mutant peptide + real MHC binding prediction.

Peptide sequence: fetched from UniProt's REST API (canonical Swiss-Prot
reviewed sequence for the gene; free, no token), substituted at the
mutated residue -- only for missense_variant, where a single amino-acid
substitution is the whole story. Nonsense/frameshift/splice variants would
need real CDS-level modeling (translate the actual shifted/read-through
sequence) to build a correct novel peptide; that's out of scope here, so
those are skipped rather than faked (Ensembl's REST API also has real
protein sequences, but its gene-symbol lookup endpoints were unreachable
from this environment while UniProt's were fast and reliable).

Binding prediction: `mhcflurry-predict`, the same real trained model
pVACtools' own MHCflurry wrapper (pvactools.lib.prediction_class.MHCflurry)
calls internally, and the one pVACtools' full pvacseq pipeline uses --
called directly here instead of through pvacseq's VEP-annotated-VCF
pipeline (needs a real VEP install with the Wildtype/Frameshift plugins,
not available here) or pVACtools' own wrapper class (which issues one
subprocess call, and therefore one from-scratch TensorFlow model load, per
allele per epitope length -- 10+ minutes for this case's variant count;
batching every peptide x allele into a single `mhcflurry-predict` call
brings that down to seconds). mhcflurry's models are pinned to an old
TF1-era Keras API (`tensorflow.compat.v1.keras...set_session`) that
current Keras 3 removed; TF_USE_LEGACY_KERAS=1 (via the tf-keras package)
restores it.

The reference residue at the mutated position is checked against the
fetched sequence before generating a peptide -- if it doesn't match (wrong
isoform, stale sequence, etc.) the variant is skipped rather than
producing a wrong peptide.
"""
import csv
import os
import re
import tempfile
from pathlib import Path
from subprocess import run

import requests

os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")

from pvactools.lib.run_utils import determine_neoepitopes

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
EPITOPE_LENGTHS = [8, 9, 10, 11]
FLANK = 10  # residues each side of the mutation -- covers every 8-11mer containing it

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


_protein_seq_cache = {}


def _fetch_protein_sequence(gene, timeout=15):
    """Canonical (Swiss-Prot reviewed) protein sequence for a gene symbol,
    via UniProt's REST API. Free, no token. Cached per gene per process."""
    if gene in _protein_seq_cache:
        return _protein_seq_cache[gene]
    seq = None
    try:
        resp = requests.get(
            UNIPROT_SEARCH_URL,
            params={
                "query": f"gene:{gene} AND organism_id:9606 AND reviewed:true",
                "fields": "accession,sequence",
                "format": "json",
            },
            timeout=timeout,
        )
        resp.raise_for_status()
        results = resp.json().get("results", [])
        if results:
            seq = results[0]["sequence"]["value"]
    except (requests.RequestException, KeyError, ValueError):
        seq = None
    _protein_seq_cache[gene] = seq
    return seq


def _parse_missense(protein_change):
    """'p.L858R' -> ('L', 858, 'R'). None for anything that isn't a plain
    single-residue substitution (stop_gained/frameshift/etc. HGVSp strings
    end in '*' or similar and won't match)."""
    m = re.match(r"^p\.([A-Z])(\d+)([A-Z])$", protein_change)
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)


def generate_mutant_peptide(gene, protein_change, flank=FLANK):
    """Real flanking peptide around a missense substitution, sliced from
    the gene's real canonical protein sequence. Returns (peptide,
    mutation_position) with mutation_position 1-based within peptide
    (matching determine_neoepitopes' window-start convention), or None if
    the sequence can't be fetched, the position is out of range, or the
    reference residue doesn't match VEP's call."""
    parsed = _parse_missense(protein_change)
    if parsed is None:
        return None
    wt_aa, pos, mut_aa = parsed

    seq = _fetch_protein_sequence(gene)
    if not seq or pos > len(seq) or seq[pos - 1] != wt_aa:
        return None

    start = max(0, pos - 1 - flank)
    end = min(len(seq), pos + flank)
    mutant_peptide = seq[start:pos - 1] + mut_aa + seq[pos:end]
    mutation_position = (pos - 1 - start) + 1
    return mutant_peptide, mutation_position


def variants_with_peptide(expressed_variants):
    """Subset of variants a real mutant peptide could be generated for --
    missense only, with a fetchable sequence and a matching reference
    residue. Used for funnel reporting."""
    return [v for v in expressed_variants if generate_mutant_peptide(v["gene"], v["protein_change"]) is not None]


def _candidate_windows(peptide, mutation_position, lengths=EPITOPE_LENGTHS):
    """All epitope windows (of the given lengths) that cover the mutated
    residue -- a window identical to the wildtype elsewhere in the flank
    isn't a neoantigen."""
    windows = []
    for length in lengths:
        for start, epitope in determine_neoepitopes(peptide, length).items():
            if start <= mutation_position < start + length:
                windows.append(epitope)
    return windows


def _run_mhcflurry(peptides, alleles, tmp_dir):
    """One batched `mhcflurry-predict` call over every peptide x allele
    pair. Returns {(allele, peptide): (ic50, percentile)}."""
    out_path = Path(tmp_dir) / "mhcflurry_out.csv"
    result = run(
        [
            "mhcflurry-predict",
            "--alleles", *alleles,
            "--peptides", *peptides,
            "--out", str(out_path),
        ],
        capture_output=True, text=True,
    )
    if result.returncode != 0 or not out_path.exists():
        return {}

    lookup = {}
    with out_path.open() as f:
        for row in csv.DictReader(f):
            ic50 = row.get("mhcflurry_affinity") or row.get("mhcflurry_prediction")
            percentile = row.get("mhcflurry_affinity_percentile") or row.get("mhcflurry_prediction_percentile")
            if ic50 is None:
                continue
            lookup[(row["allele"], row["peptide"])] = (
                float(ic50),
                float(percentile) if percentile not in (None, "") else None,
            )
    return lookup


def predict_neoantigens(expressed_variants, hla_alleles, ic50_threshold=500.0, tmp_dir=None):
    """Real neoantigen candidates: real mutant peptides, real MHC-I binding
    prediction (mhcflurry), kept only where the predicted epitope window
    actually covers the mutated residue and clears the HLA-presentation
    IC50 cutoff.

    Returns (results, evaluated_pairs) -- evaluated_pairs is the real count
    of (variant, allele, candidate-epitope-window) triples submitted to
    mhcflurry, for funnel reporting.
    """
    candidates = []  # (variant, mutation_position, [candidate epitope windows])
    all_peptides = set()
    for v in expressed_variants:
        generated = generate_mutant_peptide(v["gene"], v["protein_change"])
        if generated is None:
            continue
        flank_peptide, mutation_position = generated
        windows = _candidate_windows(flank_peptide, mutation_position)
        if not windows:
            continue
        candidates.append((v, windows))
        all_peptides.update(windows)

    evaluated_pairs = sum(len(windows) for _, windows in candidates) * len(hla_alleles)

    if not candidates or not hla_alleles or not all_peptides:
        return [], evaluated_pairs

    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        binding = _run_mhcflurry(sorted(all_peptides), hla_alleles, work_dir)

    results = []
    for v, windows in candidates:
        for allele in hla_alleles:
            for epitope in windows:
                hit = binding.get((allele, epitope))
                if hit is None:
                    continue
                ic50, percentile = hit
                if ic50 > ic50_threshold:
                    continue
                results.append({
                    "gene": v["gene"],
                    "protein_change": v["protein_change"],
                    "hla_allele": allele,
                    "peptide": epitope,
                    "ic50_nm": round(ic50, 1),
                    "percentile": round(percentile, 3) if percentile is not None else None,
                    "tumor_tpm": v["tumor_tpm"],
                    "dna_vaf": v.get("dna_vaf"),
                    "hotspot": v.get("hotspot"),
                    "functional_impact": v.get("functional_impact"),
                })
    results.sort(key=lambda r: r["ic50_nm"])
    return results, evaluated_pairs
