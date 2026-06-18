"""COSMIC Mutant file -> per-variant recurrence (distinct sample count), restricted to panel genes.

Design rationale (phase1_execution.md B-5, §6 circularity):
  - recurrence = **number of distinct COSMIC_SAMPLE_IDs** with the same protein change (missense).
  - ⚠️ ascertainment bias: MutantCensus includes targeted panel data → known hotspots are overrepresented.
    **GenomeScreensMutant** (genome/exome-wide, unbiased) is more honest as a recurrence FEATURE.
    → choose via --source. Default is genome screens.
  - ⚠️ circularity: using recurrence both as label (positive definition) and as feature is circular.
    Use this output only as a **feature candidate**, separate from labels (hotspot/CIViC). Ablation required.
  - Memory management: filter to CGC panel genes (768), then aggregate.

Output: data/interim/cosmic_recurrence__<source>.tsv  (per pkey)
Usage: python -m src.data.cosmic_recurrence --source genomescreens
"""
from __future__ import annotations
import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict

import pandas as pd

try:
    from .variant_keys import make_pkey, parse_protein_change, is_missense_change, protein_key
except ImportError:
    sys.path.insert(0, os.path.dirname(__file__))
    from variant_keys import make_pkey, parse_protein_change, is_missense_change, protein_key

SOURCES = {
    "mutantcensus": "Cosmic_MutantCensus_v104_GRCh38.tsv.gz",
    "genomescreens": "Cosmic_GenomeScreensMutant_v104_GRCh38.tsv.gz",
    "targeted": "Cosmic_CompleteTargetedScreensMutant_v104_GRCh38.tsv.gz",
}


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _ensure_extracted(cosmic_dir: str, fname: str) -> str:
    """If inside the tar, extract the .tsv.gz. If it already exists, keep it."""
    gz = os.path.join(cosmic_dir, fname)
    if os.path.exists(gz):
        return gz
    import tarfile
    tar = os.path.join(cosmic_dir, fname.replace("_v104_GRCh38.tsv.gz", "_Tsv_v104_GRCh38.tar"))
    if os.path.exists(tar):
        with tarfile.open(tar) as t:
            t.extract(fname, path=cosmic_dir)
        return gz
    raise SystemExit("COSMIC file not found: %s (or %s)" % (gz, tar))


def build_recurrence(gz_path: str, panel_genes: set[str], progress_every: int = 2_000_000):
    """⚠️ To avoid the transcript isoform numbering issue (e.g. GenomeScreens labels BRAF V600E as V640E),
    count distinct samples **per GRCh38 genome coordinate (chr,pos,ref,alt)**. Join with AM (GRCh38) by coordinate.
    Normalize chr by stripping the 'chr' prefix (AM='chr7' vs COSMIC='7')."""
    samples = defaultdict(set)     # gkey(chr,pos,ref,alt) -> set(sample_id)
    meta = {}                       # gkey -> (gene, representative protein_change)
    pchg_set = defaultdict(set)     # gkey -> observed protein_changes (consistency check)
    n = n_kept = 0
    with gzip.open(gz_path, "rt") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            n += 1
            if progress_every and n % progress_every == 0:
                print("  ...scanned %d rows, kept %d genomic variants" % (n, len(samples)), flush=True)
            gene = (row.get("GENE_SYMBOL") or "").strip().upper()
            if gene not in panel_genes:
                continue
            desc = row.get("MUTATION_DESCRIPTION", "") or ""
            if "missense" not in desc.lower():
                continue
            ch = (row.get("CHROMOSOME", "") or "").replace("chr", "")
            pos = row.get("GENOME_START", "") or ""
            wt = (row.get("GENOMIC_WT_ALLELE", "") or "").strip()
            mut = (row.get("GENOMIC_MUT_ALLELE", "") or "").strip()
            if not (ch and pos and wt and mut) or len(wt) != 1 or len(mut) != 1:
                continue  # SNV only (missense)
            gkey = (ch, pos, wt, mut)
            sid = row.get("COSMIC_SAMPLE_ID") or row.get("SAMPLE_NAME") or ""
            samples[gkey].add(sid)
            n_kept += 1
            if gkey not in meta:
                meta[gkey] = (gene, (row.get("MUTATION_AA", "") or "").replace("p.", ""))
            pchg_set[gkey].add((row.get("MUTATION_AA", "") or "").replace("p.", ""))
    rows = []
    for gkey, sset in samples.items():
        ch, ps, wt, mut = gkey
        gene, pchg = meta[gkey]
        rows.append({
            "grch38_chr": ch, "grch38_pos": int(ps), "grch38_ref": wt, "grch38_alt": mut,
            "gene": gene, "protein_change_repr": pchg,
            "n_protein_annotations": len(pchg_set[gkey]),  # >1 means transcript-specific numbering differs
            "cosmic_sample_count": len(sset),
        })
    df = pd.DataFrame(rows).sort_values("cosmic_sample_count", ascending=False).reset_index(drop=True)
    return df, n, n_kept


def main():
    ap = argparse.ArgumentParser(description="COSMIC recurrence (distinct sample count) per pkey, panel-restricted.")
    root = _repo_root()
    ap.add_argument("--source", choices=list(SOURCES), default="genomescreens")
    ap.add_argument("--cosmic-dir", default=os.path.join(root, "data", "external", "cosmic"))
    ap.add_argument("--panel", default=os.path.join(root, "data", "interim", "cancer_gene_panel.tsv"))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    panel = pd.read_csv(args.panel, sep="\t")
    panel_genes = set(panel["gene"].str.upper())
    print("[cosmic_recurrence] panel genes: %d | source: %s" % (len(panel_genes), args.source))

    gz = _ensure_extracted(args.cosmic_dir, SOURCES[args.source])
    df, n, n_kept = build_recurrence(gz, panel_genes)

    out = args.out or os.path.join(root, "data", "interim", "cosmic_recurrence__%s.tsv" % args.source)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print("[cosmic_recurrence] scanned %d rows, kept %d sample-hits" % (n, n_kept))
    print("  unique GRCh38 missense variants (panel): %d" % len(df))
    multi = int((df["n_protein_annotations"] > 1).sum())
    print("  single genomic variant with transcript-varying protein numbers(>1): %d (no effect since key is coordinate)" % multi)
    print("  top recurrence:")
    print(df.head(10)[["gene", "protein_change_repr", "grch38_chr", "grch38_pos",
                       "cosmic_sample_count"]].to_string(index=False))
    print("[cosmic_recurrence] wrote -> %s" % out)


if __name__ == "__main__":
    main()
