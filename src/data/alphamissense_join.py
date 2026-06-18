"""AlphaMissense hg38 bulk -> panel-gene missense scores (variant universe + flagship scores).

Design rationale (SOURCES.md, bridge during the dbNSFP approval delay):
  - Until dbNSFP (multi-tool) arrives, substitute the flagship germline tool (AlphaMissense) as the score backbone.
  - AM hg38: CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class.
  - Filter by the panel (CGC) uniprot → assign pkey via uniprot->gene mapping.
  - This output serves as the variant universe of 'all possible missense in the panel' + AM scores.
    (build_dataset left-joins labels/recurrence onto this universe.)

⚠️ Protein numbering consistency: AM uses UniProt canonical isoform numbering. CIViC/COSMIC use representative transcript numbering.
   For canonical proteins these usually agree (BRAF V600E etc.). Mismatched genes to be checked later via MANE alignment.
⚠️ Panel genes without uniprot (~48) drop out of this join → supplement later with ENSG/transcript mapping.

Output: data/interim/alphamissense_panel.tsv
Usage: python -m src.data.alphamissense_join
"""
from __future__ import annotations
import argparse
import csv
import gzip
import os
import sys

import pandas as pd

try:
    from .variant_keys import parse_protein_change, is_missense_change, protein_key
except ImportError:
    sys.path.insert(0, os.path.dirname(__file__))
    from variant_keys import parse_protein_change, is_missense_change, protein_key


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _build_intervals(panel: pd.DataFrame):
    """Panel-gene GRCh38 coordinate intervals -> per-chr (start-sorted) arrays. To avoid the CGC UniProt
    gaps (BRAF etc., 48 genes), assign gene by coordinate instead of uniprot."""
    import bisect
    by_chr = {}
    for _, r in panel.iterrows():
        ch = str(r["chromosome"]).replace("chr", "").strip()
        try:
            s, e = int(r["genome_start"]), int(r["genome_stop"])
        except (ValueError, TypeError):
            continue
        by_chr.setdefault(ch, []).append((s, e, str(r["gene"]).upper()))
    starts, ivs = {}, {}
    for ch, lst in by_chr.items():
        lst.sort()
        ivs[ch] = lst
        starts[ch] = [x[0] for x in lst]
    return starts, ivs, bisect


def _gene_at(ch, pos, starts, ivs, bisect):
    """Panel gene containing (chr,pos). Assumes non-overlapping intervals + checks a few neighbors."""
    arr = starts.get(ch)
    if not arr:
        return None
    i = bisect.bisect_right(arr, pos) - 1
    for j in range(i, max(-1, i - 3), -1):  # a few steps back in case of slight overlap
        if j < 0:
            break
        s, e, g = ivs[ch][j]
        if s <= pos <= e:
            return g
    return None


def build_am_panel(am_gz: str, panel: pd.DataFrame, progress_every: int = 10_000_000):
    starts, ivs, bisect = _build_intervals(panel)
    rows = []
    n = 0
    with gzip.open(am_gz, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            n += 1
            if progress_every and n % progress_every == 0:
                print("  ...scanned %d AM rows, kept %d" % (n, len(rows)), flush=True)
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            ch = p[0].replace("chr", "")
            try:
                gpos = int(p[1])
            except ValueError:
                continue
            gene = _gene_at(ch, gpos, starts, ivs, bisect)
            if gene is None:
                continue
            pc = parse_protein_change(p[7])  # protein_variant e.g. V600E
            if pc is None:
                continue
            ref, pos, alt = pc
            if not is_missense_change(ref, alt):
                continue
            rows.append((protein_key(gene, ref, pos, alt), gene, "%s%d%s" % (ref, pos, alt),
                         ch, gpos, p[2], p[3], p[5], p[6], float(p[8]), p[9]))
    df = pd.DataFrame(rows, columns=[
        "pkey", "gene", "protein_change", "grch38_chr", "grch38_pos",
        "grch38_ref", "grch38_alt", "uniprot", "transcript", "am_pathogenicity", "am_class"])
    # Same (coordinate) across multiple transcripts -> keep the max am_pathogenicity (dedup by coordinate)
    df = (df.sort_values("am_pathogenicity", ascending=False)
            .drop_duplicates(["grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt"])
            .reset_index(drop=True))
    return df, n


def main():
    ap = argparse.ArgumentParser(description="AlphaMissense hg38 -> panel missense scores.")
    root = _repo_root()
    ap.add_argument("--am", default=os.path.join(
        root, "data", "external", "alphamissense", "AlphaMissense_hg38.tsv.gz"))
    ap.add_argument("--panel", default=os.path.join(root, "data", "interim", "cancer_gene_panel.tsv"))
    ap.add_argument("--out", default=os.path.join(root, "data", "interim", "alphamissense_panel.tsv"))
    args = ap.parse_args()

    panel = pd.read_csv(args.panel, sep="\t")
    print("[am_join] panel genes (assigned by coordinate interval): %d" % len(panel))

    df, n = build_am_panel(args.am, panel)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print("[am_join] scanned %d AM rows" % n)
    print("  panel missense pkeys with AM score: %d (genes=%d)" % (len(df), df["gene"].nunique()))
    print("  am_class:"); print(df["am_class"].value_counts().to_string())
    print("[am_join] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
