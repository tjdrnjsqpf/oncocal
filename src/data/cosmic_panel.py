"""COSMIC Cancer Gene Census (CGC) -> gene panel + gene-level GOF/LOF proxy.

Design rationale (phase1_execution.md B-1, gap_analysis per-gene calibration):
  - Scope = cancer gene panel first → CGC 768 genes are the permissive panel source (instead of the OncoKB list).
  - ROLE_IN_CANCER (oncogene/TSG/fusion) = direct gene-level proxy for the 20/20 rule:
        oncogene -> GOF tendency,  TSG -> LOF tendency,  both -> ambiguous.
  - MOLECULAR_GENETICS (Dom/Rec) auxiliary: Dom=GOF/DN tendency, Rec=LOF tendency.
  - SYNONYMS includes ENSG/UniProt → used for AlphaMissense (transcript/uniprot) mapping.

⚠️ COSMIC redistribution prohibited: this panel TSV is for internal interim use. Only summary stats/gene counts in the paper.

Output: data/interim/cancer_gene_panel.tsv
Usage: python -m src.data.cosmic_panel
"""
from __future__ import annotations
import argparse
import gzip
import os
import re

import pandas as pd

ENSG_RE = re.compile(r"\bENSG\d+(?:\.\d+)?\b")
# UniProt accession (simplified official regex)
UNIPROT_RE = re.compile(r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}\b")


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _role_class(role: str) -> str:
    r = (role or "").lower()
    onc = "oncogene" in r
    tsg = "tsg" in r
    if onc and tsg:
        return "both"
    if onc:
        return "oncogene"
    if tsg:
        return "TSG"
    if "fusion" in r:
        return "fusion_only"
    return "unknown"


def _mech_proxy(role_class: str, mol_gen: str) -> str:
    """gene-level mechanism proxy (not per-variant — a weak prior)."""
    if role_class == "oncogene":
        return "GOF"
    if role_class == "TSG":
        return "LOF"
    if role_class == "both":
        return "ambiguous"
    # If role is unknown, estimate via molecular genetics as a fallback
    mg = (mol_gen or "").lower()
    if mg == "dom":
        return "GOF?"
    if mg == "rec":
        return "LOF?"
    return "unknown"


def build_panel(cgc_gz: str) -> pd.DataFrame:
    rows = []
    import csv
    with gzip.open(cgc_gz, "rt") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            syn = row.get("SYNONYMS", "") or ""
            ensg = ENSG_RE.search(syn)
            uni = UNIPROT_RE.search(syn)
            rc = _role_class(row.get("ROLE_IN_CANCER", ""))
            rows.append({
                "gene": (row.get("GENE_SYMBOL") or "").strip().upper(),
                "cosmic_gene_id": row.get("COSMIC_GENE_ID", ""),
                "chromosome": row.get("CHROMOSOME", ""),
                "genome_start": row.get("GENOME_START", ""),
                "genome_stop": row.get("GENOME_STOP", ""),
                "role_in_cancer": row.get("ROLE_IN_CANCER", ""),
                "role_class": rc,
                "molecular_genetics": row.get("MOLECULAR_GENETICS", ""),
                "mech_proxy": _mech_proxy(rc, row.get("MOLECULAR_GENETICS", "")),
                "tier": row.get("TIER", ""),
                "mutation_types": row.get("MUTATION_TYPES", ""),
                "ensembl_gene": ensg.group(0) if ensg else "",
                "uniprot": uni.group(0) if uni else "",
            })
    return pd.DataFrame(rows).drop_duplicates("gene").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description="COSMIC CGC -> cancer gene panel + GOF/LOF proxy.")
    root = _repo_root()
    ap.add_argument("--cgc", default=os.path.join(
        root, "data", "external", "cosmic", "Cosmic_CancerGeneCensus_v104_GRCh38.tsv.gz"))
    ap.add_argument("--out", default=os.path.join(root, "data", "interim", "cancer_gene_panel.tsv"))
    args = ap.parse_args()

    df = build_panel(args.cgc)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print("[cosmic_panel] genes: %d" % len(df))
    print("role_class:"); print(df["role_class"].value_counts().to_string())
    print("\nmech_proxy:"); print(df["mech_proxy"].value_counts().to_string())
    print("\nENSG mappable: %d / UniProt mappable: %d"
          % (int((df["ensembl_gene"] != "").sum()), int((df["uniprot"] != "").sum())))
    print("[cosmic_panel] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
