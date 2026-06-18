"""COSMIC Cancer Gene Census (CGC) -> 유전자 패널 + GOF/LOF gene-level proxy.

설계 근거 (phase1_execution.md B-1, gap_analysis per-gene 보정):
  - Scope = cancer gene panel 먼저 → CGC 768 유전자가 permissive 패널 출처(OncoKB 리스트 대신).
  - ROLE_IN_CANCER(oncogene/TSG/fusion) = 20/20 rule 의 직접적 gene-level proxy:
        oncogene -> GOF 경향,  TSG -> LOF 경향,  둘 다 -> ambiguous.
  - MOLECULAR_GENETICS(Dom/Rec) 보조: Dom=GOF/DN 경향, Rec=LOF 경향.
  - SYNONYMS 에 ENSG/UniProt 포함 → AlphaMissense(transcript/uniprot) 매핑에 사용.

⚠️ COSMIC 재배포 금지: 이 패널 TSV는 내부 interim 용. 논문엔 요약통계/유전자 수만.

출력: data/interim/cancer_gene_panel.tsv
사용: python -m src.data.cosmic_panel
"""
from __future__ import annotations
import argparse
import gzip
import os
import re

import pandas as pd

ENSG_RE = re.compile(r"\bENSG\d+(?:\.\d+)?\b")
# UniProt accession (단순화된 공식 정규식)
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
    """gene-level 기전 proxy (변이 단위 아님 — 약한 prior)."""
    if role_class == "oncogene":
        return "GOF"
    if role_class == "TSG":
        return "LOF"
    if role_class == "both":
        return "ambiguous"
    # role 불명 시 분자유전학으로 보조 추정
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
    print("\nENSG 매핑 가능: %d / UniProt 매핑 가능: %d"
          % (int((df["ensembl_gene"] != "").sum()), int((df["uniprot"] != "").sum())))
    print("[cosmic_panel] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
