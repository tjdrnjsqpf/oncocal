"""COSMIC Mutant 파일 -> 변이별 recurrence (distinct sample 수), 패널 유전자 한정.

설계 근거 (phase1_execution.md B-5, §6 circularity):
  - recurrence = 같은 protein change(missense)를 가진 **distinct COSMIC_SAMPLE_ID 수**.
  - ⚠️ ascertainment bias: MutantCensus 는 타깃 패널 데이터 포함 → 알려진 hotspot 과대.
    **GenomeScreensMutant**(genome/exome-wide, unbiased)가 recurrence FEATURE 로 더 정직.
    → --source 로 선택. 기본은 genome screens.
  - ⚠️ circularity: recurrence 를 label(positive 정의)과 feature 양쪽에 쓰면 순환.
    이 산출물은 **feature 후보**로만 쓰고, label(hotspot/CIViC)과 분리. ablation 필수.
  - 메모리 관리: CGC 패널 유전자(768)로 필터 후 집계.

출력: data/interim/cosmic_recurrence__<source>.tsv  (pkey 단위)
사용: python -m src.data.cosmic_recurrence --source genomescreens
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
    """tar 안에 있으면 .tsv.gz 추출. 이미 있으면 그대로."""
    gz = os.path.join(cosmic_dir, fname)
    if os.path.exists(gz):
        return gz
    import tarfile
    tar = os.path.join(cosmic_dir, fname.replace("_v104_GRCh38.tsv.gz", "_Tsv_v104_GRCh38.tar"))
    if os.path.exists(tar):
        with tarfile.open(tar) as t:
            t.extract(fname, path=cosmic_dir)
        return gz
    raise SystemExit("COSMIC 파일 없음: %s (또는 %s)" % (gz, tar))


def build_recurrence(gz_path: str, panel_genes: set[str], progress_every: int = 2_000_000):
    """⚠️ transcript isoform 번호 문제(예: GenomeScreens가 BRAF V600E를 V640E로 표기) 회피를 위해
    **GRCh38 유전체 좌표(chr,pos,ref,alt) 단위**로 distinct sample 카운트. AM(GRCh38)과 좌표로 join.
    chr 은 'chr' 접두어 제거해 정규화(AM='chr7' vs COSMIC='7')."""
    samples = defaultdict(set)     # gkey(chr,pos,ref,alt) -> set(sample_id)
    meta = {}                       # gkey -> (gene, 대표 protein_change)
    pchg_set = defaultdict(set)     # gkey -> 관측된 protein_change 들(정합성 점검)
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
                continue  # SNV 만 (missense)
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
            "n_protein_annotations": len(pchg_set[gkey]),  # >1 이면 transcript별 번호 상이
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
    print("  유전체변이 1개가 transcript별로 다른 protein번호(>1): %d (좌표 키라 무영향)" % multi)
    print("  top recurrence:")
    print(df.head(10)[["gene", "protein_change_repr", "grch38_chr", "grch38_pos",
                       "cosmic_sample_count"]].to_string(index=False))
    print("[cosmic_recurrence] wrote -> %s" % out)


if __name__ == "__main__":
    main()
