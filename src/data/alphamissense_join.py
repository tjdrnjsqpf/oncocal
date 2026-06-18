"""AlphaMissense hg38 bulk -> 패널 유전자 missense 점수 (변이 universe + flagship 점수).

설계 근거 (SOURCES.md, dbNSFP 승인 지연 브리지):
  - dbNSFP(다도구) 도착 전까지 flagship germline 도구(AlphaMissense)로 score 백본 대체.
  - AM hg38: CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class.
  - 패널(CGC) 의 uniprot 으로 필터 → uniprot->gene 매핑으로 pkey 부여.
  - 이 산출물이 '패널 내 가능한 모든 missense'의 변이 universe + AM 점수 역할.
    (build_dataset 에서 이 universe 에 라벨/recurrence 를 left-join.)

⚠️ 단백질 번호 정합성: AM 은 UniProt canonical isoform 번호. CIViC/COSMIC 는 대표 transcript 번호.
   canonical 단백질에선 대개 일치(BRAF V600E 등). 불일치 유전자는 추후 MANE 정렬로 점검.
⚠️ uniprot 미보유 패널 유전자(~48)는 이 join 에서 빠짐 → 추후 ENSG/transcript 매핑 보강.

출력: data/interim/alphamissense_panel.tsv
사용: python -m src.data.alphamissense_join
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
    """패널 유전자 GRCh38 좌표구간 -> chr별 (start정렬) 배열. CGC UniProt 결손(BRAF 등 48개)
    회피를 위해 uniprot 대신 좌표로 gene 배정."""
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
    """(chr,pos) 를 포함하는 패널 유전자. 구간 비중첩 가정 + 인접 몇 개 확인."""
    arr = starts.get(ch)
    if not arr:
        return None
    i = bisect.bisect_right(arr, pos) - 1
    for j in range(i, max(-1, i - 3), -1):  # 약간의 중첩 대비 뒤로 몇 개
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
    # 동일 (좌표) 다중 transcript -> am_pathogenicity 최댓값 유지(좌표 기준 dedup)
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
    print("[am_join] panel genes (좌표구간 배정): %d" % len(panel))

    df, n = build_am_panel(args.am, panel)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print("[am_join] scanned %d AM rows" % n)
    print("  panel missense pkeys with AM score: %d (genes=%d)" % (len(df), df["gene"].nunique()))
    print("  am_class:"); print(df["am_class"].value_counts().to_string())
    print("[am_join] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
