"""cancerhotspots.org v2 -> 잔기 hotspot 기반 positive 보조 라벨.

설계 근거 (phase1_execution.md C, civic-label 한계):
  - CIViC 'Oncogenic' 근거는 희소(미세변이 positive ~46) + canonical driver 미플래그.
    cancerhotspots 의 통계적으로 유의한 recurrent hotspot residue 가 positive 의 주 기둥.
  - SNV-hotspots 시트: 1행 = (residue, 특정 alt AA). 컬럼:
      Reference_Amino_Acid='Q:422'(잔기 ref:잔기총카운트), Variant_Amino_Acid='R:204'(alt:해당 alt 카운트),
      Mutation_Count=잔기 총 변이 수, qvalue/log10_pvalue=hotspot 유의성.
  - ⚠️ Amino_Acid_Position 은 문자열 dtype → int 변환 필요.
  - ⚠️ 함정: hotspot recurrence 는 feature/label 양쪽에 들어갈 수 있어 순환 위험.
    여기서는 '잔기 단위 hotspot 라벨'만 만들고, recurrence 카운트를 label에 직접 넣지 않는다
    (카운트는 분석/보조용으로만 보존). build_dataset 에서 ablation 으로 점검.

출력: data/interim/cancerhotspots_labels.tsv (pkey 단위, alt-specific missense)
사용: python -m src.data.cancerhotspots_labels
"""
from __future__ import annotations
import argparse
import os
import sys

import pandas as pd

try:
    from .variant_keys import AA1, protein_key
except ImportError:
    sys.path.insert(0, os.path.dirname(__file__))
    from variant_keys import AA1, protein_key


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _split_aa_count(s):
    """'Q:422' -> ('Q', 422). 비정상 입력은 (None, None)."""
    if not isinstance(s, str) or ":" not in s:
        return None, None
    aa, _, cnt = s.partition(":")
    aa = aa.strip()
    try:
        return aa, int(cnt)
    except ValueError:
        return aa, None


def build_hotspot_labels(xls_path: str, qvalue_max: float = 0.05) -> pd.DataFrame:
    df = pd.read_excel(xls_path, sheet_name="SNV-hotspots")
    rows = []
    for _, r in df.iterrows():
        gene = str(r["Hugo_Symbol"]).strip()
        try:
            pos = int(str(r["Amino_Acid_Position"]).strip())
        except (ValueError, TypeError):
            continue
        ref_aa, residue_count = _split_aa_count(r["Reference_Amino_Acid"])
        alt_aa, alt_count = _split_aa_count(r["Variant_Amino_Acid"])
        # 단일 missense 만: ref/alt 표준 20 AA, alt != ref, nonsense(*) 제외
        if not ref_aa or not alt_aa:
            continue
        if ref_aa not in AA1 or alt_aa not in AA1 or alt_aa == ref_aa:
            continue
        qv = pd.to_numeric(r.get("qvalue"), errors="coerce")
        rows.append({
            "pkey": protein_key(gene, ref_aa, pos, alt_aa),
            "gene": gene, "ref_aa": ref_aa, "pos": pos, "alt_aa": alt_aa,
            "protein_change": "%s%d%s" % (ref_aa, pos, alt_aa),
            "hotspot_residue_count": residue_count,   # 잔기 총 변이 수(보조, label 직접사용 금지)
            "hotspot_alt_count": alt_count,            # 해당 alt 변이 수
            "hotspot_qvalue": qv,
            "hotspot_log10p": pd.to_numeric(r.get("log10_pvalue"), errors="coerce"),
            "is_hotspot": bool(pd.notna(qv) and qv <= qvalue_max),
        })
    out = pd.DataFrame(rows).drop_duplicates("pkey").reset_index(drop=True)
    out = out.sort_values("hotspot_alt_count", ascending=False).reset_index(drop=True)
    return out


def main():
    ap = argparse.ArgumentParser(description="cancerhotspots v2 -> residue hotspot positive labels.")
    root = _repo_root()
    ap.add_argument("--xls", default=os.path.join(root, "data", "external", "cancerhotspots", "hotspots_v2.xls"))
    ap.add_argument("--out", default=os.path.join(root, "data", "interim", "cancerhotspots_labels.tsv"))
    ap.add_argument("--qvalue-max", type=float, default=0.05)
    args = ap.parse_args()

    df = build_hotspot_labels(args.xls, args.qvalue_max)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print("[cancerhotspots] missense alt-specific rows: %d" % len(df))
    print("  is_hotspot (q<=%.3g): %d" % (args.qvalue_max, int(df["is_hotspot"].sum())))
    print("  unique genes: %d" % df["gene"].nunique())
    print("[cancerhotspots] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
