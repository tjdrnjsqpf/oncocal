"""ClinVar variant_summary -> benign negative 라벨 (+ 보조 oncogenicity 신호).

설계 근거 (phase1_execution.md C):
  - negative = ClinVar Benign / Likely benign (germline 임상 분류) — gnomAD common 과 함께 negative 정의.
  - GRCh38 행만, Type='single nucleotide variant', missense(Name 의 p.change)만.
  - ⚠️ leakage 주의: 많은 germline 도구(REVEL 등)가 ClinVar benign 으로 학습됨 →
    이 negative 는 그 도구들과 변이를 공유할 수 있다. 평가 split(gene/position) + ablation 으로 대응.
  - 보너스: variant_summary 에 'Oncogenicity' 컬럼(체세포 oncogenicity 분류)이 있어
    'Oncogenic/Likely oncogenic' 를 보조 positive 신호, 'Benign' 을 보조 negative 신호로 보존.
    (단 1차 positive 백본은 CIViC+hotspots+COSMIC; ClinVar oncogenicity 는 cross-check 용.)

Name 예: 'NM_000355.4(TCN2):c.776G>C (p.Arg259Pro)'  -> 3-letter p.change.
'p.Asn378=' 같은 synonymous 는 parse 실패로 자동 제외(alt AA 없음).

출력: data/interim/clinvar_labels.tsv (pkey 단위)
사용: python -m src.data.clinvar_labels
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys
from collections import defaultdict

import pandas as pd

try:
    from .variant_keys import parse_protein_change, is_missense_change, protein_key
except ImportError:
    sys.path.insert(0, os.path.dirname(__file__))
    from variant_keys import parse_protein_change, is_missense_change, protein_key

BENIGN_TOKENS = ("Benign", "Likely benign")  # ClinicalSignificance 부분일치
# review status -> star (대표적인 것만; 그 외는 0성으로 보수적 처리)
STAR = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
}


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _is_benign(cs: str) -> bool:
    if not cs:
        return False
    # 'Pathogenic' 류가 섞인 conflicting 은 제외, benign 토큰만 포함되어야 함
    low = cs.lower()
    if "pathogenic" in low or "uncertain" in low or "oncogenic" in low:
        return False
    return any(t.lower() in low for t in BENIGN_TOKENS)


def build_clinvar_labels(path: str, min_star: int = 1) -> pd.DataFrame:
    import csv
    agg = {}  # pkey -> dict
    onc_signal = defaultdict(set)  # pkey -> set of oncogenicity classifications
    with gzip.open(path, "rt") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("Assembly") != "GRCh38":
                continue
            if row.get("Type") != "single nucleotide variant":
                continue
            gene = (row.get("GeneSymbol") or "").strip()
            if not gene or gene == "-":
                continue
            pc = parse_protein_change(row.get("Name", ""))
            if pc is None:
                continue
            ref, pos, alt = pc
            if not is_missense_change(ref, alt):
                continue
            pkey = protein_key(gene, ref, pos, alt)

            # 보조 oncogenicity 신호 보존
            onc = (row.get("Oncogenicity") or "").strip()
            if onc and onc not in ("-", ""):
                onc_signal[pkey].add(onc)

            cs = row.get("ClinicalSignificance", "")
            if not _is_benign(cs):
                continue
            star = STAR.get((row.get("ReviewStatus") or "").strip().lower(), 0)
            cur = agg.get(pkey)
            if cur is None or star > cur["clinvar_review_star"]:
                agg[pkey] = {
                    "pkey": pkey, "gene": gene, "ref_aa": ref, "pos": pos, "alt_aa": alt,
                    "protein_change": "%s%d%s" % (ref, pos, alt),
                    "clinvar_significance": cs,
                    "clinvar_review_star": star,
                    "clinvar_variation_id": row.get("VariationID", ""),
                }
    df = pd.DataFrame(agg.values())
    if df.empty:
        raise SystemExit("ClinVar: benign missense 0건 — 입력/필터 확인.")
    df["clinvar_oncogenicity"] = df["pkey"].map(
        lambda k: "|".join(sorted(onc_signal.get(k, set()))) or "")
    df["clinvar_benign"] = df["clinvar_review_star"] >= min_star
    df = df.sort_values("clinvar_review_star", ascending=False).reset_index(drop=True)
    return df


def main():
    ap = argparse.ArgumentParser(description="ClinVar -> benign negative labels (GRCh38 missense).")
    root = _repo_root()
    ap.add_argument("--clinvar", default=os.path.join(root, "data", "external", "clinvar", "variant_summary.txt.gz"))
    ap.add_argument("--out", default=os.path.join(root, "data", "interim", "clinvar_labels.tsv"))
    ap.add_argument("--min-star", type=int, default=1, help="clinvar_benign=True 로 칠 최소 review star")
    args = ap.parse_args()

    df = build_clinvar_labels(args.clinvar, args.min_star)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print("[clinvar] benign missense pkey rows: %d" % len(df))
    print("  clinvar_benign (>=%d star): %d" % (args.min_star, int(df["clinvar_benign"].sum())))
    print("  with somatic Oncogenicity annotation: %d" % int((df["clinvar_oncogenicity"] != "").sum()))
    print("[clinvar] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
