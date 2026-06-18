"""통합 데이터셋 빌드 — 패널 missense universe에 점수·라벨·recurrence를 pkey로 결합.

아키텍처 (2026-06-16, COSMIC 도착 후):
  변이 universe = AlphaMissense 패널 missense (data/interim/alphamissense_panel.tsv)
    → 패널 내 '가능한 모든 missense' + flagship 점수. dbNSFP 도착 시 다도구 점수 컬럼 확장.
  여기에 left-join:
    - CIViC      : civic_oncogenic, civic_mechanism
    - hotspots   : is_hotspot, hotspot_qvalue
    - ClinVar    : clinvar_benign
    - COSMIC rec : cosmic_sample_count (MutantCensus, 라벨 정의용)
    - panel      : role_class, mech_proxy (gene-level)

라벨 규칙 (civic-label 한계 + JAK2 V617F 구제 반영):
  positive: is_hotspot | civic_oncogenic | (cosmic_sample_count >= --recurrence-min)
            └ JAK2 V617F 는 hotspot 누락이지만 COSMIC recurrence(42,930)로 구제됨.
  negative: clinvar_benign  且 positive 신호 없음
  conflict: positive·negative 동시 → 학습 제외
  unknown : 신호 없음 → 학습 제외 (절대 negative 금지 — 핵심 함정)

기전: civic_mechanism(변이 단위) 우선, 없으면 panel mech_proxy(유전자 단위, 약한 prior).

⚠️ circularity(§6): cosmic_sample_count 는 **라벨 정의**에 사용 → 기본 feature 에서 제외.
   feature 로 쓸 땐 GenomeScreens(unbiased) recurrence 로 ablation arm 만 사용.
⚠️ leakage: gene/position group-split 은 make_group_split.py 에서 별도 생성.

출력: data/processed/dataset.parquet (전체 universe) + labels_summary.tsv (라벨된 부분집합)
사용: python -m src.data.build_dataset --recurrence-min 10
"""
from __future__ import annotations
import argparse
import os

import pandas as pd

INTERIM = "data/interim"
PROCESSED = "data/processed"


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _load(name: str, root: str):
    p = os.path.join(root, INTERIM, name)
    if not os.path.exists(p):
        print("[build] (없음, skip): %s" % p)
        return None
    return pd.read_csv(p, sep="\t")


def build(root: str, recurrence_min: int, gnomad_common_af: float = 1e-3) -> pd.DataFrame:
    am = _load("alphamissense_panel.tsv", root)
    if am is None:
        raise SystemExit("alphamissense_panel.tsv 없음 — alphamissense_join 먼저 실행.")
    base = am.copy()

    civic = _load("civic_labels.tsv", root)
    hot = _load("cancerhotspots_labels.tsv", root)
    clin = _load("clinvar_labels.tsv", root)
    rec = _load("cosmic_recurrence__mutantcensus.tsv", root)
    panel = _load("cancer_gene_panel.tsv", root)
    dbnsfp = _load("dbnsfp_panel.tsv", root)  # 도착 시 자동 결합(없으면 skip)

    if civic is not None:
        base = base.merge(civic[["pkey", "civic_oncogenic", "civic_mechanism", "civic_onc_support"]],
                          on="pkey", how="left")
    if hot is not None:
        base = base.merge(hot[["pkey", "is_hotspot", "hotspot_qvalue", "hotspot_alt_count"]],
                          on="pkey", how="left")
    if clin is not None:
        base = base.merge(clin[["pkey", "clinvar_benign", "clinvar_significance"]],
                          on="pkey", how="left")
    if rec is not None:
        # ⚠️ recurrence 는 GRCh38 유전체 좌표 단위(isoform 번호 문제 회피) → AM(GRCh38)과 좌표로 join.
        gk = ["grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt"]
        for c in ["grch38_chr", "grch38_ref", "grch38_alt"]:
            base[c] = base[c].astype(str)
            rec[c] = rec[c].astype(str)
        base["grch38_pos"] = base["grch38_pos"].astype(int)
        rec["grch38_pos"] = rec["grch38_pos"].astype(int)
        base = base.merge(rec[gk + ["cosmic_sample_count"]], on=gk, how="left")
    if panel is not None:
        base = base.merge(panel[["gene", "role_class", "mech_proxy", "tier"]], on="gene", how="left")
    if dbnsfp is not None:
        # 다도구 점수/보존도/gnomAD AF — GRCh38 좌표 join (40+ 도구 컬럼 추가, 구조 변경 없음)
        gk = ["grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt"]
        for c in ["grch38_chr", "grch38_ref", "grch38_alt"]:
            base[c] = base[c].astype(str); dbnsfp[c] = dbnsfp[c].astype(str)
        base["grch38_pos"] = base["grch38_pos"].astype(int)
        dbnsfp["grch38_pos"] = dbnsfp["grch38_pos"].astype(int)
        drop = [c for c in ["aaref", "aaalt", "dbnsfp_gene"] if c in dbnsfp.columns]
        base = base.merge(dbnsfp.drop(columns=drop), on=gk, how="left")
        print("[build] dbNSFP 점수 컬럼 %d개 결합" % sum(c.startswith("dbnsfp_") for c in base.columns))

    for c, d in [("civic_oncogenic", False), ("is_hotspot", False), ("clinvar_benign", False)]:
        base[c] = base.get(c, d)
        base[c] = base[c].fillna(False).astype(bool)
    base["cosmic_sample_count"] = base.get("cosmic_sample_count", 0)
    base["cosmic_sample_count"] = base["cosmic_sample_count"].fillna(0).astype(int)

    # ---- 라벨 ----
    base["recurrent_driver"] = base["cosmic_sample_count"] >= recurrence_min
    base["pos_signal"] = base["is_hotspot"] | base["civic_oncogenic"] | base["recurrent_driver"]
    # gnomAD common variant = benign negative 보강 (dbNSFP 있을 때만)
    af_col = "dbnsfp_gnomAD4.1_joint_AF"
    if af_col in base.columns:
        base["gnomad_common"] = pd.to_numeric(base[af_col], errors="coerce").fillna(0) >= gnomad_common_af
    else:
        base["gnomad_common"] = False
    base["neg_signal"] = base["clinvar_benign"] | base["gnomad_common"]
    base["label_conflict"] = base["pos_signal"] & base["neg_signal"]

    def _lab(r):
        if r["label_conflict"]:
            return "conflict"
        if r["pos_signal"]:
            return "positive"
        if r["neg_signal"]:
            return "negative"
        return "unknown"
    base["oncogenic_label"] = base.apply(_lab, axis=1)
    base["y"] = base["oncogenic_label"].map({"positive": 1, "negative": 0})

    # ---- 기전 ----
    base["mechanism"] = base.get("civic_mechanism")
    base["mechanism"] = base["mechanism"].fillna(base.get("mech_proxy"))
    return base


def summarize(df: pd.DataFrame, recurrence_min: int) -> None:
    print("[build] universe (panel missense w/ AM): %d  (genes=%d)" % (len(df), df["gene"].nunique()))
    print(df["oncogenic_label"].value_counts().to_string())
    tr = df[df["y"].notna()]
    print("\n학습가능 pos=%d / neg=%d  (recurrence_min=%d)"
          % (int((tr["y"] == 1).sum()), int((tr["y"] == 0).sum()), recurrence_min))
    print("  positive 출처: hotspot=%d, CIViC=%d, recurrent(>=%d)=%d"
          % (int(df["is_hotspot"].sum()), int(df["civic_oncogenic"].sum()),
             recurrence_min, int(df["recurrent_driver"].sum())))
    # 핵심 driver 점검: 과소평가 미리보기
    canon = ["JAK2:p.V617F", "BRAF:p.V600E", "EGFR:p.L858R", "KRAS:p.G12D", "IDH1:p.R132H"]
    sub = df[df["pkey"].isin(canon)][["pkey", "am_pathogenicity", "am_class",
                                      "cosmic_sample_count", "is_hotspot", "oncogenic_label"]]
    if len(sub):
        print("\n[과소평가 미리보기] canonical driver — AM 점수 vs recurrence:")
        print(sub.to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description="Build unified panel dataset (AM universe + labels + recurrence).")
    root = _repo_root()
    ap.add_argument("--recurrence-min", type=int, default=10,
                    help="COSMIC distinct-sample 수 >= 이 값이면 recurrent driver(positive) 후보")
    ap.add_argument("--gnomad-common-af", type=float, default=1e-3,
                    help="gnomAD joint AF >= 이 값이면 common=benign negative(dbNSFP 있을 때)")
    ap.add_argument("--out-parquet", default=os.path.join(root, PROCESSED, "dataset.parquet"))
    ap.add_argument("--out-labeled", default=os.path.join(root, PROCESSED, "labels_summary.tsv"))
    args = ap.parse_args()

    df = build(root, args.recurrence_min, args.gnomad_common_af)
    os.makedirs(os.path.join(root, PROCESSED), exist_ok=True)
    df.to_parquet(args.out_parquet, index=False)
    df[df["y"].notna()].to_csv(args.out_labeled, sep="\t", index=False)
    summarize(df, args.recurrence_min)
    print("\n[build] wrote -> %s (전체), %s (라벨된 부분)" % (args.out_parquet, args.out_labeled))


if __name__ == "__main__":
    main()
