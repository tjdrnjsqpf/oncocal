"""도구 supervised/unsupervised 2-way 재분류 + driver 과소평가율 재계산 (R3).

원칙(사용자 정의):
  supervised        = 임상/병원성 라벨(ClinVar/HGMD/질병변이) 또는 proxy 판별 라벨로 학습한 판별모델
  unsupervised/self = 임상라벨 없이 보존도·진화·생성형 언어모델·집단유전학으로 점수 산출
경계사례는 실제 학습방식 근거로 판정(아래 RATIONALE 로그). 추정 금지 — 문헌상 학습방식 근거.
출력: driver_miss_rate.tsv 에 type2 추가, 두 그룹 평균 재계산.
사용: python -m src.evaluation.reclassify_tools
"""
from __future__ import annotations
import os
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FILES = [os.path.join(ROOT, "results", "tables", "driver_miss_rate.tsv"),
         os.path.join(ROOT, "paper_handoff", "tables", "driver_miss_rate.tsv")]

# 임상라벨 없이 점수 산출 → unsupervised/self-supervised (파일의 short tool명 기준)
UNSUP = {
    # 보존도/진화 (라벨 없음)
    "GERPpp_RS", "GERP_92_mammals", "phyloP17way_primate", "phyloP100way_vertebrate",
    "phyloP470way_mammalian", "phastCons17way_primate", "phastCons100way_vertebrate",
    "phastCons470way_mammalian", "SIFT", "SIFT4G", "PROVEAN", "MutationAssessor",
    # 생성형/자기지도 언어·진화 모델
    "ESM1b", "popEVE", "AlphaMissense",
    # 집단유전학(빈도/선택계수)·비지도 통합
    "bStatistic", "MisFit_D", "MisFit_S", "Eigen-raw_coding", "Eigen-PC-raw_coding",
    # 집단(영장류 common) 라벨 기반 (임상라벨 아님)
    "PrimateAI",
}
# 경계사례 판정 근거 로그
RATIONALE = {
    "AlphaMissense": "unsup — 단백질 LM + AlphaFold 구조, 집단 빈도(common=benign) 약지도. 임상/ClinVar 라벨 미사용(Cheng 2023, Science).",
    "PrimateAI": "unsup — 영장류/집단 common 변이를 benign으로 쓰는 판별 CNN. 임상 병원성 라벨(ClinVar/HGMD) 미사용(Sundaram 2018, Nat Genet). [판별모델이나 임상라벨 아님 → 원칙상 unsup]",
    "Eigen-raw_coding": "unsup — 라벨 없이 주석들을 스펙트럼/고유분해로 통합하는 비지도 메타점수(Ionita-Laza 2016, Nat Genet).",
    "Eigen-PC-raw_coding": "unsup — Eigen의 PC 변형, 동일하게 비지도.",
    "MisFit_D": "unsup — 집단 대립유전자빈도로 선택계수/해로움을 추정하는 집단유전학 모델. 임상라벨 미사용.",
    "MisFit_S": "unsup — 위와 동일(선택계수 S).",
    "DANN": "supervised — CADD와 동일한 판별 학습셋(관측 proxy-benign vs 시뮬 proxy-deleterious)으로 학습한 신경망(Quang 2015). CADD를 supervised로 두는 본 분류와 일관. (엄격히 '임상라벨'만 기준이면 CADD/DANN은 unsup으로 이동 — 그 경우 재계산 필요.)",
}


def classify():
    df = pd.read_csv(FILES[0] if os.path.exists(FILES[0]) else FILES[1], sep="\t")
    df["type2"] = df["tool"].map(lambda t: "unsupervised" if t in UNSUP else "supervised")
    return df


def main():
    df = classify()
    for f in FILES:
        if os.path.exists(f):
            df.to_csv(f, sep="\t", index=False)

    print("=== 경계사례 판정 근거 ===")
    for k, v in RATIONALE.items():
        print("  [%s] %s" % (k, v))

    print("\n=== type2 그룹별 도구 ===")
    for t2 in ["unsupervised", "supervised"]:
        g = df[df.type2 == t2].sort_values("miss_rate(<0.5)", ascending=False)
        print("\n-- %s (n=%d) --" % (t2, len(g)))
        print("  도구:", ", ".join(g.tool.tolist()))
        print("  평균 under-scoring(rec≥50 driver, rankscore<0.5): %.3f (%.1f%%)"
              % (g["miss_rate(<0.5)"].mean(), 100 * g["miss_rate(<0.5)"].mean()))

    u = df[df.type2 == "unsupervised"]["miss_rate(<0.5)"].mean()
    s = df[df.type2 == "supervised"]["miss_rate(<0.5)"].mean()
    print("\n=== R3 새 수치(현 14.3/9.1 대체) ===")
    print("  unsupervised %.1f%% vs supervised %.1f%% (차 %.1f%%p)" % (100 * u, 100 * s, 100 * (u - s)))
    print("\n[done] type2 추가 -> driver_miss_rate.tsv (results/ + paper_handoff/)")


if __name__ == "__main__":
    main()
