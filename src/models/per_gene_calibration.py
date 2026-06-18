"""Per-gene 보정 — 점수→oncogenic 임계가 유전자마다 다름을 정량화 + 보정 맵 (Tejura 2024 암 버전, 차별점).

핵심 주장: 전역 단일 임계는 유전자 이질성을 가린다. 암 oncogenicity는 유전자별로 보정해야 한다.
  - gene-group OOF 예측(baseline_oof.tsv, 누수 없음) + AlphaMissense raw 각각에 대해
    유전자별 최적 임계(Youden J)를 구하고, 임계가 얼마나 흩어지는지 보인다.
  - 전역 임계 vs per-gene 임계의 balanced accuracy 비교 → per-gene 보정의 임상적 이득.
  - 산출: 유전자별 "이 유전자에선 점수 ≥ T 면 oncogenic" 표(임상 해석).

사용: python -m src.models.per_gene_calibration
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBL = os.path.join(ROOT, "results", "tables")


def youden_threshold(y, s):
    """Youden J 최대화 임계."""
    order = np.argsort(s)
    ys, ss = y[order], s[order]
    P, N = ys.sum(), len(ys) - ys.sum()
    if P == 0 or N == 0:
        return np.nan, np.nan
    # 각 후보 임계(고유 점수)에서 TPR-FPR
    best_j, best_t = -1, np.nan
    tp = P; fp = N
    prev = None
    for i in range(len(ss)):
        # 임계를 ss[i]로 올리면 ss[i] 이하를 음성 예측
        if prev is not None and ss[i] != prev:
            tpr = tp / P; fpr = fp / N
            j = tpr - fpr
            if j > best_j:
                best_j, best_t = j, prev
        if ys[i] == 1:
            tp -= 1
        else:
            fp -= 1
        prev = ss[i]
    return best_t, best_j


def main():
    df = pd.read_csv(os.path.join(ROOT, "data", "processed", "baseline_oof.tsv"), sep="\t")
    df = df.dropna(subset=["model_pred", "am_pathogenicity"])
    # 전역 임계
    gt_model, _ = youden_threshold(df.y.values, df.model_pred.values)
    gt_am, _ = youden_threshold(df.y.values, df.am_pathogenicity.values)
    print("전역 Youden 임계: model=%.3f, AlphaMissense=%.3f" % (gt_model, gt_am))

    rows = []
    for gene, g in df.groupby("gene"):
        npos, nneg = int((g.y == 1).sum()), int((g.y == 0).sum())
        if npos < 5 or nneg < 5:
            continue
        tm, _ = youden_threshold(g.y.values, g.model_pred.values)
        ta, _ = youden_threshold(g.y.values, g.am_pathogenicity.values)
        try:
            au_m = roc_auc_score(g.y, g.model_pred)
            au_a = roc_auc_score(g.y, g.am_pathogenicity)
        except ValueError:
            au_m = au_a = np.nan
        rows.append({"gene": gene, "role_class": g.role_class.iloc[0], "n_pos": npos, "n_neg": nneg,
                     "model_threshold": tm, "am_threshold": ta,
                     "model_auroc": au_m, "am_auroc": au_a})
    cal = pd.DataFrame(rows).sort_values("am_threshold")
    cal.to_csv(os.path.join(TBL, "per_gene_calibration.tsv"), sep="\t", index=False)

    print("\n[per-gene 보정] 충분 라벨 유전자: %d" % len(cal))
    print("  AlphaMissense per-gene 임계 분포: 중앙 %.3f, 범위 [%.3f, %.3f], IQR [%.3f, %.3f]"
          % (cal.am_threshold.median(), cal.am_threshold.min(), cal.am_threshold.max(),
             cal.am_threshold.quantile(.25), cal.am_threshold.quantile(.75)))
    print("  → 전역 단일 임계(%.3f)로는 가릴 수 없는 유전자별 이질성 (이게 핵심 주장)" % gt_am)
    print("\n  임계 낮은 유전자(저점수도 oncogenic) Top5:")
    print(cal.head(5)[["gene", "role_class", "n_pos", "am_threshold", "model_threshold"]].to_string(index=False))
    print("  임계 높은 유전자(고점수만 oncogenic) Top5:")
    print(cal.tail(5)[["gene", "role_class", "n_pos", "am_threshold", "model_threshold"]].to_string(index=False))

    # 전역 vs per-gene 임계: balanced accuracy
    def bacc(score_col, thr_col):
        ok = cal.set_index("gene")[thr_col].to_dict()
        sub = df[df.gene.isin(ok)]
        pred_global = (sub[score_col] >= (gt_model if score_col == "model_pred" else gt_am)).astype(int)
        pred_gene = (sub[score_col] >= sub.gene.map(ok)).astype(int)
        return (balanced_accuracy_score(sub.y, pred_global),
                balanced_accuracy_score(sub.y, pred_gene))
    for col, thr in [("am_pathogenicity", "am_threshold"), ("model_pred", "model_threshold")]:
        bg, bp = bacc(col, thr)
        print("\n  %s: balanced-acc 전역임계 %.3f → per-gene임계 %.3f (Δ%+.3f)"
              % (col, bg, bp, bp - bg))
    print("\n[per-gene] -> results/tables/per_gene_calibration.tsv")


if __name__ == "__main__":
    main()
