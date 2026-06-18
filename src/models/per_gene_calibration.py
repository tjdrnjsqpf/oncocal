"""Per-gene calibration — quantify that the score→oncogenic threshold differs per gene + a calibration map (cancer version of Tejura 2024, a differentiator).

Core claim: a single global threshold masks gene heterogeneity. Cancer oncogenicity should be calibrated per gene.
  - For both the gene-group OOF predictions (baseline_oof.tsv, leak-free) and AlphaMissense raw,
    compute the optimal per-gene threshold (Youden J) and show how widely thresholds scatter.
  - Compare balanced accuracy of global vs per-gene thresholds → clinical benefit of per-gene calibration.
  - Output: a per-gene table "in this gene, score ≥ T means oncogenic" (clinical interpretation).

Usage: python -m src.models.per_gene_calibration
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBL = os.path.join(ROOT, "results", "tables")


def youden_threshold(y, s):
    """Threshold maximizing Youden J."""
    order = np.argsort(s)
    ys, ss = y[order], s[order]
    P, N = ys.sum(), len(ys) - ys.sum()
    if P == 0 or N == 0:
        return np.nan, np.nan
    # TPR-FPR at each candidate threshold (unique score)
    best_j, best_t = -1, np.nan
    tp = P; fp = N
    prev = None
    for i in range(len(ss)):
        # Raising the threshold to ss[i] predicts everything ≤ ss[i] as negative
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
    # Global threshold
    gt_model, _ = youden_threshold(df.y.values, df.model_pred.values)
    gt_am, _ = youden_threshold(df.y.values, df.am_pathogenicity.values)
    print("Global Youden threshold: model=%.3f, AlphaMissense=%.3f" % (gt_model, gt_am))

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

    print("\n[per-gene calibration] genes with sufficient labels: %d" % len(cal))
    print("  AlphaMissense per-gene threshold distribution: median %.3f, range [%.3f, %.3f], IQR [%.3f, %.3f]"
          % (cal.am_threshold.median(), cal.am_threshold.min(), cal.am_threshold.max(),
             cal.am_threshold.quantile(.25), cal.am_threshold.quantile(.75)))
    print("  → per-gene heterogeneity that a single global threshold(%.3f) cannot capture (this is the core claim)" % gt_am)
    print("\n  Genes with low threshold (low scores still oncogenic) Top5:")
    print(cal.head(5)[["gene", "role_class", "n_pos", "am_threshold", "model_threshold"]].to_string(index=False))
    print("  Genes with high threshold (only high scores oncogenic) Top5:")
    print(cal.tail(5)[["gene", "role_class", "n_pos", "am_threshold", "model_threshold"]].to_string(index=False))

    # Global vs per-gene threshold: balanced accuracy
    def bacc(score_col, thr_col):
        ok = cal.set_index("gene")[thr_col].to_dict()
        sub = df[df.gene.isin(ok)]
        pred_global = (sub[score_col] >= (gt_model if score_col == "model_pred" else gt_am)).astype(int)
        pred_gene = (sub[score_col] >= sub.gene.map(ok)).astype(int)
        return (balanced_accuracy_score(sub.y, pred_global),
                balanced_accuracy_score(sub.y, pred_gene))
    for col, thr in [("am_pathogenicity", "am_threshold"), ("model_pred", "model_threshold")]:
        bg, bp = bacc(col, thr)
        print("\n  %s: balanced-acc global threshold %.3f → per-gene threshold %.3f (Δ%+.3f)"
              % (col, bg, bp, bp - bg))
    print("\n[per-gene] -> results/tables/per_gene_calibration.tsv")


if __name__ == "__main__":
    main()
