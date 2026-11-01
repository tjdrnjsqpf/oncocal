"""R1-F — Gene-clustered bootstrap CI on the per-gene calibration gain.

Question this answers
---------------------
With the corrected threshold function, within-gene leave-one-variant-out calibration
improves balanced accuracy by only +0.0096 (AlphaMissense) / +0.0044 (OncoCal) on core
roles. Is that distinguishable from zero? The manuscript's wording depends on the answer:
a CI excluding zero supports "a small but consistent gain"; a CI spanning zero means we
must write "no measurable gain".

Method
------
Gene is the independent unit (same convention as scripts/revision/s2_clustered.py and
s5_logo_ci.py), so we resample GENES with replacement, B times.

Within each replicate:
  - the global Youden threshold is REFITTED on the resampled pool, because it is a
    pooled estimate and must be subject to the same resampling as the per-gene arm;
  - the within-gene LOVO prediction for a variant is a property of its own gene's data,
    which gene-level resampling leaves intact, so those predictions are precomputed once;
  - balanced accuracy is computed for both arms on the resampled pool and differenced,
    giving a paired delta per replicate.

Reported: percentile 95% CI of the delta, plus P(delta > 0).

Run:  python r1f_pergene_ci.py
Out:  r1_out/r1f_pergene_ci.tsv
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)

CORE = ["oncogene", "TSG", "both"]
MIN_PER_CLASS = 5
B = 2000
SEED = 20261031
SCORES = [("am_pathogenicity", "AlphaMissense"), ("model_pred", "OncoCal")]


def youden_threshold(y, s):
    """Score value T maximising Youden's J for the rule: predict positive iff score >= T."""
    y = np.asarray(y)
    s = np.asarray(s)
    P, N = y.sum(), len(y) - y.sum()
    if P == 0 or N == 0:
        return np.nan
    order = np.argsort(s, kind="mergesort")
    ys, ss = y[order], s[order]
    best_j, best_t, tp, fp, i, n = -np.inf, np.nan, P, N, 0, len(ss)
    while i < n:
        t = ss[i]
        j = tp / P - fp / N
        if j > best_j:
            best_j, best_t = j, t
        while i < n and ss[i] == t:
            if ys[i] == 1:
                tp -= 1
            else:
                fp -= 1
            i += 1
    return best_t


def eligible(df):
    ok = df.groupby("gene").y.agg(lambda x: ((x == 1).sum() >= MIN_PER_CLASS
                                             and (x == 0).sum() >= MIN_PER_CLASS))
    return set(ok[ok].index)


def main():
    df = pd.read_csv(os.path.join(ROOT, "data", "processed", "baseline_oof.tsv"), sep="\t")
    df = df.dropna(subset=["model_pred", "am_pathogenicity"])
    rows = []

    for set_nm, sub in [("core_roles", df[df.role_class.isin(CORE)].copy()),
                        ("all_roles", df)]:
        genes = eligible(sub)
        sub_e = sub[sub.gene.isin(genes)].copy()
        print("\n[%s] genes=%d variants=%d" % (set_nm, len(genes), len(sub_e)))

        for col, nm in SCORES:
            gt_full = youden_threshold(sub.y.values, sub[col].values)

            # ---- precompute per-variant LOVO predictions (gene-local, resampling-invariant)
            per_gene = {}
            for gene, g in sub_e.groupby("gene"):
                y, s = g.y.values, g[col].values
                lovo = np.empty(len(g), dtype=int)
                for i in range(len(g)):
                    m = np.ones(len(g), bool)
                    m[i] = False
                    t = youden_threshold(y[m], s[m])
                    if np.isnan(t):
                        t = gt_full
                    lovo[i] = int(s[i] >= t)
                per_gene[gene] = dict(y=y, s=s, lovo=lovo)

            gene_list = np.array(sorted(per_gene))
            y_all = np.concatenate([per_gene[g]["y"] for g in gene_list])
            s_all = np.concatenate([per_gene[g]["s"] for g in gene_list])
            lovo_all = np.concatenate([per_gene[g]["lovo"] for g in gene_list])
            point = (balanced_accuracy_score(y_all, lovo_all)
                     - balanced_accuracy_score(y_all, (s_all >= gt_full).astype(int)))

            # ---- gene-clustered bootstrap ----
            rng = np.random.RandomState(SEED)
            deltas = []
            for _ in range(B):
                pick = rng.choice(gene_list, len(gene_list), replace=True)
                y_b = np.concatenate([per_gene[g]["y"] for g in pick])
                s_b = np.concatenate([per_gene[g]["s"] for g in pick])
                l_b = np.concatenate([per_gene[g]["lovo"] for g in pick])
                if len(np.unique(y_b)) < 2:
                    continue
                gt_b = youden_threshold(y_b, s_b)        # refit global on the replicate
                if np.isnan(gt_b):
                    continue
                deltas.append(balanced_accuracy_score(y_b, l_b)
                              - balanced_accuracy_score(y_b, (s_b >= gt_b).astype(int)))
            deltas = np.array(deltas)
            lo, hi = np.percentile(deltas, [2.5, 97.5])
            p_gt0 = float((deltas > 0).mean())
            excl = "YES" if lo > 0 else "no"
            rows.append({"role_set": set_nm, "predictor": nm, "n_gene": len(gene_list),
                         "n_variant": len(y_all), "delta_point": round(point, 4),
                         "boot_median": round(float(np.median(deltas)), 4),
                         "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
                         "P_delta_gt_0": round(p_gt0, 4), "excludes_zero": excl,
                         "n_boot": len(deltas)})
            print("   %-13s delta %+.4f  95%% CI [%+.4f, %+.4f]  P(>0)=%.3f  excludes 0: %s"
                  % (nm, point, lo, hi, p_gt0, excl))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "r1f_pergene_ci.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n--- summary ---")
    print(res.to_string(index=False))
    print("\n[saved] r1_out/r1f_pergene_ci.tsv")


if __name__ == "__main__":
    main()
