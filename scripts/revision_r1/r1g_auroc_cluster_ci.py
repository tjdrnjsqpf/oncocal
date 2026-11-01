"""R1-G — Gene-clustered bootstrap CI for AUROC, replacing the variant-level CI.

Why
---
The submitted manuscript reports OncoCal AUROC 0.928 [0.923, 0.933] vs best single tool
0.870 [0.862, 0.878], dAUROC +0.058 [0.052, 0.064]. Those CIs come from
scripts/revision/s5_logo_ci.py, which resamples VARIANTS (`RNG.randint(0, n, n)`).

But the manuscript's mechanism analysis deliberately switched to GENE-clustered inference
to remove pseudo-replication (naive p = 1e-219 -> gene-clustered p = 4.8e-20). Mixing
inference units within one paper is a reviewable inconsistency, and R1-F showed the per-gene
calibration gain evaporates once genes are the unit. This script applies the same gene
clustering to the AUROC comparison so the paper is internally consistent.

Both bootstraps are reported side by side so the difference is explicit and auditable.

Run:  python r1g_auroc_cluster_ci.py
Out:  r1_out/r1g_auroc_cluster_ci.tsv
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)

B = 2000
SEED = 20261031
CORE = ["oncogene", "TSG", "both"]
PREDS = ["p_model_tools", "p_model_full", "AlphaMissense", "VARITY_R"]


def find_table(name):
    tried = []
    for sub in ("results/tables", "paper_handoff/tables", "tables"):
        p = os.path.join(ROOT, sub, name)
        tried.append(p)
        if os.path.exists(p):
            return p
    raise FileNotFoundError("%s not found; tried:\n  %s" % (name, "\n  ".join(tried)))


def auroc(y, s):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, s)


def main():
    oof = pd.read_csv(find_table("oof_predictions.tsv"), sep="\t")
    oof["gene"] = oof.pkey.str.split(":").str[0]
    for c in PREDS:
        if c in oof.columns:
            oof[c] = oof[c].fillna(0.5)   # same NA handling as the published s5 script
    print("oof n=%d pos=%d genes=%d" % (len(oof), int(oof.y_true.sum()), oof.gene.nunique()))

    rows = []
    for set_nm, sub in [("all_roles", oof),
                        ("core_roles", oof[oof.role_class.isin(CORE)].copy())]:
        y = sub.y_true.values
        avail = [c for c in PREDS if c in sub.columns]
        point = {c: auroc(y, sub[c].values) for c in avail}
        best = max(("AlphaMissense", "VARITY_R"), key=lambda c: point.get(c, -1))
        print("\n[%s] n=%d pos=%d genes=%d | best single tool = %s (%.4f)"
              % (set_nm, len(sub), int(y.sum()), sub.gene.nunique(), best, point[best]))

        # index variants by gene for cluster resampling
        gene_list = np.array(sorted(sub.gene.unique()))
        idx_by_gene = {g: np.flatnonzero((sub.gene == g).values) for g in gene_list}
        mat = {c: sub[c].values for c in avail}

        for scheme in ("variant", "gene_clustered"):
            rng = np.random.RandomState(SEED)
            boot = {c: [] for c in avail}
            gaps = []
            n = len(sub)
            for _ in range(B):
                if scheme == "variant":
                    idx = rng.randint(0, n, n)
                else:
                    pick = rng.choice(gene_list, len(gene_list), replace=True)
                    idx = np.concatenate([idx_by_gene[g] for g in pick])
                yb = y[idx]
                if len(np.unique(yb)) < 2:
                    continue
                vals = {}
                for c in avail:
                    vals[c] = auroc(yb, mat[c][idx])
                    boot[c].append(vals[c])
                gaps.append(vals["p_model_tools"] - vals[best])
            gaps = np.array(gaps)
            glo, ghi = np.percentile(gaps, [2.5, 97.5])
            for c in avail:
                arr = np.array(boot[c])
                lo, hi = np.percentile(arr, [2.5, 97.5])
                rows.append({"role_set": set_nm, "scheme": scheme, "quantity": c,
                             "point": round(point[c], 4), "ci_lo": round(float(lo), 4),
                             "ci_hi": round(float(hi), 4), "n_boot": len(arr)})
            rows.append({"role_set": set_nm, "scheme": scheme,
                         "quantity": "dAUROC_OncoCal_minus_%s" % best,
                         "point": round(point["p_model_tools"] - point[best], 4),
                         "ci_lo": round(float(glo), 4), "ci_hi": round(float(ghi), 4),
                         "n_boot": len(gaps),
                         "P_gt_0": round(float((gaps > 0).mean()), 4),
                         "excludes_zero": "YES" if glo > 0 else "no"})
            print("   %-15s OncoCal %.4f [%.4f, %.4f] | %s %.4f [%.4f, %.4f] | dAUROC %+.4f [%+.4f, %+.4f] excl0=%s"
                  % (scheme, point["p_model_tools"],
                     np.percentile(boot["p_model_tools"], 2.5), np.percentile(boot["p_model_tools"], 97.5),
                     best, point[best],
                     np.percentile(boot[best], 2.5), np.percentile(boot[best], 97.5),
                     point["p_model_tools"] - point[best], glo, ghi,
                     "YES" if glo > 0 else "no"))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "r1g_auroc_cluster_ci.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 230)
    print("\n--- summary ---")
    print(res.to_string(index=False))
    print("\n[saved] r1_out/r1g_auroc_cluster_ci.tsv")


if __name__ == "__main__":
    main()
