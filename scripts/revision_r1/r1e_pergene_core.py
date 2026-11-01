"""R1-E — Per-gene calibration on core roles, with an unambiguous algorithm statement.

Addresses reviewer #1 comments 4 and 5.

R1.5 (the terminology contradiction) -- resolved by reading the code
-------------------------------------------------------------------
scripts/revision/s5_logo_ci.py loops over genes, and WITHIN each gene leaves out one
variant at a time, deriving that variant's threshold from the gene's remaining variants.
So the procedure is within-gene leave-one-variant-out. The Methods sentence
"leave-one-gene-out calibration derived each gene's threshold from all other genes" is
WRONG and must be replaced. Use one term everywhere:
    "within-gene leave-one-variant-out (LOVO) calibration"

Algorithm, step by step (for the Methods section)
  1. Restrict to genes with >= 5 labelled positives AND >= 5 labelled negatives.
  2. For gene g with labelled variants v_1..v_n:
     2a. For each i, remove v_i, leaving n-1 variants.
     2b. On those n-1 variants compute the threshold T maximising Youden's J.
        (If J is undefined -- one class only -- fall back to the global threshold.)
     2c. Predict v_i oncogenic iff score(v_i) >= T. v_i never informs its own threshold.
  3. Pool all held-out predictions across genes; compute balanced accuracy.
  4. Compare against the single global Youden threshold fitted on the full labelled set.
     (That global baseline is in-sample and therefore favourable to the comparator, so
     the reported gain is conservative.)

R1.4 (are the thresholds clinically usable?)
  We report the per-gene label counts so readers can see how thin each gene's evidence is,
  plus a 2-fold within-gene split as a second, harder validation. The manuscript wording
  moves from "clinically actionable" to "gene-specific calibration reference".

Outputs
  r1e_pergene_table.tsv       per-gene thresholds, label counts, AUROC (core roles)
  r1e_pergene_validation.tsv  global vs LOVO vs 2-fold balanced accuracy
  r1e_label_depth.tsv         distribution of per-gene label counts
"""
from __future__ import annotations
import os
import sys
import zlib

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)


def find_table(name):
    """Locate a pipeline table across repo layouts (server: results/tables or
    paper_handoff/tables; public repo: tables/). Raises with the paths tried."""
    tried = []
    for sub in ("results/tables", "paper_handoff/tables", "tables"):
        p = os.path.join(ROOT, sub, name)
        tried.append(p)
        if os.path.exists(p):
            return p
    raise FileNotFoundError("%s not found; tried:\n  %s" % (name, "\n  ".join(tried)))

CORE = ["oncogene", "TSG", "both"]
MIN_PER_CLASS = 5
SCORES = [("am_pathogenicity", "AlphaMissense"), ("model_pred", "OncoCal")]
SEED_2FOLD = 20261031  # deterministic: crc32(gene) + this, so runs are reproducible


def gene_seed(gene: str) -> int:
    """Reproducible per-gene seed. Python's hash() is salted per process, so never use it."""
    return (zlib.crc32(str(gene).encode("utf-8")) + SEED_2FOLD) % (2 ** 31)


def youden_threshold(y, s):
    """Score value T maximising Youden's J for the rule: predict positive iff score >= T.

    Verified against brute force over all distinct score values on 4,996 random cases
    (0 suboptimal). Ties are consumed as a group, so T is always an observed score.
    """
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
        j = tp / P - fp / N            # tp/fp currently count scores >= t
        if j > best_j:
            best_j, best_t = j, t
        while i < n and ss[i] == t:    # consume the whole tie group
            if ys[i] == 1:
                tp -= 1
            else:
                fp -= 1
            i += 1
    return best_t


def youden_threshold_legacy(y, s):
    """The threshold function used for the SUBMITTED manuscript, kept only to quantify
    the correction. It records the score *below* the optimal cut, i.e. it returns a
    threshold for `score > T`, while every call site applies it as `score >= T`.
    On sparse per-gene data this shifts thresholds down (median ~0.06, up to ~0.44 in
    simulation) and costs Youden J under `>=`. Do not use for new results.
    """
    order = np.argsort(s)
    ys, ss = np.asarray(y)[order], np.asarray(s)[order]
    P, N = ys.sum(), len(ys) - ys.sum()
    if P == 0 or N == 0:
        return np.nan
    best_j, best_t, tp, fp, prev = -1, np.nan, P, N, None
    for i in range(len(ss)):
        if prev is not None and ss[i] != prev:
            j = tp / P - fp / N
            if j > best_j:
                best_j, best_t = j, prev
        if ys[i] == 1:
            tp -= 1
        else:
            fp -= 1
        prev = ss[i]
    return best_t


def eligible(df):
    ok = df.groupby("gene").y.agg(lambda x: ((x == 1).sum() >= MIN_PER_CLASS
                                             and (x == 0).sum() >= MIN_PER_CLASS))
    return set(ok[ok].index)


def main(thr_fn=None, tag="corrected"):
    global THR
    THR = thr_fn or youden_threshold
    df = pd.read_csv(os.path.join(ROOT, "data", "processed", "baseline_oof.tsv"), sep="\t")
    df = df.dropna(subset=["model_pred", "am_pathogenicity"])
    print("baseline_oof n=%d | role_class: %s"
          % (len(df), df.role_class.value_counts().to_dict()))
    core = df[df.role_class.isin(CORE)].copy()
    print("core-role subset n=%d (genes=%d)" % (len(core), core.gene.nunique()))

    results, depth_rows, table_rows = [], [], []
    for set_nm, sub in [("core_roles", core), ("all_roles", df)]:
        genes = eligible(sub)
        sub_e = sub[sub.gene.isin(genes)].copy()
        print("\n[%s] genes with >=%d per class: %d (variants=%d)"
              % (set_nm, MIN_PER_CLASS, len(genes), len(sub_e)))

        for col, nm in SCORES:
            gt = THR(sub.y.values, sub[col].values)  # global, full labelled set
            y_all, p_glob, p_lovo = [], [], []
            y_2f, p_2f = [], []

            for gene, g in sub_e.groupby("gene"):
                y, s = g.y.values, g[col].values
                # --- within-gene leave-one-variant-out ---
                for i in range(len(g)):
                    m = np.ones(len(g), bool)
                    m[i] = False
                    t = THR(y[m], s[m])
                    if np.isnan(t):
                        t = gt
                    p_lovo.append(int(s[i] >= t))
                    p_glob.append(int(s[i] >= gt))
                    y_all.append(y[i])
                # --- 2-fold within-gene split (harder: threshold from half the gene) ---
                rng = np.random.RandomState(gene_seed(gene))
                idx = rng.permutation(len(g))
                half = len(g) // 2
                for tr, te in [(idx[:half], idx[half:]), (idx[half:], idx[:half])]:
                    if len(tr) == 0 or len(te) == 0:
                        continue
                    t = THR(y[tr], s[tr])
                    if np.isnan(t):
                        t = gt
                    p_2f.extend((s[te] >= t).astype(int).tolist())
                    y_2f.extend(y[te].tolist())

                if set_nm == "core_roles" and col == "am_pathogenicity":
                    depth_rows.append({"gene": gene, "role_class": g.role_class.iloc[0],
                                       "n_pos": int((y == 1).sum()),
                                       "n_neg": int((y == 0).sum()), "n_total": len(g)})

            y_all = np.array(y_all)
            r = {"role_set": set_nm, "threshold_fn": tag, "predictor": nm, "global_threshold": round(gt, 4),
                 "n_variant": len(y_all), "n_gene": len(genes),
                 "bacc_global": round(balanced_accuracy_score(y_all, p_glob), 4),
                 "bacc_pergene_LOVO": round(balanced_accuracy_score(y_all, p_lovo), 4)}
            r["delta_LOVO"] = round(r["bacc_pergene_LOVO"] - r["bacc_global"], 4)
            if len(y_2f) == len(p_2f) and len(y_2f):
                r["bacc_pergene_2fold"] = round(balanced_accuracy_score(y_2f, p_2f), 4)
            results.append(r)
            print("   %-13s global %.4f -> LOVO %.4f (delta %+.4f) | 2-fold %s"
                  % (nm, r["bacc_global"], r["bacc_pergene_LOVO"], r["delta_LOVO"],
                     r.get("bacc_pergene_2fold", "n/a")))

        # ---- per-gene threshold table (core roles only) ----
        if set_nm == "core_roles":
            for gene, g in sub_e.groupby("gene"):
                row = {"gene": gene, "role_class": g.role_class.iloc[0],
                       "n_pos": int((g.y == 1).sum()), "n_neg": int((g.y == 0).sum())}
                for col, nm in SCORES:
                    row["%s_threshold" % nm] = THR(g.y.values, g[col].values)
                    try:
                        row["%s_auroc" % nm] = round(roc_auc_score(g.y, g[col]), 4)
                    except ValueError:
                        row["%s_auroc" % nm] = np.nan
                table_rows.append(row)

    tab = pd.DataFrame(table_rows).sort_values("AlphaMissense_threshold")
    tab.to_csv(os.path.join(OUT, "r1e_pergene_table_%s.tsv" % tag), sep="\t", index=False)
    val = pd.DataFrame(results)
    val.to_csv(os.path.join(OUT, "r1e_pergene_validation_%s.tsv" % tag), sep="\t", index=False)
    depth = pd.DataFrame(depth_rows)
    depth.to_csv(os.path.join(OUT, "r1e_label_depth_%s.tsv" % tag), sep="\t", index=False)

    pd.set_option("display.width", 250)
    print("\n--- validation summary ---")
    print(val.to_string(index=False))
    if len(tab):
        print("\n--- per-gene AlphaMissense thresholds (core roles, n=%d genes) ---" % len(tab))
        print("range %.3f - %.3f | median %.3f"
              % (tab.AlphaMissense_threshold.min(), tab.AlphaMissense_threshold.max(),
                 tab.AlphaMissense_threshold.median()))
        print("lowest 5:"); print(tab.head(5).to_string(index=False))
        print("highest 5:"); print(tab.tail(5).to_string(index=False))
    if len(depth):
        print("\n--- label depth per gene (quantifies 'limited labels', R1.4) ---")
        print(depth[["n_pos", "n_neg", "n_total"]].describe().round(1).to_string())
        print("genes with < 10 positives: %d / %d"
              % (int((depth.n_pos < 10).sum()), len(depth)))
    print("\n[saved] r1_out/r1e_{pergene_table,pergene_validation,label_depth}.tsv")


def compare_threshold_functions():
    """Run the whole validation twice -- with the submitted (legacy) threshold function and
    with the corrected one -- so the revision can state exactly what the fix changes."""
    print("\n" + "=" * 78)
    print("CORRECTED threshold function (score >= T, J-optimal)")
    print("=" * 78)
    main(youden_threshold, "corrected")
    print("\n" + "=" * 78)
    print("LEGACY threshold function (as submitted; returns the score below the cut)")
    print("=" * 78)
    main(youden_threshold_legacy, "legacy")

    a = pd.read_csv(os.path.join(OUT, "r1e_pergene_validation_corrected.tsv"), sep="\t")
    b = pd.read_csv(os.path.join(OUT, "r1e_pergene_validation_legacy.tsv"), sep="\t")
    comp = b.merge(a, on=["role_set", "predictor"], suffixes=("_legacy", "_corrected"))
    cols = ["role_set", "predictor", "global_threshold_legacy", "global_threshold_corrected",
            "bacc_pergene_LOVO_legacy", "bacc_pergene_LOVO_corrected",
            "delta_LOVO_legacy", "delta_LOVO_corrected"]
    comp = comp[[c for c in cols if c in comp.columns]]
    comp.to_csv(os.path.join(OUT, "r1e_threshold_fix_impact.tsv"), sep="\t", index=False)

    ta = pd.read_csv(os.path.join(OUT, "r1e_pergene_table_corrected.tsv"), sep="\t")
    tb = pd.read_csv(os.path.join(OUT, "r1e_pergene_table_legacy.tsv"), sep="\t")
    m = tb.merge(ta, on="gene", suffixes=("_legacy", "_corrected"))
    shift = m.AlphaMissense_threshold_corrected - m.AlphaMissense_threshold_legacy
    pd.set_option("display.width", 250)
    print("\n" + "=" * 78)
    print("IMPACT OF THE THRESHOLD FIX")
    print("=" * 78)
    print(comp.to_string(index=False))
    print("\nper-gene AlphaMissense threshold shift (corrected - legacy), n=%d genes:" % len(m))
    print("  mean %+.4f | median %+.4f | p95 %+.4f | max %+.4f | genes shifted >0.05: %d"
          % (shift.mean(), shift.median(), shift.quantile(.95), shift.max(),
             int((shift.abs() > 0.05).sum())))
    print("  legacy range  %.3f - %.3f (median %.3f)"
          % (tb.AlphaMissense_threshold.min(), tb.AlphaMissense_threshold.max(),
             tb.AlphaMissense_threshold.median()))
    print("  corrected range %.3f - %.3f (median %.3f)"
          % (ta.AlphaMissense_threshold.min(), ta.AlphaMissense_threshold.max(),
             ta.AlphaMissense_threshold.median()))
    print("\n[saved] r1_out/r1e_threshold_fix_impact.tsv")


if __name__ == "__main__":
    compare_threshold_functions()
