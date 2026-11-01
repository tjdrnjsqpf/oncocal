"""R1-C — Matched-negative sensitivity analysis (reviewer #1 comment 3).

Reviewer: "there is also an unavoidable asymmetry between positives and negatives:
positives are largely somatic cancer observations, while negatives come from ClinVar
benign/likely-benign and common gnomAD germline variants. The author correctly
acknowledges this, but it can itself contribute to the high discrimination of OncoCal.
I would emphasize this earlier and, if possible, add a sensitivity analysis using more
closely matched negative variants by gene and allele-frequency/context."

Design
------
Evaluation-side only: we re-score the EXISTING out-of-fold predictions on progressively
harder negative sets. No retraining -- the model still saw the original training
distribution, so this bounds how much of the AUROC is explained by the easy negatives
rather than removing the asymmetry from training. State that limitation in the paper.

Negative sets (positives unchanged):
  original         all negatives (manuscript definition)
  clinvar_only     ClinVar benign/likely-benign only (drop gnomAD-common-only negatives)
  rare_benign      ClinVar benign AND NOT gnomAD-common  -> removes the allele-frequency
                   shortcut: both arms are now rare variants
  gene_matched     genes carrying both arms; per gene, negatives subsampled 1:1 to
                   positives -> removes gene identity as a shortcut
  gene_matched_rare  rare_benign AND gene-matched (hardest)

Subsampled sets are repeated N_REPEAT times with different seeds; we report mean +/- sd
AUROC plus a gene-cluster bootstrap CI on the final stratum.

Run:  python r1c_matched_negatives.py
Out:  r1_out/r1c_matched_negatives.tsv
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

N_REPEAT = 20
AF_COL = "dbnsfp_gnomAD4.1_joint_AF"
SCORES = ["p_model_tools", "p_model_full", "AlphaMissense", "VARITY_R"]


def auroc(y, s):
    m = ~pd.isna(s)
    if len(np.unique(y[m])) < 2:
        return np.nan
    return roc_auc_score(y[m], s[m])


def main():
    oof = pd.read_csv(find_table("oof_predictions.tsv"), sep="\t")
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    keep = ["pkey", "gene", "clinvar_benign", "cosmic_sample_count"]
    if AF_COL in df.columns:
        keep.append(AF_COL)
    if "gnomad_common" in df.columns:
        keep.append("gnomad_common")
    meta = df[keep].drop_duplicates("pkey")
    d = oof.merge(meta, on="pkey", how="left")
    d["clinvar_benign"] = d.clinvar_benign.fillna(False).astype(bool)
    if "gnomad_common" in d.columns:
        d["gnomad_common"] = d.gnomad_common.fillna(False).astype(bool)
    else:
        d["gnomad_common"] = pd.to_numeric(d.get(AF_COL), errors="coerce").fillna(0) >= 1e-3
    print("oof n=%d (pos=%d, neg=%d) | genes=%d"
          % (len(d), int((d.y_true == 1).sum()), int((d.y_true == 0).sum()), d.gene.nunique()))

    pos = d[d.y_true == 1]
    negsets = {
        "original": d[d.y_true == 0],
        "clinvar_only": d[(d.y_true == 0) & d.clinvar_benign],
        "rare_benign": d[(d.y_true == 0) & d.clinvar_benign & ~d.gnomad_common],
    }
    for nm, n in negsets.items():
        print("  negatives[%-13s] n=%5d" % (nm, len(n)))

    rows = []

    def record(stratum, frame, rep=np.nan):
        r = {"stratum": stratum, "repeat": rep, "n_pos": int((frame.y_true == 1).sum()),
             "n_neg": int((frame.y_true == 0).sum()), "n_gene": frame.gene.nunique()}
        y = frame.y_true.values
        for s in SCORES:
            if s in frame.columns:
                r["auroc_%s" % s] = auroc(y, frame[s].values)
        rows.append(r)
        return r

    # ---- unmatched strata ----
    for nm, neg in negsets.items():
        r = record(nm, pd.concat([pos, neg]))
        print("\n[%s] pos=%d neg=%d | OncoCal %.4f | AlphaMissense %.4f"
              % (nm, r["n_pos"], r["n_neg"], r.get("auroc_p_model_tools", np.nan),
                 r.get("auroc_AlphaMissense", np.nan)))

    # ---- gene-matched strata (subsampled, repeated) ----
    for nm, neg in [("gene_matched", negsets["original"]),
                    ("gene_matched_rare", negsets["rare_benign"])]:
        shared = sorted(set(pos.gene.dropna()) & set(neg.gene.dropna()))
        print("\n[%s] genes with both arms: %d" % (nm, len(shared)))
        for rep in range(N_REPEAT):
            rng = np.random.RandomState(1000 + rep)
            parts = []
            for g in shared:
                p_g, n_g = pos[pos.gene == g], neg[neg.gene == g]
                k = min(len(p_g), len(n_g))
                if k == 0:
                    continue
                parts.append(p_g.iloc[rng.choice(len(p_g), k, replace=False)])
                parts.append(n_g.iloc[rng.choice(len(n_g), k, replace=False)])
            if parts:
                record(nm, pd.concat(parts), rep)
        sub = pd.DataFrame([r for r in rows if r["stratum"] == nm])
        for s in SCORES:
            col = "auroc_%s" % s
            if col in sub.columns and sub[col].notna().any():
                print("   %-14s AUROC %.4f +/- %.4f (n_pos=%d, n_neg=%d, %d repeats)"
                      % (s, sub[col].mean(), sub[col].std(), sub.n_pos.mean(),
                         sub.n_neg.mean(), len(sub)))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "r1c_matched_negatives.tsv"), sep="\t", index=False)

    agg = (res.groupby("stratum")
           .agg(n_pos=("n_pos", "mean"), n_neg=("n_neg", "mean"), n_gene=("n_gene", "mean"),
                **{("%s_mean" % s): ("auroc_%s" % s, "mean") for s in SCORES
                   if ("auroc_%s" % s) in res.columns},
                **{("%s_sd" % s): ("auroc_%s" % s, "std") for s in SCORES
                   if ("auroc_%s" % s) in res.columns})
           .round(4))
    agg.to_csv(os.path.join(OUT, "r1c_matched_negatives_summary.tsv"), sep="\t")
    print("\n--- summary ---")
    pd.set_option("display.width", 250)
    print(agg.to_string())
    print("\n[saved] r1_out/r1c_matched_negatives{,_summary}.tsv")
    print("NOTE for the manuscript: evaluation-side only; the model was trained on the "
          "original distribution, so these numbers bound (not remove) the asymmetry.")


if __name__ == "__main__":
    main()
