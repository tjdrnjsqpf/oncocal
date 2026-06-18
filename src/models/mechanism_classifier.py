"""GOF/LOF mechanism classification — a differentiator Tran 2025 lacks ("predict the mechanism too").

Labels: CIViC variant-level mechanism (civic_mechanism: GOF / LOF / DN). Gold but scarce.
  - Binary: GOF(1) vs LOF(0). DN (dominant-negative) reported separately (mostly TSG, LOF family).
  - ⚠️ Gene-role confounder: oncogene→GOF, TSG→LOF is a strong prior → the question is whether variant-specific features predict beyond it.
    Evaluate honestly with **leave-one-gene-out (LOGO) CV** + comparison to a gene-role-only baseline.

Output: (1) GOF vs LOF mechanism signatures (feature differences, Mann-Whitney) (2) LOGO classification AUROC.
Usage: python -m src.models.mechanism_classifier
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd
import lightgbm as lgb
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut

from ..features.mechanism_features import add_mechanism_features, MECH_COLS
from ..features.structure_features import add_structure_features, STRUCT_COLS

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBL = os.path.join(ROOT, "results", "tables")


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    # civic_mechanism is already joined into dataset.parquet → use directly
    d = df.drop_duplicates("pkey")
    d = d[d["civic_mechanism"].notna()].copy()
    print("CIViC mechanism-labeled variants: %d" % len(d))
    print(d["civic_mechanism"].value_counts().to_string())

    d = add_structure_features(add_mechanism_features(d))
    bina = d[d.civic_mechanism.isin(["GOF", "LOF"])].copy()
    bina["gof"] = (bina.civic_mechanism == "GOF").astype(int)
    print("\nBinary (GOF vs LOF): GOF=%d, LOF=%d, genes=%d"
          % (int(bina.gof.sum()), int((bina.gof == 0).sum()), bina.gene.nunique()))

    # ---- (1) Mechanism signatures ----
    print("\n[mechanism signatures] GOF vs LOF medians (Mann-Whitney p):")
    sig_feats = ["struct_rsa", "struct_buried", "dbnsfp_phyloP100way_vertebrate", "dbnsfp_GERPpp_RS",
                 "mech_d_charge", "mech_abs_d_volume", "am_pathogenicity", "dbnsfp_AlphaMissense_rankscore"]
    rows = []
    for f in sig_feats:
        gof = bina.loc[bina.gof == 1, f].dropna()
        lof = bina.loc[bina.gof == 0, f].dropna()
        if len(gof) < 5 or len(lof) < 5:
            continue
        p = mannwhitneyu(gof, lof).pvalue
        rows.append({"feature": f.replace("dbnsfp_", ""), "GOF_median": round(gof.median(), 3),
                     "LOF_median": round(lof.median(), 3), "p": p})
    sig = pd.DataFrame(rows)
    sig.to_csv(os.path.join(TBL, "gof_lof_signatures.tsv"), sep="\t", index=False)
    print(sig.to_string(index=False))

    # ---- (2) LOGO classification ----
    feats = [c for c in d.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")] + MECH_COLS + STRUCT_COLS
    X, y, grp = bina[feats].astype(float).values, bina.gof.values, bina.gene.values
    logo = LeaveOneGroupOut()
    oof = np.full(len(bina), np.nan)
    for tr, te in logo.split(X, y, grp):
        if len(np.unique(y[tr])) < 2:
            continue
        m = lgb.LGBMClassifier(objective="binary", n_estimators=200, learning_rate=0.05,
                               num_leaves=15, min_child_samples=10, reg_lambda=1.0, verbose=-1)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    ok = ~np.isnan(oof)
    auc_model = roc_auc_score(y[ok], oof[ok]) if len(np.unique(y[ok])) > 1 else np.nan
    # gene-role baseline: oncogene→GOF
    role = bina.role_class.values
    role_pred = np.where(role == "oncogene", 1.0, np.where(role == "TSG", 0.0, 0.5))
    auc_role = roc_auc_score(y, role_pred)
    print("\n[LOGO classification] variant-specific model AUROC=%.3f (n=%d) vs gene-role baseline AUROC=%.3f"
          % (auc_model, ok.sum(), auc_role))
    print("  → If the model beats the role baseline, it demonstrates 'variant-specific mechanism signal beyond gene role'.")
    print("\n[done] -> results/tables/gof_lof_signatures.tsv")


if __name__ == "__main__":
    main()
