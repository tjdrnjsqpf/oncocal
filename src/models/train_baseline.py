"""Phase 2 — tool stack vs tools+mechanism features (gene-group CV, compared to best single tool).

Goal: (RQ3) does the tool stack improve over the best single tool + (differentiator) **does adding mechanism features improve it further**.
  - tool feature = dbNSFP 49 *_rankscore (all label-defining columns excluded — block leakage).
  - mech feature = sequence biochemical changes (src/features/mechanism_features): cysteine-gain, hydrophobicity/volume/charge changes, etc.
  - training unit = **unique pkey** (dedup codon duplicates). **gene-group 5-fold OOF**.
  - comparison: tools-only vs tools+mech vs best single tool. rescue (recovery of AM-missed positives) change. SHAP.

Output: results/tables/model_vs_tools.tsv, data/processed/baseline_oof.tsv, results/tables/shap_top.tsv
Usage: python -m src.models.train_baseline
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.metrics import roc_auc_score, average_precision_score

from ..features.mechanism_features import add_mechanism_features, MECH_COLS
from ..features.structure_features import add_structure_features, STRUCT_COLS

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBL = os.path.join(ROOT, "results", "tables")
PROC = os.path.join(ROOT, "data", "processed")

PARAMS = dict(objective="binary", metric="auc", learning_rate=0.03, num_leaves=31,
              min_child_samples=40, subsample=0.8, subsample_freq=1, colsample_bytree=0.7,
              reg_lambda=1.0, n_estimators=2000, verbose=-1, random_state=42)


def _tool_feats(df):
    return [c for c in df.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]


def run_oof(lab, feature_cols, y, fold):
    X = lab[feature_cols].astype(float).values
    folds = sorted(np.unique(fold))
    oof = np.full(len(lab), np.nan)
    for i, f in enumerate(folds):
        vf = folds[(i + 1) % len(folds)]
        te, va = fold == f, fold == folds[(i + 1) % len(folds)]
        tr = ~te & ~va
        m = lgb.LGBMClassifier(**PARAMS)
        m.fit(X[tr], y[tr], eval_set=[(X[va], y[va])], eval_metric="auc",
              callbacks=[lgb.early_stopping(80, verbose=False), lgb.log_evaluation(0)])
        oof[te] = m.predict_proba(X[te])[:, 1]
    return oof


def _auc(y, p, mask):
    s = pd.DataFrame({"y": y[mask], "p": p[mask]}).dropna()
    if s["y"].nunique() < 2:
        return np.nan, np.nan
    return roc_auc_score(s.y, s.p), average_precision_score(s.y, s.p)


def main():
    os.makedirs(TBL, exist_ok=True)
    df = pd.read_parquet(os.path.join(PROC, "dataset.parquet"))
    split = pd.read_csv(os.path.join(PROC, "split_gene.tsv"), sep="\t")
    lab = df[df.y.notna()].drop_duplicates("pkey").merge(
        split[["pkey", "fold"]].drop_duplicates("pkey"), on="pkey", how="inner").copy()
    lab["y"] = lab["y"].astype(int)
    lab = add_mechanism_features(lab)
    lab = add_structure_features(lab)
    y, fold = lab.y.values, lab.fold.values
    tool_feats = _tool_feats(lab)
    print("labeled=%d (pos=%d) | tool=%d mech=%d struct=%d | struct missing rate=%.2f"
          % (len(lab), int(y.sum()), len(tool_feats), len(MECH_COLS), len(STRUCT_COLS),
             lab["struct_rsa"].isna().mean()))

    oof_tools = run_oof(lab, tool_feats, y, fold)
    oof_mech = run_oof(lab, tool_feats + MECH_COLS, y, fold)
    oof_full = run_oof(lab, tool_feats + MECH_COLS + STRUCT_COLS, y, fold)
    lab["model_pred"] = oof_full

    # ---- Hypothesis test: RSA/burial of AM-missed vs AM-caught positives ----
    from scipy.stats import mannwhitneyu
    pp = lab[lab.y == 1].copy()
    pp["am_missed"] = pp.am_pathogenicity < 0.5
    mm, cc = pp.loc[pp.am_missed, "struct_rsa"].dropna(), pp.loc[~pp.am_missed, "struct_rsa"].dropna()
    if len(mm) > 10 and len(cc) > 10:
        print("[hypothesis] positive RSA (surface exposure): AM missed median %.3f vs caught %.3f (Mann-Whitney p=%.1e)"
              % (mm.median(), cc.median(), mannwhitneyu(mm, cc).pvalue))

    groups = {"all": np.ones(len(lab), bool),
              "oncogene": (lab.role_class == "oncogene").values,
              "TSG": (lab.role_class == "TSG").values}
    # best single tool
    best_au, best_name = -1, None
    for c in tool_feats:
        au, _ = _auc(y, lab[c].values, groups["all"])
        if au and au > best_au:
            best_au, best_name = au, c
    rows = []
    for name, p in [("MODEL_tools", oof_tools), ("MODEL_tools+mech", oof_mech),
                    ("MODEL_tools+mech+struct", oof_full),
                    ("best_tool:" + best_name.replace("dbnsfp_", "").replace("_rankscore", ""), lab[best_name].values),
                    ("AlphaMissense", lab["dbnsfp_AlphaMissense_rankscore"].values)]:
        r = {"predictor": name}
        for g, mask in groups.items():
            au, ap = _auc(y, p, mask)
            r["AUROC_" + g], r["AUPRC_" + g] = au, ap
        rows.append(r)
    comp = pd.DataFrame(rows)
    comp.to_csv(os.path.join(TBL, "model_vs_tools.tsv"), sep="\t", index=False)
    print("\n[comparison] gene-group OOF:")
    print(comp[["predictor", "AUROC_all", "AUPRC_all", "AUROC_oncogene", "AUROC_TSG"]].round(4).to_string(index=False))
    a0, a1, a2 = comp.loc[0, "AUROC_all"], comp.loc[1, "AUROC_all"], comp.loc[2, "AUROC_all"]
    print("\n  +mech(sequence) effect: AUROC %+.4f | +struct additional effect: AUROC %+.4f" % (a1 - a0, a2 - a1))
    print("  final tools+mech+struct vs best tool(%.4f): AUROC %+.4f" % (best_au, a2 - best_au))
    print("  oncogene(GOF) AUROC: tools %.4f → +mech %.4f → +struct %.4f"
          % (comp.loc[0, "AUROC_oncogene"], comp.loc[1, "AUROC_oncogene"], comp.loc[2, "AUROC_oncogene"]))

    # rescue comparison
    pos = lab[lab.y == 1]
    am_missed = pos[pos.am_pathogenicity < 0.5]
    for tag, col in [("tools", oof_tools), ("tools+mech", oof_mech), ("tools+mech+struct", oof_full)]:
        lab["_p"] = col
        r = am_missed.merge(lab[["pkey", "_p"]], on="pkey")
        print("  [rescue/%s] of %d AM-missed positives, rescued(pred>=0.5)=%d (%.1f%%)"
              % (tag, len(am_missed), int((r["_p"] >= 0.5).sum()), 100 * (r["_p"] >= 0.5).mean()))
    for k in ["JAK2:p.V617F", "MPL:p.W515L", "FGFR3:p.Y373C"]:
        if k in set(lab.pkey):
            row = lab[lab.pkey == k].iloc[0]
            i = lab.index.get_loc(row.name)
            print("    %-16s AM=%.3f | tools=%.3f | +mech=%.3f | +struct=%.3f (RSA=%.2f)"
                  % (k, row.am_pathogenicity, oof_tools[i], oof_mech[i], oof_full[i],
                     row.struct_rsa if row.struct_rsa == row.struct_rsa else -1))

    lab[["pkey", "gene", "role_class", "y", "fold", "model_pred", "am_pathogenicity",
         "cosmic_sample_count"]].to_csv(os.path.join(PROC, "baseline_oof.tsv"), sep="\t", index=False)

    # SHAP (retrain on the full tools+mech set)
    try:
        import shap
        allf = tool_feats + MECH_COLS + STRUCT_COLS
        mfull = lgb.LGBMClassifier(**{**PARAMS, "n_estimators": 300}).fit(lab[allf].astype(float).values, y)
        sv = shap.TreeExplainer(mfull).shap_values(lab[allf].astype(float).values)
        sv = sv[1] if isinstance(sv, list) else sv
        imp = pd.DataFrame({"feature": [f.replace("dbnsfp_", "").replace("_rankscore", "") for f in allf],
                            "is_mech": [f in MECH_COLS or f in STRUCT_COLS for f in allf],
                            "mean_abs_shap": np.abs(sv).mean(0)}).sort_values("mean_abs_shap", ascending=False)
        imp.to_csv(os.path.join(TBL, "shap_top.tsv"), sep="\t", index=False)
        print("\n[SHAP] top 15 (mech=mechanism/structure features):")
        print(imp.head(15).to_string(index=False))
        print("  top mechanism/structure features:")
        print(imp[imp.is_mech].head(8).to_string(index=False))
    except Exception as e:
        print("[SHAP] skip:", e)
    print("\n[done] -> results/tables/{model_vs_tools,shap_top}.tsv, data/processed/baseline_oof.tsv")


if __name__ == "__main__":
    main()
