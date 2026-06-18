"""Figure 입력 export (새 분석 없음, 기존 파이프라인 재사용).
  [1] fig2_rsa_input.tsv      — R2 집단(oncogenic positive)의 AM_status + struct_rsa
  [2] oof_predictions.tsv     — gene-group 5-fold OOF(tools / tools+mech+struct) + AM raw + VARITY_R
사용: python -m src.models.export_figure_data
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

from .train_baseline import run_oof, _tool_feats
from ..features.mechanism_features import add_mechanism_features, MECH_COLS
from ..features.structure_features import add_structure_features, STRUCT_COLS

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "paper_handoff", "tables")


def main():
    os.makedirs(OUT, exist_ok=True)
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))

    # ---------- [1] fig2_rsa_input.tsv ----------
    pos = df[df.oncogenic_label == "positive"].drop_duplicates("pkey").copy()
    pos = add_structure_features(pos)
    pos["AM_status"] = np.where(pos.am_pathogenicity < 0.5, "AM missed", "AM caught")
    rsa = pos[["pkey", "AM_status", "struct_rsa"]].copy()
    rsa.to_csv(os.path.join(OUT, "fig2_rsa_input.tsv"), sep="\t", index=False)
    mm = rsa.loc[rsa.AM_status == "AM missed", "struct_rsa"].dropna()
    cc = rsa.loc[rsa.AM_status == "AM caught", "struct_rsa"].dropna()
    p = mannwhitneyu(mm, cc).pvalue
    print("[1] fig2_rsa_input.tsv  (n=%d, struct_rsa 결측 %d)" % (len(rsa), rsa.struct_rsa.isna().sum()))
    print("    AM missed: n=%d 중앙 RSA=%.3f | AM caught: n=%d 중앙 RSA=%.3f | Mann-Whitney p=%.2e"
          % (len(mm), mm.median(), len(cc), cc.median(), p))

    # ---------- [2] oof_predictions.tsv ----------
    split = pd.read_csv(os.path.join(ROOT, "data", "processed", "split_gene.tsv"), sep="\t")
    lab = df[df.y.notna()].drop_duplicates("pkey").merge(
        split[["pkey", "fold"]].drop_duplicates("pkey"), on="pkey", how="inner").copy()
    lab["y"] = lab["y"].astype(int)
    lab = add_structure_features(add_mechanism_features(lab))
    y, fold = lab.y.values, lab.fold.values
    tool_feats = _tool_feats(lab)
    oof_tools = run_oof(lab, tool_feats, y, fold)
    oof_full = run_oof(lab, tool_feats + MECH_COLS + STRUCT_COLS, y, fold)

    out = pd.DataFrame({
        "pkey": lab.pkey.values, "y_true": y, "role_class": lab.role_class.values,
        "p_model": np.round(oof_full, 6), "p_model_tools": np.round(oof_tools, 6),
        "AlphaMissense": lab.am_pathogenicity.values,
        "VARITY_R": lab.get("dbnsfp_VARITY_R_rankscore"),
    })
    out.to_csv(os.path.join(OUT, "oof_predictions.tsv"), sep="\t", index=False)
    au = lambda s, m: roc_auc_score(y[m], s[m]) if len(np.unique(y[m])) > 1 else float("nan")
    allm = np.ones(len(out), bool)
    print("\n[2] oof_predictions.tsv  (unique pkey n=%d, pos=%d)" % (len(out), int(y.sum())))
    print("    OOF AUROC(전체): p_model=%.4f, p_model_tools=%.4f, AlphaMissense=%.4f, VARITY_R=%.4f"
          % (au(oof_full, allm), au(oof_tools, allm),
             roc_auc_score(y, lab.am_pathogenicity.fillna(0.5)),
             roc_auc_score(y, lab["dbnsfp_VARITY_R_rankscore"].fillna(0.5))))
    print("    VARITY_R 결측: %d / AM 결측: %d"
          % (out.VARITY_R.isna().sum(), out.AlphaMissense.isna().sum()))
    print("\n[done] -> paper_handoff/tables/{fig2_rsa_input,oof_predictions}.tsv")


if __name__ == "__main__":
    main()
