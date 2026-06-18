"""GOF/LOF 기전 분류 — Tran 2025 미보유 차별점("기전을 함께 예측").

라벨: CIViC 변이-수준 기전(civic_mechanism: GOF / LOF / DN). gold지만 소량.
  - 이진: GOF(1) vs LOF(0). DN(dominant-negative)은 별도 보고(주로 TSG, LOF 계열).
  - ⚠️ 유전자 역할 교란: oncogene→GOF, TSG→LOF 가 강한 prior → 변이특이 feature가 그 너머를 예측하는지가 관건.
    **leave-one-gene-out(LOGO) CV** + gene-role-only 베이스라인 대비로 정직하게 평가.

산출: (1) GOF vs LOF 기전 시그니처(feature 차이, Mann-Whitney) (2) LOGO 분류 AUROC.
사용: python -m src.models.mechanism_classifier
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
    # dataset.parquet 에 이미 civic_mechanism 결합됨 → 직접 사용
    d = df.drop_duplicates("pkey")
    d = d[d["civic_mechanism"].notna()].copy()
    print("CIViC 기전 라벨 변이: %d" % len(d))
    print(d["civic_mechanism"].value_counts().to_string())

    d = add_structure_features(add_mechanism_features(d))
    bina = d[d.civic_mechanism.isin(["GOF", "LOF"])].copy()
    bina["gof"] = (bina.civic_mechanism == "GOF").astype(int)
    print("\n이진(GOF vs LOF): GOF=%d, LOF=%d, 유전자=%d"
          % (int(bina.gof.sum()), int((bina.gof == 0).sum()), bina.gene.nunique()))

    # ---- (1) 기전 시그니처 ----
    print("\n[기전 시그니처] GOF vs LOF 중앙값 (Mann-Whitney p):")
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

    # ---- (2) LOGO 분류 ----
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
    # gene-role 베이스라인: oncogene→GOF
    role = bina.role_class.values
    role_pred = np.where(role == "oncogene", 1.0, np.where(role == "TSG", 0.0, 0.5))
    auc_role = roc_auc_score(y, role_pred)
    print("\n[LOGO 분류] 변이특이 모델 AUROC=%.3f (n=%d) vs gene-role 베이스라인 AUROC=%.3f"
          % (auc_model, ok.sum(), auc_role))
    print("  → 모델이 role 베이스라인을 넘으면 '유전자 역할 너머의 변이특이 기전 신호' 입증.")
    print("\n[done] -> results/tables/gof_lof_signatures.tsv")


if __name__ == "__main__":
    main()
