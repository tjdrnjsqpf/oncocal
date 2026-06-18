"""Fig 4 regeneration — consistent output over the same pipeline across 10 seeds (0..9).
  [1] model_vs_tools.tsv : AUROC/AUPRC mean±sd per predictor (all, oncogene, TSG). Tool baselines have sd=0.
  [2] oof_predictions.tsv: per-variant predictions from the representative seed (median full AUROC).
  [3] rescue mean±sd + JAK2 V617F.
gene-group 5-fold OOF, unique pkey, recurrence feature excluded (only tool rankscore+mech+struct).
Usage: python -m src.models.fig4_multiseed
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
OUTS = [os.path.join(ROOT, "paper_handoff", "tables"), os.path.join(ROOT, "results", "tables")]
PARAMS = dict(objective="binary", metric="auc", learning_rate=0.03, num_leaves=31,
              min_child_samples=40, subsample=0.8, subsample_freq=1, colsample_bytree=0.7,
              reg_lambda=1.0, n_estimators=2000, verbose=-1)
SEEDS = list(range(10))
GROUPS = ["all", "oncogene", "TSG"]


def run_oof(X, y, fold, seed):
    folds = sorted(np.unique(fold))
    oof = np.full(len(y), np.nan)
    for i, f in enumerate(folds):
        te = fold == f
        va = fold == folds[(i + 1) % len(folds)]
        tr = ~te & ~va
        m = lgb.LGBMClassifier(random_state=seed, **PARAMS)
        m.fit(X[tr], y[tr], eval_set=[(X[va], y[va])], eval_metric="auc",
              callbacks=[lgb.early_stopping(80, verbose=False), lgb.log_evaluation(0)])
        oof[te] = m.predict_proba(X[te])[:, 1]
    return oof


def auc_ap(score, y, mask):
    s = pd.DataFrame({"y": y[mask], "p": np.asarray(score)[mask]}).dropna()
    if s.y.nunique() < 2:
        return np.nan, np.nan
    return roc_auc_score(s.y, s.p), average_precision_score(s.y, s.p)


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    split = pd.read_csv(os.path.join(ROOT, "data", "processed", "split_gene.tsv"), sep="\t")
    lab = df[df.y.notna()].drop_duplicates("pkey").merge(
        split[["pkey", "fold"]].drop_duplicates("pkey"), on="pkey", how="inner").copy()
    lab["y"] = lab.y.astype(int)
    lab = add_structure_features(add_mechanism_features(lab))
    y, fold = lab.y.values, lab.fold.values
    tool_feats = [c for c in lab.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]
    FSETS = {"tools": tool_feats, "mech": tool_feats + MECH_COLS, "full": tool_feats + MECH_COLS + STRUCT_COLS}
    masks = {"all": np.ones(len(lab), bool),
             "oncogene": (lab.role_class == "oncogene").values,
             "TSG": (lab.role_class == "TSG").values}
    Xs = {k: lab[v].astype(float).values for k, v in FSETS.items()}
    print("labeled=%d pos=%d | feats tools=%d mech=%d full=%d | recurrence excluded check=%s"
          % (len(lab), int(y.sum()), len(tool_feats), len(FSETS["mech"]), len(FSETS["full"]),
             "cosmic_sample_count" not in FSETS["full"]))

    # OOF + metric per seed
    rec = {m: {g: {"AUROC": [], "AUPRC": []} for g in GROUPS} for m in FSETS}
    oofs = {m: {} for m in FSETS}   # model -> seed -> oof array
    for seed in SEEDS:
        for m in FSETS:
            oof = run_oof(Xs[m], y, fold, seed)
            oofs[m][seed] = oof
            for g in GROUPS:
                au, ap = auc_ap(oof, y, masks[g])
                rec[m][g]["AUROC"].append(au); rec[m][g]["AUPRC"].append(ap)
        print("  seed %d done (full AUROC_all=%.4f)" % (seed, rec["full"]["all"]["AUROC"][-1]), flush=True)

    # Tool baselines (seed-independent)
    base = {"best_tool:VARITY_R": lab["dbnsfp_VARITY_R_rankscore"].values,
            "AlphaMissense": lab["am_pathogenicity"].values}

    # [1] model_vs_tools.tsv (mean±sd)
    rows = []
    name_map = {"tools": "MODEL_tools", "mech": "MODEL_tools+mech", "full": "MODEL_tools+mech+struct"}
    for m in ["tools", "mech", "full"]:
        r = {"predictor": name_map[m]}
        for g in GROUPS:
            for met in ["AUROC", "AUPRC"]:
                v = np.array(rec[m][g][met])
                r["%s_%s_mean" % (met, g)] = round(np.nanmean(v), 4)
                r["%s_%s_sd" % (met, g)] = round(np.nanstd(v, ddof=1), 4)
        rows.append(r)
    for nm, sc in base.items():
        r = {"predictor": nm}
        for g in GROUPS:
            au, ap = auc_ap(sc, y, masks[g])
            r["AUROC_%s_mean" % g] = round(au, 4); r["AUROC_%s_sd" % g] = 0.0
            r["AUPRC_%s_mean" % g] = round(ap, 4); r["AUPRC_%s_sd" % g] = 0.0
        rows.append(r)
    mvt = pd.DataFrame(rows)

    # [2] Representative seed = seed closest to the median full AUROC_all
    fa = np.array(rec["full"]["all"]["AUROC"])
    med = np.median(fa)
    rep = SEEDS[int(np.argmin(np.abs(fa - med)))]
    oofp = pd.DataFrame({
        "pkey": lab.pkey.values, "y_true": y, "role_class": lab.role_class.values,
        "p_model_tools": np.round(oofs["tools"][rep], 6),
        "p_model_mech": np.round(oofs["mech"][rep], 6),
        "p_model_full": np.round(oofs["full"][rep], 6),
        "AlphaMissense": lab.am_pathogenicity.values,
        "VARITY_R": lab.get("dbnsfp_VARITY_R_rankscore"),
    })

    for d in OUTS:
        os.makedirs(d, exist_ok=True)
        mvt.to_csv(os.path.join(d, "model_vs_tools.tsv"), sep="\t", index=False)
        oofp.to_csv(os.path.join(d, "oof_predictions.tsv"), sep="\t", index=False)

    # ---- Console summary ----
    print("\n=== [1] model_vs_tools (10 seeds mean±sd, all) ===")
    for _, r in mvt.iterrows():
        print("  %-26s AUROC_all=%.4f±%.4f  AUPRC_all=%.4f±%.4f  AUROC_onc=%.4f  AUROC_TSG=%.4f"
              % (r.predictor, r.AUROC_all_mean, r.AUROC_all_sd, r.AUPRC_all_mean, r.AUPRC_all_sd,
                 r.AUROC_oncogene_mean, r.AUROC_TSG_mean))

    print("\n=== seed noise vs feature effect ===")
    dt = np.array(rec["full"]["all"]["AUROC"]) - np.array(rec["tools"]["all"]["AUROC"])
    print("  full−tools AUROC_all: mean %+.4f, seed sd(per model)=tools±%.4f/full±%.4f"
          % (dt.mean(), np.nanstd(rec["tools"]["all"]["AUROC"], ddof=1),
             np.nanstd(rec["full"]["all"]["AUROC"], ddof=1)))
    print("  → |full−tools| (%.4f) vs seed sd range → decide whether the mech/struct effect is within noise" % abs(dt.mean()))

    print("\n=== [2] representative seed = %d (full AUROC_all=%.4f, median=%.4f) ===" % (rep, fa[rep], med))
    for m in ["tools", "mech", "full"]:
        au, ap = auc_ap(oofs[m][rep], y, masks["all"])
        print("    p_model_%s: AUROC=%.4f AUPRC=%.4f" % (m, au, ap))
    for nm, sc in base.items():
        au, _ = auc_ap(sc, y, masks["all"]); print("    %s: AUROC=%.4f" % (nm, au))

    # [3] rescue (10 seeds)
    pos = lab.y.values == 1
    am_missed = pos & (lab.am_pathogenicity.values < 0.5)
    nmiss = int(am_missed.sum())
    resc_full, resc_tools, jak2 = [], [], []
    ji = lab.index[lab.pkey == "JAK2:p.V617F"]
    jpos = lab.reset_index().index[lab.reset_index().pkey == "JAK2:p.V617F"]
    for seed in SEEDS:
        resc_full.append(100 * np.mean(oofs["full"][seed][am_missed] >= 0.5))
        resc_tools.append(100 * np.mean(oofs["tools"][seed][am_missed] >= 0.5))
        if len(jpos):
            jak2.append(oofs["full"][seed][jpos[0]])
    print("\n=== [3] rescue (10 seeds mean±sd) — AM raw<0.5 oncogenic positive n=%d ===" % nmiss)
    print("  full rescue rate: %.1f%% ± %.1f%%" % (np.mean(resc_full), np.std(resc_full, ddof=1)))
    print("  tools rescue rate: %.1f%% ± %.1f%%" % (np.mean(resc_tools), np.std(resc_tools, ddof=1)))
    if jak2:
        print("  JAK2 V617F p_model_full: %.3f ± %.3f (AM raw=0.334)" % (np.mean(jak2), np.std(jak2, ddof=1)))
    print("\n[done] -> paper_handoff/tables/ + results/tables/ {model_vs_tools, oof_predictions}.tsv")


if __name__ == "__main__":
    main()
