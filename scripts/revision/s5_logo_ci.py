"""S5 — (1) leave-one-out per-gene threshold balanced-acc (remove in-sample bias),
        (2) OncoCal vs best tool AUROC bootstrap 95% CI.
Usage: python scripts/revision/s5_logo_ci.py
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from src.paths import find_table, revision_out_dir

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RNG = np.random.RandomState(2024)

def youden_threshold(y, s):
    order = np.argsort(s); ys, ss = y[order], s[order]
    P, N = ys.sum(), len(ys)-ys.sum()
    if P == 0 or N == 0: return np.nan
    best_j, best_t, tp, fp, prev = -1, np.nan, P, N, None
    for i in range(len(ss)):
        if prev is not None and ss[i] != prev:
            j = tp/P - fp/N
            if j > best_j: best_j, best_t = j, prev
        if ys[i] == 1: tp -= 1
        else: fp -= 1
        prev = ss[i]
    return best_t

# ---------- (1) leave-one-out (variant) per-gene threshold ----------
df = pd.read_csv(os.path.join(ROOT, "data", "processed", "baseline_oof.tsv"), sep="\t")
df = df.dropna(subset=["model_pred", "am_pathogenicity"])
gt_model = youden_threshold(df.y.values, df.model_pred.values)
gt_am = youden_threshold(df.y.values, df.am_pathogenicity.values)

def loo_pergene(score_col, gt):
    """Per gene: for each variant, compute Youden threshold on the rest (leaving that variant out) -> held-out prediction."""
    preds_global, preds_loo, ys = [], [], []
    for gene, g in df.groupby("gene"):
        if (g.y==1).sum() < 5 or (g.y==0).sum() < 5:
            continue
        y = g.y.values; s = g[score_col].values
        for i in range(len(g)):
            mask = np.ones(len(g), bool); mask[i] = False
            t = youden_threshold(y[mask], s[mask])
            if np.isnan(t): t = gt
            preds_loo.append(int(s[i] >= t))
            preds_global.append(int(s[i] >= gt))
            ys.append(y[i])
    ys = np.array(ys)
    return (balanced_accuracy_score(ys, preds_global),
            balanced_accuracy_score(ys, preds_loo), len(ys),
            len(set(df[( (df.groupby('gene').y.transform(lambda x:(x==1).sum())>=5) &
                         (df.groupby('gene').y.transform(lambda x:(x==0).sum())>=5))].gene)))

print("=== (1) per-gene threshold: in-sample vs leave-one-out(variant) ===")
print("Global Youden threshold: model=%.3f, AM=%.3f" % (gt_model, gt_am))
rows=[]
for col, gt, nm in [("am_pathogenicity", gt_am, "AlphaMissense"), ("model_pred", gt_model, "OncoCal")]:
    bg, bl, n, ng = loo_pergene(col, gt)
    print("  %-13s: balanced-acc global %.3f -> per-gene LOO %.3f (delta%+.3f) [n_variant=%d, gene=%d]"
          % (nm, bg, bl, bl-bg, n, ng))
    rows.append({"predictor":nm,"global_bacc":round(bg,3),"pergene_LOO_bacc":round(bl,3),
                 "delta":round(bl-bg,3),"n_variant":n,"n_gene":ng})
pd.DataFrame(rows).to_csv(os.path.join(revision_out_dir(),"s5_logo_bacc.tsv"),sep="\t",index=False)

# ---------- (2) AUROC bootstrap 95% CI ----------
oof = pd.read_csv(find_table("oof_predictions.tsv"), sep="\t")
print("\n=== (2) OncoCal vs best tool AUROC bootstrap 95%% CI (oof n=%d, pos=%d) ===" % (len(oof), int(oof.y_true.sum())))
y = oof.y_true.values
preds = {"OncoCal(tools)": oof.p_model_tools.values, "OncoCal(full)": oof.p_model_full.values,
         "AlphaMissense": oof.AlphaMissense.fillna(0.5).values,
         "VARITY_R": oof.VARITY_R.fillna(0.5).values}
def boot_auroc(yv, sv, B=2000):
    n=len(yv); aus=[]
    for _ in range(B):
        idx = RNG.randint(0, n, n)
        if len(np.unique(yv[idx]))<2: continue
        aus.append(roc_auc_score(yv[idx], sv[idx]))
    return np.percentile(aus,[2.5,50,97.5])
base = {}
for nm, s in preds.items():
    lo,md,hi = boot_auroc(y, s)
    base[nm]=s
    print("  %-15s AUROC=%.4f  95%%CI [%.4f, %.4f]" % (nm, roc_auc_score(y,s), lo, hi))
# gap CI: OncoCal(tools) - best single tool (paired bootstrap)
best_tool = "VARITY_R" if roc_auc_score(y,preds["VARITY_R"])>=roc_auc_score(y,preds["AlphaMissense"]) else "AlphaMissense"
print("\n  best single tool = %s" % best_tool)
n=len(y); gaps=[]
sm_model, sm_best = preds["OncoCal(tools)"], preds[best_tool]
for _ in range(5000):
    idx=RNG.randint(0,n,n)
    if len(np.unique(y[idx]))<2: continue
    gaps.append(roc_auc_score(y[idx],sm_model[idx])-roc_auc_score(y[idx],sm_best[idx]))
glo,gmd,ghi=np.percentile(gaps,[2.5,50,97.5])
print("  deltaAUROC OncoCal(tools) - %s = %.4f  95%%CI [%.4f, %.4f]  (P(delta>0)=%.3f)"
      % (best_tool, gmd, glo, ghi, np.mean(np.array(gaps)>0)))
pd.DataFrame([{"model":"OncoCal_tools","best_tool":best_tool,"delta_auroc":round(gmd,4),
              "ci_lo":round(glo,4),"ci_hi":round(ghi,4),"P_delta_gt0":round(np.mean(np.array(gaps)>0),4)}]
             ).to_csv(os.path.join(revision_out_dir(),"s5_auroc_ci.tsv"),sep="\t",index=False)
print("\n[saved] s5_logo_bacc.tsv, s5_auroc_ci.tsv")
