"""S3 — unsupervised vs supervised: threshold-free re-validation + 17.8/14.3 reconciliation.
(A) Explain 17.8 vs 14.3: recompute group means of driver_miss_rate.tsv by recurrence cut and classification (type vs type2).
(B) threshold-free: compare per-tool oncogene AUROC and oncogene-TSG gap (GOF penalty) by group (Mann-Whitney + permutation).
Usage: python scripts/revision/s3_threshold_free.py
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RNG = np.random.RandomState(999)
from src.paths import find_table, revision_out_dir

# ---------- (A) 17.8 vs 14.3 reconciliation ----------
mr = pd.read_csv(find_table("driver_miss_rate.tsv"), sep="\t")
print("driver_miss_rate.tsv: %d tools, columns=%s" % (len(mr), list(mr.columns)))
print("\n(A) group mean miss_rate(rankscore<0.5, rec>=50 driver n=%d):" % mr.n_drivers.max())
for grpcol in ["type", "type2"]:
    print("  -- classification=%s --" % grpcol)
    for g, sub in mr.groupby(grpcol):
        print("     %-12s n=%2d  mean miss=%.4f (%.1f%%)" % (g, len(sub), sub["miss_rate(<0.5)"].mean(),
                                                            100*sub["miss_rate(<0.5)"].mean()))
# unsup vs sup gap (type2 = final principle-based classification)
u = mr[mr.type2=="unsupervised"]["miss_rate(<0.5)"]; s = mr[mr.type2=="supervised"]["miss_rate(<0.5)"]
umean, smean = 100*u.mean(), 100*s.mean()
print("\n  > type2 final: unsup(n=%d) %.1f%% vs sup(n=%d) %.1f%% -> gap %.1f%%p"
      % (len(u), umean, len(s), smean, umean-smean))
mwp = mannwhitneyu(u, s, alternative="greater").pvalue
print("    Mann-Whitney(unsup>sup) p=%.4f" % mwp)

# ---------- (B) threshold-free: oncogene AUROC + GOF gap ----------
perf = pd.read_csv(find_table("tool_performance.tsv"), sep="\t")
og = perf[perf.group=="oncogene"].set_index("tool")["AUROC"].rename("onco_auroc")
ts = perf[perf.group=="TSG"].set_index("tool")["AUROC"].rename("tsg_auroc")
alla = perf[perf.group=="all"].set_index("tool")["AUROC"].rename("all_auroc")
tf = pd.concat([og, ts, alla], axis=1).dropna(subset=["onco_auroc","tsg_auroc"])
tf["gof_gap"] = tf["tsg_auroc"] - tf["onco_auroc"]  # positive = worse at catching oncogene(GOF)

# type2 mapping (consistent with tool names in driver_miss_rate)
type2 = mr.set_index("tool")["type2"]
tf["type2"] = tf.index.map(type2)
tf = tf.dropna(subset=["type2"])
print("\n(B) threshold-free — %d matched tools (unsup %d / sup %d)"
      % (len(tf), (tf.type2=="unsupervised").sum(), (tf.type2=="supervised").sum()))

def grp_compare(metric, greater):
    a = tf[tf.type2=="unsupervised"][metric]; b = tf[tf.type2=="supervised"][metric]
    p_mw = mannwhitneyu(a, b, alternative=greater).pvalue
    # permutation (label shuffle, gene-free: tool is the unit)
    obs = a.mean()-b.mean(); pooled = tf[metric].values; na=len(a); B=20000
    cnt=0
    for _ in range(B):
        idx = RNG.permutation(len(pooled)); dif = pooled[idx[:na]].mean()-pooled[idx[na:]].mean()
        if (dif>=obs) if greater=="greater" else (dif<=obs): cnt+=1
    p_perm=(cnt+1)/(B+1)
    return a.mean(), b.mean(), obs, p_mw, p_perm

print("\n  per-metric group comparison (unsup vs sup):")
rows=[]
for metric, greater, desc in [("gof_gap","greater","GOF penalty(TSG-onco AUROC), larger positive = more GOF underestimation"),
                              ("onco_auroc","less","oncogene AUROC, smaller = worse at catching GOF"),
                              ("all_auroc","less","overall AUROC")]:
    am,bm,obs,pmw,pperm = grp_compare(metric, greater)
    print("   [%s] unsup=%.4f sup=%.4f diff=%+.4f | MW p=%.4f perm p=%.4f\n      (%s)"
          % (metric, am, bm, obs, pmw, pperm, desc))
    rows.append({"metric":metric,"unsup_mean":round(am,4),"sup_mean":round(bm,4),
                 "diff":round(obs,4),"mannwhitney_p":round(pmw,4),"perm_p":round(pperm,4)})

pd.DataFrame(rows).to_csv(os.path.join(revision_out_dir(),"s3_threshold_free.tsv"),
                          sep="\t", index=False)
tf.reset_index().rename(columns={"index":"tool"}).to_csv(
    os.path.join(revision_out_dir(),"s3_tool_thresholdfree.tsv"), sep="\t", index=False)
print("\n[saved] s3_threshold_free.tsv, s3_tool_thresholdfree.tsv")
