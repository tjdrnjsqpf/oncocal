"""S4 — ascertainment control: do missed vs caught phyloP/RSA differences persist after controlling for label source/recurrence count?
Set = S1 unified (unique oncogenic positive, dedup, missed=RAW<0.5).
(1) diagnose source skew, (2) gene-clustered GEE within a single source (hotspot-only, cosmic-only),
(3) tests within recurrence-matched strata.
Usage: python scripts/revision/s4_ascertainment.py
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm, statsmodels.formula.api as smf
from src.features.structure_features import add_structure_features

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
pos = df[df.oncogenic_label == "positive"].drop_duplicates("pkey").copy()
pos = add_structure_features(pos)
pos["missed"] = (pos["am_pathogenicity"] < 0.5).astype(int)
pos["src_hotspot"] = pos.is_hotspot.fillna(False).astype(bool)
pos["src_cosmic10"] = pos.cosmic_sample_count.fillna(0) >= 10
PHY = "dbnsfp_phyloP100way_vertebrate"

# (1) source skew: are missed concentrated in a particular source?
print("(1) missed rate by source (ascertainment skew diagnostic):")
for nm, m in [("hotspot", pos.src_hotspot), ("cosmic>=10", pos.src_cosmic10),
              ("hotspot-only", pos.src_hotspot & ~pos.src_cosmic10),
              ("cosmic-only", pos.src_cosmic10 & ~pos.src_hotspot),
              ("both", pos.src_hotspot & pos.src_cosmic10)]:
    s = pos[m]
    print("   %-14s n=%4d | missed %.1f%%" % (nm, len(s), 100*s.missed.mean()))

def gee_test(sub, col, label):
    d = sub[["gene", "missed", col]].dropna(subset=[col]).rename(columns={col: "x"})
    if d.missed.nunique() < 2 or d.gene.nunique() < 10:
        return "   [%s] insufficient n (n=%d, gene=%d)" % (label, len(d), d.gene.nunique())
    d["z"] = (d.x - d.x.mean())/d.x.std()
    miss = d.loc[d.missed==1,"x"]; caught = d.loc[d.missed==0,"x"]
    u,_ = mannwhitneyu(miss, caught, alternative="two-sided")
    r = 1 - 2*u/(len(miss)*len(caught))
    try:
        g = smf.gee("missed ~ z", groups="gene", data=d, family=sm.families.Binomial(),
                    cov_struct=sm.cov_struct.Exchangeable()).fit()
        orr = np.exp(g.params["z"]); ci = np.exp(g.conf_int().loc["z"]); p = g.pvalues["z"]
        return ("   [%s] n=%d gene=%d | missed med=%.3f caught med=%.3f | r=%.3f | "
                "GEE OR/SD=%.3f[%.3f,%.3f] p=%.2e"
                % (label, len(d), d.gene.nunique(), miss.median(), caught.median(), r,
                   orr, ci[0], ci[1], p))
    except Exception as e:
        return "   [%s] GEE failed: %s" % (label, e)

print("\n(2) gene-clustered GEE within single-source strata (persistence after control):")
for col, cname in [(PHY, "phyloP"), ("struct_rsa", "RSA")]:
    print("  -- %s --" % cname)
    print(gee_test(pos, col, "all"))
    print(gee_test(pos[pos.src_hotspot], col, "hotspot-positive only"))
    print(gee_test(pos[pos.src_cosmic10], col, "cosmic>=10-positive only"))
    print(gee_test(pos[pos.src_cosmic10 & ~pos.src_hotspot], col, "cosmic-only(non-hotspot)"))

# (3) recurrence-matched strata: tests within cosmic_sample_count bins
print("\n(3) tests within recurrence-matched strata (cosmic_sample_count quantile bins):")
pos["rec"] = pos.cosmic_sample_count.fillna(0)
pos["rec_bin"] = pd.qcut(pos.rec.rank(method="first"), 4, labels=["Q1","Q2","Q3","Q4"])
for col, cname in [(PHY, "phyloP"), ("struct_rsa", "RSA")]:
    print("  -- %s --" % cname)
    for b in ["Q1","Q2","Q3","Q4"]:
        print(gee_test(pos[pos.rec_bin==b], col, "rec %s"%b))
