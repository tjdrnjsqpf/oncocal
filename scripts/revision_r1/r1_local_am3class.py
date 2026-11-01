"""R1 local sensitivity on fig2_master_3922.tsv (no server data needed).
Missed definitions: raw<0.5 (original) vs AM official non-LP (<0.564) vs likely-benign (<0.34).
Role sets: all 3,922 vs core roles (oncogene/TSG/both; excludes fusion_only & unknown, e.g. MUC6, HLA-A).
Tests: medians, Mann-Whitney, gene-clustered GEE (cluster = gene)."""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm, statsmodels.formula.api as smf
M = os.path.expanduser("~/oncocal/scripts/revision/fig2_master_3922.tsv")
d = pd.read_csv(M, sep="\t")
rows = []
for set_nm, sub in [("all", d), ("core_roles", d[d.role_class.isin(["oncogene", "TSG", "both"])])]:
    for def_nm, thr in [("raw<0.5", 0.5), ("notLP<0.564", 0.564), ("LB<0.34", 0.34)]:
        s = sub.copy(); s["missed"] = (s.am_pathogenicity < thr).astype(int)
        for ax, col in [("phyloP", "phyloP100way"), ("RSA", "struct_rsa")]:
            x = s[["gene", "missed", col]].dropna().rename(columns={col: "x"})
            mi, ca = x[x.missed == 1].x, x[x.missed == 0].x
            u, p = mannwhitneyu(mi, ca)
            x["z"] = (x.x - x.x.mean()) / x.x.std()
            g = smf.gee("missed ~ z", groups="gene", data=x, family=sm.families.Binomial(),
                        cov_struct=sm.cov_struct.Exchangeable()).fit()
            rows.append(dict(set=set_nm, missed_def=def_nm, axis=ax, n=len(x), n_missed=len(mi),
                             pct_missed=round(100 * len(mi) / len(x), 1),
                             med_missed=round(mi.median(), 3), med_caught=round(ca.median(), 3),
                             rank_biserial=round(1 - 2 * u / (len(mi) * len(ca)), 3), mwu_p=p,
                             gee_OR_perSD=round(np.exp(g.params["z"]), 3), gee_p=g.pvalues["z"],
                             n_gene=x.gene.nunique()))
out = pd.DataFrame(rows)
out.to_csv("r1_local_am3class.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(out.to_string(index=False))
# missed counts by role, 3 definitions
d["cls"] = pd.cut(d.am_pathogenicity, [-1, .34, .564, 2], right=False, labels=["LB", "ambiguous", "LP"])
t = pd.crosstab(d.role_class, d.cls, margins=True); t.to_csv("r1_am3class_by_role.tsv", sep="\t"); print(t)
