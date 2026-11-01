"""S2 — Recompute phyloP/RSA statistics (remove pseudo-replication).
Set = S1 unified definition (unique oncogenic positive, dedup pkey), missed = RAW<0.5.
Three methods: (1) gene-collapsed paired Wilcoxon (gene = unit)
        (2) GEE cluster-robust logistic (cluster = gene)
        (3) gene-cluster bootstrap 95% CI (rank-biserial r, median diff)
        (+) mixed-effects logistic (random intercept = gene) as auxiliary
Usage: python scripts/revision/s2_clustered.py
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu
import statsmodels.api as sm
import statsmodels.formula.api as smf
from src.features.structure_features import add_structure_features
from src.paths import revision_out_dir

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RNG = np.random.RandomState(12345)  # fixed seed instead of Math.random
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
pos = df[df.oncogenic_label == "positive"].drop_duplicates("pkey").copy()
pos = add_structure_features(pos)
pos["missed"] = (pos["am_pathogenicity"] < 0.5).astype(int)

AXES = [("phyloP", "dbnsfp_phyloP100way_vertebrate"), ("RSA", "struct_rsa")]
out_lines = []

def rank_biserial(a, b):
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return 1 - 2*u/(len(a)*len(b)), p

for name, col in AXES:
    d = pos[["gene", "missed", col]].dropna(subset=[col]).rename(columns={col: "x"})
    miss = d.loc[d.missed == 1, "x"]; caught = d.loc[d.missed == 0, "x"]
    r_naive, p_naive = rank_biserial(miss, caught)
    print("\n" + "="*66)
    print("[%s] naive: missed n=%d med=%.3f | caught n=%d med=%.3f | r=%.3f p=%.2e"
          % (name, len(miss), miss.median(), len(caught), caught.median(), r_naive, p_naive))

    # ---- (1) gene-collapsed paired Wilcoxon: within-gene missed median vs caught median ----
    g = d.groupby(["gene", "missed"])["x"].median().unstack("missed")
    both = g.dropna()  # genes having both missed(1) and caught(0)
    both.columns = ["caught_med", "missed_med"]
    if len(both) >= 6:
        w, pw = wilcoxon(both.missed_med, both.caught_med)
        diff = (both.missed_med - both.caught_med)
        print("  (1) gene-collapsed paired Wilcoxon: n_gene(both)=%d | median within-gene diff=%.3f | p=%.2e"
              % (len(both), diff.median(), pw))
    else:
        pw, both = float("nan"), both
        print("  (1) insufficient both-group genes (%d)" % len(both))
    # gene-level simple comparison (gene = unit, missed-dominant vs caught-dominant)
    gmed = d.groupby("gene").agg(missed_frac=("missed", "mean"), xmed=("x", "median"))

    # ---- (2) GEE cluster-robust logistic: missed ~ z(x), cluster=gene ----
    d2 = d.copy(); d2["z"] = (d2.x - d2.x.mean())/d2.x.std()
    gee = smf.gee("missed ~ z", groups="gene", data=d2,
                  family=sm.families.Binomial(),
                  cov_struct=sm.cov_struct.Exchangeable()).fit()
    coef = gee.params["z"]; se = gee.bse["z"]; pgee = gee.pvalues["z"]
    orr = np.exp(coef); ci = np.exp(gee.conf_int().loc["z"])
    print("  (2) GEE cluster-robust: OR/SD(z)=%.3f [%.3f,%.3f] | Wald p=%.2e | n_gene=%d"
          % (orr, ci[0], ci[1], pgee, d2.gene.nunique()))

    # ---- (3) gene-cluster bootstrap 95% CI (rank-biserial r, median diff) ----
    genes = d.gene.unique(); B = 2000
    rs, mds = [], []
    gd = {gname: sub for gname, sub in d.groupby("gene")}
    for _ in range(B):
        samp = pd.concat([gd[g] for g in RNG.choice(genes, len(genes), replace=True)])
        mm = samp.loc[samp.missed == 1, "x"]; cc = samp.loc[samp.missed == 0, "x"]
        if len(mm) < 5 or len(cc) < 5:
            continue
        u, _ = mannwhitneyu(mm, cc, alternative="two-sided")
        rs.append(1 - 2*u/(len(mm)*len(cc)))
        mds.append(cc.median() - mm.median())
    rlo, rhi = np.percentile(rs, [2.5, 97.5]); mlo, mhi = np.percentile(mds, [2.5, 97.5])
    print("  (3) gene-cluster bootstrap(B=%d): rank-biserial r=%.3f [%.3f,%.3f] | median diff(caught-missed) [%.3f,%.3f]"
          % (len(rs), np.median(rs), rlo, rhi, mlo, mhi))

    # ---- (+) mixed-effects logistic (random intercept=gene), variational Bayes ----
    try:
        vc = {"gene": "0 + C(gene)"}
        bm = sm.BinomialBayesMixedGLM.from_formula("missed ~ z", vc, d2).fit_vb()
        # fixed effect z: posterior mean / sd
        idx = list(bm.model.exog_names).index("z")
        fe, fsd = bm.fe_mean[idx], bm.fe_sd[idx]
        zval = fe/fsd
        print("  (+) mixed logistic(VB, RE=gene): z coef=%.3f (sd=%.3f, |z|=%.2f) -> %s"
              % (fe, fsd, abs(zval), "significant" if abs(zval) > 1.96 else "not significant"))
    except Exception as e:
        print("  (+) mixed logistic failed: %s" % e)

    out_lines.append({"axis": name, "naive_r": round(r_naive,3), "naive_p": p_naive,
                      "gee_OR_per_SD": round(orr,3), "gee_CI_lo": round(ci[0],3),
                      "gee_CI_hi": round(ci[1],3), "gee_wald_p": pgee,
                      "wilcoxon_gene_p": pw, "n_gene_both": len(both) if isinstance(both, pd.DataFrame) else np.nan,
                      "boot_r_med": round(np.median(rs),3), "boot_r_lo": round(rlo,3), "boot_r_hi": round(rhi,3),
                      "boot_mdiff_lo": round(mlo,3), "boot_mdiff_hi": round(mhi,3)})

pd.DataFrame(out_lines).to_csv(
    os.path.join(revision_out_dir(), "s2_clustered_stats.tsv"), sep="\t", index=False)
print("\n[saved] s2_clustered_stats.tsv")
