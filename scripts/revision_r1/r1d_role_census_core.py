"""R1-D — Gene-role census and core-role primary analysis.

Addresses reviewer #1 minor comment 2 ("clarify how genes annotated as both oncogene and
TSG were handled in mechanism-stratified analyses, and whether fusion-only CGC genes
contribute") and reviewer #2 comment 4 (AlphaMissense threshold is 0.564, not 0.5).

What the local audit found
--------------------------
Of 3,922 unique oncogenic positives: oncogene 1,389, TSG 1,194, both 904, fusion_only 137,
unknown 298. Fig 1A (the 42/49 result) already used only role_class oncogene vs TSG, so it
is unaffected. But Fig 1B (per-tool miss rate) and Fig 2 (phyloP/RSA mechanism analysis)
pooled ALL positives, including fusion_only and unknown -- which are 77% and 74%
AlphaMissense-missed and dominated by long repetitive genes (MUC6 142, HLA-A 31,
PDE4DIP 27). Those inflate the headline miss rate.

Decision (2026-10-31): core roles = oncogene + TSG + both become the PRIMARY set; the full
3,922 set is reported as a sensitivity analysis. AlphaMissense's official class boundaries
(0.34 / 0.564) replace the ad-hoc 0.5 cut.

Outputs
  r1d_role_census.tsv        positives/negatives by role_class x AlphaMissense class
  r1d_suspect_genes.tsv      fusion_only/unknown genes ranked by missed count (artefact audit)
  r1d_core_missed.tsv        miss rate under 3 definitions x 2 role sets
  r1d_core_mechanism.tsv     phyloP/RSA medians, rank-biserial, gene-clustered GEE
  r1d_core_corner.tsv        low-conservation/exposed corner enrichment
  r1d_core_tool_missrate.tsv Fig 1B rebuild: per-tool miss rate, core roles, official cut
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.features.structure_features import add_structure_features  # noqa: E402

OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)


def find_table(name):
    """Locate a pipeline table across repo layouts (server: results/tables or
    paper_handoff/tables; public repo: tables/). Raises with the paths tried."""
    tried = []
    for sub in ("results/tables", "paper_handoff/tables", "tables"):
        p = os.path.join(ROOT, sub, name)
        tried.append(p)
        if os.path.exists(p):
            return p
    raise FileNotFoundError("%s not found; tried:\n  %s" % (name, "\n  ".join(tried)))

AM_LB, AM_LP = 0.34, 0.564
CORE = ["oncogene", "TSG", "both"]
PHY = "dbnsfp_phyloP100way_vertebrate"
RSA = "struct_rsa"
RECURRENCE_FIG1B = 50  # Fig 1B used high-recurrence drivers (>= 50 samples)

# Sign convention, matching scripts/revision/s2_clustered.py:
#   rb_missed_lower = 1 - 2U/(n1*n2) with U = mannwhitneyu(missed, caught).statistic
#   POSITIVE  -> missed variants have the LOWER value (the claim for phyloP)
#   NEGATIVE  -> missed variants have the HIGHER value (the claim for RSA)

DEFS = [("official_notLP", AM_LP, "AlphaMissense below likely-pathogenic (official 0.564)"),
        ("official_LB", AM_LB, "AlphaMissense likely-benign (official 0.34)"),
        ("legacy_0.5", 0.50, "original manuscript cut (0.5)")]


def rank_cols(df):
    return [c for c in df.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]


def short(c):
    return c.replace("dbnsfp_", "").replace("_rankscore", "").replace("_converted", "")


def gee_or(frame, col):
    """Cluster-robust logistic regression of missed on z(col), cluster = gene."""
    d = frame[["gene", "missed", col]].dropna().rename(columns={col: "x"})
    if d.missed.nunique() < 2 or d.gene.nunique() < 10:
        return dict(n=len(d), n_gene=d.gene.nunique(), or_per_sd=np.nan, p=np.nan)
    d["z"] = (d.x - d.x.mean()) / d.x.std()
    g = smf.gee("missed ~ z", groups="gene", data=d, family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Exchangeable()).fit()
    ci = np.exp(g.conf_int().loc["z"])
    return dict(n=len(d), n_gene=d.gene.nunique(), or_per_sd=float(np.exp(g.params["z"])),
                ci_lo=float(ci.iloc[0]), ci_hi=float(ci.iloc[1]), p=float(g.pvalues["z"]))


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    pos = df[df.oncogenic_label == "positive"].drop_duplicates("pkey").copy()
    pos = add_structure_features(pos)
    pos["am_class"] = pd.cut(pos.am_pathogenicity, [-np.inf, AM_LB, AM_LP, np.inf], right=False,
                             labels=["likely_benign", "ambiguous", "likely_pathogenic"])
    print("unique oncogenic positives = %d" % len(pos))

    # ---- (1) census ----
    census = pd.crosstab(pos.role_class, pos.am_class, margins=True)
    census.to_csv(os.path.join(OUT, "r1d_role_census.tsv"), sep="\t")
    print("\n(1) positives by gene role x AlphaMissense official class:")
    print(census.to_string())
    neg = df[df.oncogenic_label == "negative"].drop_duplicates("pkey")
    print("\n    negatives by role_class:")
    print(neg.role_class.value_counts(dropna=False).to_string())

    # ---- (2) suspect genes in the excluded roles ----
    excl = pos[~pos.role_class.isin(CORE)].copy()
    excl["missed_official"] = excl.am_pathogenicity < AM_LP
    susp = (excl.groupby(["gene", "role_class"])
            .agg(n_positive=("pkey", "size"), n_missed=("missed_official", "sum"),
                 median_cosmic=("cosmic_sample_count", "median"))
            .reset_index().sort_values("n_missed", ascending=False))
    susp["frac_missed"] = (susp.n_missed / susp.n_positive).round(3)
    susp.to_csv(os.path.join(OUT, "r1d_suspect_genes.tsv"), sep="\t", index=False)
    print("\n(2) fusion_only/unknown positives = %d across %d genes; top contributors:"
          % (len(excl), excl.gene.nunique()))
    print(susp.head(12).to_string(index=False))

    # ---- (3) miss rate: 3 definitions x 2 role sets ----
    rows = []
    for set_nm, sub in [("core_roles", pos[pos.role_class.isin(CORE)]), ("all_roles", pos)]:
        for def_nm, thr, desc in DEFS:
            s = sub.am_pathogenicity.dropna()
            rows.append({"role_set": set_nm, "missed_def": def_nm, "threshold": thr,
                         "n_scored": len(s), "n_missed": int((s < thr).sum()),
                         "pct_missed": round(100 * (s < thr).mean(), 2), "description": desc})
    miss = pd.DataFrame(rows)
    miss.to_csv(os.path.join(OUT, "r1d_core_missed.tsv"), sep="\t", index=False)
    print("\n(3) AlphaMissense miss rate by definition and role set:")
    pd.set_option("display.width", 250)
    print(miss.to_string(index=False))

    # ---- (4) mechanism analysis (phyloP / RSA) on both role sets ----
    rows = []
    for set_nm, sub0 in [("core_roles", pos[pos.role_class.isin(CORE)]), ("all_roles", pos)]:
        for def_nm, thr, _ in DEFS:
            sub = sub0.copy()
            sub["missed"] = (sub.am_pathogenicity < thr).astype(int)
            for ax, col in [("phyloP", PHY), ("RSA", RSA)]:
                d = sub[["gene", "missed", col]].dropna()
                mi, ca = d.loc[d.missed == 1, col], d.loc[d.missed == 0, col]
                if len(mi) < 10 or len(ca) < 10:
                    continue
                u, p = mannwhitneyu(mi, ca, alternative="two-sided")
                g = gee_or(sub, col)
                rows.append({"role_set": set_nm, "missed_def": def_nm, "axis": ax,
                             "n": len(d), "n_missed": len(mi),
                             "median_missed": round(mi.median(), 3),
                             "median_caught": round(ca.median(), 3),
                             "rb_missed_lower": round(1 - 2 * u / (len(mi) * len(ca)), 3),
                             "mwu_p": p, "gee_or_per_sd": round(g["or_per_sd"], 3),
                             "gee_ci_lo": round(g.get("ci_lo", np.nan), 3),
                             "gee_ci_hi": round(g.get("ci_hi", np.nan), 3),
                             "gee_p": g["p"], "n_gene": g["n_gene"]})
    mech = pd.DataFrame(rows)
    mech.to_csv(os.path.join(OUT, "r1d_core_mechanism.tsv"), sep="\t", index=False)
    print("\n(4) mechanism analysis (gene-clustered GEE):")
    print(mech.to_string(index=False))

    # ---- (5) corner enrichment (phyloP < 4 & RSA > 0.4) ----
    rows = []
    for set_nm, sub0 in [("core_roles", pos[pos.role_class.isin(CORE)]), ("all_roles", pos)]:
        for def_nm, thr, _ in DEFS:
            sub = sub0.dropna(subset=[PHY, RSA]).copy()
            sub["missed"] = sub.am_pathogenicity < thr
            corner = (sub[PHY] < 4) & (sub[RSA] > 0.4)
            rows.append({"role_set": set_nm, "missed_def": def_nm,
                         "n_corner": int(corner.sum()),
                         "pct_missed_in_corner": round(100 * sub.loc[corner, "missed"].mean(), 1),
                         "n_elsewhere": int((~corner).sum()),
                         "pct_missed_elsewhere": round(100 * sub.loc[~corner, "missed"].mean(), 1)})
    corner_df = pd.DataFrame(rows)
    corner_df.to_csv(os.path.join(OUT, "r1d_core_corner.tsv"), sep="\t", index=False)
    print("\n(5) low-conservation / exposed corner (phyloP < 4 & RSA > 0.4):")
    print(corner_df.to_string(index=False))

    # ---- (6) Fig 1B rebuild: per-tool miss rate on core roles ----
    full = df[df.oncogenic_label == "positive"].copy()
    drivers = full[(full.cosmic_sample_count >= RECURRENCE_FIG1B)
                   & full.role_class.isin(CORE)]
    rows = []
    for c in rank_cols(df):
        s = drivers[c].dropna()
        if len(s) < 20:
            continue
        rows.append({"tool": short(c), "n_drivers": len(s),
                     "miss_rate_rank_lt_0.5": round(float((s < 0.5).mean()), 4),
                     "median_rankscore": round(float(s.median()), 4)})
    tool_mr = pd.DataFrame(rows).sort_values("miss_rate_rank_lt_0.5", ascending=False)
    tool_mr.to_csv(os.path.join(OUT, "r1d_core_tool_missrate.tsv"), sep="\t", index=False)
    print("\n(6) Fig 1B rebuild on core roles: %d high-recurrence drivers (>=%d samples)"
          % (len(drivers), RECURRENCE_FIG1B))
    print("    NOTE: merge with the published type/type2 columns of tables/driver_miss_rate.tsv")
    print("    to recompute the unsupervised-vs-supervised means (17.8% vs 9.1%).")
    print(tool_mr.head(8).to_string(index=False))
    print("\n[saved] 6 tables -> r1_out/r1d_*.tsv")


if __name__ == "__main__":
    main()
