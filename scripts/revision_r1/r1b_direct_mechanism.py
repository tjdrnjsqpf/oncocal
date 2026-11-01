"""R1-B — Direct mechanism labels vs gene-role proxy (reviewer #1 comment 1).

Reviewer: "the main conceptual issue is the use of oncogene/TSG status as a proxy for
GOF/LOF mechanism ... only 49 variants have direct mechanism annotation ... statements
comparing 'GOF variants' with 'LOF variants' should distinguish direct mechanism labels
from gene-role proxies throughout the abstract, Results and figures."

Purpose
-------
Quantify, for the subset that HAS a variant-level CIViC mechanism annotation, whether the
under-scoring pattern seen with the gene-role proxy also holds with direct labels. With
~49 variants this is underpowered by design; the point is to report it honestly and to
state how much of the main analysis rests on the proxy.

Outputs
  r1b_mechanism_census.tsv   how many positives have direct GOF / LOF / DN labels, by gene
                             and by gene-role; agreement between direct label and proxy
  r1b_direct_scores.tsv      per-tool median score GOF vs LOF on direct labels only,
                             Mann-Whitney p (uncorrected; exploratory)
  r1b_proxy_agreement.tsv    cross-tab direct mechanism x role_class
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)

AM_LB, AM_LP = 0.34, 0.564  # AlphaMissense official class boundaries

# rb_GOF_lower = 1 - 2U/(n1*n2) with U = mannwhitneyu(GOF, LOF).statistic
#   POSITIVE -> GOF scored LOWER than LOF (the manuscript's claim)


def rank_cols(df):
    return [c for c in df.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]


def short(c):
    return c.replace("dbnsfp_", "").replace("_rankscore", "").replace("_converted", "")


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    pos = df[df.oncogenic_label == "positive"].drop_duplicates("pkey").copy()
    print("unique oncogenic positives = %d" % len(pos))

    mech_col = "civic_mechanism" if "civic_mechanism" in pos.columns else "mechanism"
    pos["direct_mech"] = pos[mech_col]
    has = pos[pos.direct_mech.notna()].copy()

    # Two denominators circulate; report both so the manuscript number is unambiguous.
    # The manuscript's "49 CIViC mechanism annotations ... across 17 genes" counts GOF+LOF
    # among ALL deduplicated variants carrying a mechanism (label-agnostic), excluding DN.
    alld = df[df[mech_col].notna()].drop_duplicates("pkey")
    cnt = alld[mech_col].value_counts().to_dict()
    gof, lof, dn = cnt.get("GOF", 0), cnt.get("LOF", 0), cnt.get("DN", 0)
    print("ALL variants with a CIViC mechanism (pkey-dedup): %d across %d genes"
          % (len(alld), alld.gene.nunique()))
    print("   GOF %d | LOF %d | DN %d   -> GOF+LOF = %d  <-- the manuscript's '49'"
          % (gof, lof, dn, gof + lof))
    print("   NOTE: DN (%d) is omitted from that count, so '49' reads as if it were the "
          "total (%d). State all three classes in the revision." % (dn, len(alld)))
    print("POSITIVES ONLY with a CIViC mechanism: %d across %d genes"
          % (len(has), has.gene.nunique()))
    print(has.direct_mech.value_counts().to_string())

    # ---- census: direct labels by gene and role ----
    census = (has.groupby(["gene", "role_class", "direct_mech"]).size()
              .rename("n").reset_index().sort_values(["gene", "direct_mech"]))
    census.to_csv(os.path.join(OUT, "r1b_mechanism_census.tsv"), sep="\t", index=False)

    # ---- agreement: direct mechanism vs gene-role proxy ----
    xt = pd.crosstab(has.direct_mech, has.role_class, margins=True)
    xt.to_csv(os.path.join(OUT, "r1b_proxy_agreement.tsv"), sep="\t")
    print("\ndirect mechanism x gene-role proxy:")
    print(xt.to_string())
    # concordance for the arms the manuscript equates (GOF~oncogene, LOF~TSG)
    gof_onc = int(((has.direct_mech == "GOF") & (has.role_class == "oncogene")).sum())
    lof_tsg = int(((has.direct_mech == "LOF") & (has.role_class == "TSG")).sum())
    n_arm = int(has.direct_mech.isin(["GOF", "LOF"]).sum())
    print("\nproxy concordance on GOF/LOF-labelled variants: (GOF&oncogene)+(LOF&TSG) = %d/%d = %.1f%%"
          % (gof_onc + lof_tsg, n_arm, 100 * (gof_onc + lof_tsg) / n_arm if n_arm else np.nan))

    # ---- per-tool GOF vs LOF on direct labels only ----
    gof = has[has.direct_mech == "GOF"]
    lof = has[has.direct_mech == "LOF"]
    rows = []
    for c in rank_cols(df) + ["am_pathogenicity"]:
        if c not in has.columns:
            continue
        a, b = gof[c].dropna(), lof[c].dropna()
        if len(a) < 5 or len(b) < 5:
            continue
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"tool": short(c), "n_GOF": len(a), "n_LOF": len(b),
                     "median_GOF": round(a.median(), 4), "median_LOF": round(b.median(), 4),
                     "median_diff_GOF_minus_LOF": round(a.median() - b.median(), 4),
                     "rb_GOF_lower": round(1 - 2 * u / (len(a) * len(b)), 3),
                     "mwu_p_uncorrected": p})
    res = pd.DataFrame(rows).sort_values("median_diff_GOF_minus_LOF")
    res.to_csv(os.path.join(OUT, "r1b_direct_scores.tsv"), sep="\t", index=False)
    print("\nper-tool GOF vs LOF on DIRECT labels (n_GOF=%d, n_LOF=%d) -- most GOF-penalising first:"
          % (len(gof), len(lof)))
    pd.set_option("display.width", 220)
    print(res.head(12).to_string(index=False))

    # ---- AlphaMissense class distribution on direct labels ----
    if "am_pathogenicity" in has.columns:
        has["am_class"] = pd.cut(has.am_pathogenicity, [-np.inf, AM_LB, AM_LP, np.inf],
                                 right=False,
                                 labels=["likely_benign", "ambiguous", "likely_pathogenic"])
        print("\nAlphaMissense official class by direct mechanism:")
        print(pd.crosstab(has.direct_mech, has.am_class, margins=True).to_string())
    print("\n[saved] r1_out/r1b_{mechanism_census,proxy_agreement,direct_scores}.tsv")


if __name__ == "__main__":
    main()
