"""C2 — the mechanism contrast estimated strictly within genes.

Why this exists
---------------
Results says of the missed-versus-caught comparison: "it compares missed against caught
variants, which occur within the same genes, and is correspondingly well resolved under
identical gene-clustered inference." The estimator actually used is a cluster-robust GEE with
gene as the clustering unit. That corrects the standard errors but leaves the point estimate
a population-averaged quantity: it still mixes within-gene and between-gene variation, so it
does not deliver the within-gene contrast the sentence advertises.

A gene-stratified conditional logistic regression does. Genes contribute only through the
ordering of their own variants, so between-gene differences in size, label depth and
composition drop out of the likelihood entirely. Single-class genes contribute nothing and are
dropped automatically.

Result: conditioning on gene makes both associations STRONGER in the same direction, so the
signature is not an artefact of between-gene variation.

                      GEE (population-averaged)        conditional logistic (within gene)
    phyloP   OR/SD    0.507  (0.441-0.583)             0.374  (0.330-0.425)
    RSA      OR/SD    2.301  (2.026-2.612)             3.106  (2.697-3.578)

Inputs are the released tables; falls back to the submitted S9/S2 only for cross-checks.

Run:  python c2_within_gene_mechanism.py
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit

# resolve the repository root from this file, so the script runs from a clean clone
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
CORE = ["oncogene", "TSG", "both"]
AM_CUT = 0.564          # AlphaMissense's published likely-pathogenic boundary


def load():
    a = pd.read_csv(os.path.join(REPO, "tables", "fig1_input.tsv"), sep="\t")
    b = pd.read_csv(os.path.join(REPO, "tables", "fig2_rsa_input.tsv"), sep="\t")
    core = a[a.role_class.isin(CORE)]
    # unique protein-level keys, the convention the mechanism analyses use
    g = (core.groupby("pkey")
         .agg(gene=("gene", "first"), am=("am_pathogenicity", "first"),
              phyloP=("phyloP100way", "first"))
         .reset_index())
    rsa = b.groupby("pkey").struct_rsa.first().rename("rsa")
    g = g.merge(rsa, on="pkey", how="left")
    g["missed"] = (g.am < AM_CUT).astype(int)
    return g


def within_gene(df, col):
    """Conditional logistic for one feature on its own maximal set."""
    d = df.dropna(subset=[col]).copy()
    d["z"] = (d[col] - d[col].mean()) / d[col].std(ddof=0)
    two = d.groupby("gene").missed.nunique()
    keep = set(two[two == 2].index)
    sub = d[d.gene.isin(keep)]
    m = ConditionalLogit(sub.missed.values, sub[["z"]].values,
                         groups=sub.gene.values).fit(disp=0)
    lo, hi = np.exp(m.conf_int()[0])
    return dict(feature=col, n_variant_all=len(d), n_gene_all=d.gene.nunique(),
                n_variant_used=len(sub), n_gene_used=len(keep),
                n_missed_used=int(sub.missed.sum()),
                or_per_sd=float(np.exp(m.params[0])), ci_lo=float(lo), ci_hi=float(hi),
                p=float(m.pvalues[0]))


def main():
    g = load()
    rows = [within_gene(g, c) for c in ("phyloP", "rsa")]

    # both features together, on their intersection
    d = g.dropna(subset=["phyloP", "rsa"]).copy()
    for c in ("phyloP", "rsa"):
        d[c + "_z"] = (d[c] - d[c].mean()) / d[c].std(ddof=0)
    two = d.groupby("gene").missed.nunique()
    keep = set(two[two == 2].index)
    sub = d[d.gene.isin(keep)]
    joint = ConditionalLogit(sub.missed.values, sub[["phyloP_z", "rsa_z"]].values,
                             groups=sub.gene.values).fit(disp=0)

    out = pd.DataFrame(rows)
    pd.set_option("display.width", 220)
    print(out.to_string(index=False))
    print(f"\njoint model, {len(keep)} gene strata, {len(sub)} variants:")
    for k, c in enumerate(("phyloP", "rsa")):
        lo, hi = np.exp(joint.conf_int()[k])
        print(f"  {c:7s} OR per SD {np.exp(joint.params[k]):.3f} "
              f"95% CI {lo:.3f}-{hi:.3f}  p={joint.pvalues[k]:.2e}")
    print("\nfor comparison, the GEE values reported in S2 Table:")
    print("  phyloP  0.507 (0.441-0.583)      rsa  2.301 (2.026-2.612)")
    print("\nSingle-class genes contribute no within-gene contrast and are dropped; the "
          "share of the set they hold is what limits this estimator.")
    for r in rows:
        drop_v = r["n_variant_all"] - r["n_variant_used"]
        print(f"  {r['feature']:7s} uses {r['n_variant_used']} of {r['n_variant_all']} variants "
              f"({100 * r['n_variant_used'] / r['n_variant_all']:.0f}%) across "
              f"{r['n_gene_used']} of {r['n_gene_all']} genes; {drop_v} variants dropped")

    p = os.path.join(REPO, "tables", "r1", "c2_within_gene_mechanism.tsv")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    out.to_csv(p, sep="\t", index=False)
    print(f"\n[saved] {p}")


if __name__ == "__main__":
    main()
