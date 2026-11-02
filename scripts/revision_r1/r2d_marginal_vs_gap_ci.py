"""R2-D — Are separate x and y error bars on Fig 1A honest?

Fig 1A plots AUROC_TSG (x) against AUROC_oncogene (y) per tool. The revision reports a
gene-clustered CI on the GAP (TSG - oncogene). Drawing that as a vertical bar implies the
x value is exact; drawing a 2D cross from the two marginal CIs implies the gap is as
uncertain as the marginals. Both are wrong if the two AUROCs are correlated across gene
resamples, which they almost certainly are: a resample containing easy genes lifts both.

This measures it. For every tool it stores, per gene-bootstrap replicate, AUROC_TSG,
AUROC_oncogene and their difference, then reports the marginal CI widths, the gap CI
width, and the correlation between the two AUROCs across replicates.

Run:  python r2d_marginal_vs_gap_ci.py
Out:  r1_out/r2d_marginal_vs_gap.tsv
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)
B, SEED = 2000, 20261031


def auc(y, s):
    P = int(y.sum()); N = len(y) - P
    if P == 0 or N == 0:
        return np.nan
    r = rankdata(s)
    return (r[y == 1].sum() - P * (P + 1) / 2.0) / (P * N)


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    lab = df[df.y.notna()].copy(); lab["y"] = lab.y.astype(int)
    tools = [c for c in lab.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]
    short = [c.replace("dbnsfp_", "").replace("_rankscore", "").replace("_converted", "")
             for c in tools]
    sub = lab[lab.role_class.isin(["oncogene", "TSG"])].copy()
    is_onc = (sub.role_class == "oncogene").values
    y = sub.y.values; S = sub[tools].astype(float).values
    genes = sub.gene.values
    gl = np.array(sorted(pd.unique(genes)))
    idx = {g: np.flatnonzero(genes == g) for g in gl}
    print("rows=%d genes=%d tools=%d" % (len(sub), len(gl), len(tools)))

    def both(rows):
        yy, oo = y[rows], is_onc[rows]
        a = np.full(len(tools), np.nan); b = np.full(len(tools), np.nan)
        for j in range(len(tools)):
            col = S[rows, j]; ok = ~np.isnan(col)
            a[j] = auc(yy[ok & oo], col[ok & oo])        # oncogene
            b[j] = auc(yy[ok & ~oo], col[ok & ~oo])      # TSG
        return a, b

    pa, pb = both(np.arange(len(sub)))
    rng = np.random.RandomState(SEED)
    BA = np.full((B, len(tools)), np.nan); BB = np.full((B, len(tools)), np.nan)
    for i in range(B):
        pick = rng.choice(gl, len(gl), replace=True)
        rows = np.concatenate([idx[g] for g in pick])
        BA[i], BB[i] = both(rows)
        if (i + 1) % 500 == 0:
            print("  %d/%d" % (i + 1, B), flush=True)

    rows_out = []
    for j, nm in enumerate(short):
        a, b = BA[:, j], BB[:, j]
        ok = ~np.isnan(a) & ~np.isnan(b)
        if ok.sum() < 100:
            continue
        a, b = a[ok], b[ok]; g = b - a
        alo, ahi = np.percentile(a, [2.5, 97.5])
        blo, bhi = np.percentile(b, [2.5, 97.5])
        glo, ghi = np.percentile(g, [2.5, 97.5])
        rows_out.append({"tool": nm,
                         "auroc_oncogene": round(pa[j], 4), "onc_ci_lo": round(alo, 4),
                         "onc_ci_hi": round(ahi, 4), "onc_width": round(ahi - alo, 4),
                         "auroc_TSG": round(pb[j], 4), "tsg_ci_lo": round(blo, 4),
                         "tsg_ci_hi": round(bhi, 4), "tsg_width": round(bhi - blo, 4),
                         "gap": round(pb[j] - pa[j], 4), "gap_ci_lo": round(glo, 4),
                         "gap_ci_hi": round(ghi, 4), "gap_width": round(ghi - glo, 4),
                         "corr_onc_tsg": round(float(np.corrcoef(a, b)[0, 1]), 3)})
    t = pd.DataFrame(rows_out)
    t["width_ratio_marginal_over_gap"] = (t[["onc_width", "tsg_width"]].mean(axis=1)
                                          / t.gap_width).round(2)
    t.to_csv(os.path.join(OUT, "r2d_marginal_vs_gap.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 260)
    print("\n================ MARGINAL vs GAP ================")
    print("correlation between the two AUROCs across gene resamples: median %.3f (IQR %.3f-%.3f)"
          % (t.corr_onc_tsg.median(), t.corr_onc_tsg.quantile(.25), t.corr_onc_tsg.quantile(.75)))
    print("CI width  oncogene median %.4f | TSG median %.4f | GAP median %.4f"
          % (t.onc_width.median(), t.tsg_width.median(), t.gap_width.median()))
    print("marginal width / gap width: median %.2fx" % t.width_ratio_marginal_over_gap.median())
    print()
    print(t.nlargest(6, "gap")[["tool", "auroc_oncogene", "onc_width", "auroc_TSG",
                                "tsg_width", "gap", "gap_width", "corr_onc_tsg"]].to_string(index=False))
    print("\n[saved] tables/r1/r2d_marginal_vs_gap.tsv")


if __name__ == "__main__":
    main()
