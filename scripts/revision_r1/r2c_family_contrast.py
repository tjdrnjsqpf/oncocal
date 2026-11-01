"""R2-C — Pre-specified family contrast for the oncogene-vs-TSG gap, plus honest per-tool statistics.

Why this exists
---------------
r2b showed that the panel-average gap has a gene-clustered CI spanning zero, while ten
individual tools had CIs excluding zero. Reporting "ten tools are significant" is not
defensible: an independent statistical review found that none of the ten survives
Benjamini-Hochberg adjustment (smallest q ~ 0.17), and the count itself moves with the
deduplication convention (10 / 6 / 4) and with the CI method (percentile vs BCa).

The fix is to stop running 49 marginal tests and instead ask ONE pre-specified question:

    Is the oncogene-vs-TSG gap larger in predictors that were NOT trained on clinical
    labels (conservation, evolutionary and protein-language models) than in supervised
    meta-predictors?

That is exactly the claim the paper makes in the next section, and the classification is
already fixed and published in S3 Table (type2: supervised n=28, unsupervised n=21), so it
is pre-specified rather than chosen after seeing these results. One contrast, one CI, no
multiplicity problem.

Also emitted, so nothing is hidden:
  - per-tool percentile two-sided p and Benjamini-Hochberg q (expected: no tool survives)
  - P(count of oncogene-worse tools >= 25), i.e. how often the direction is even a majority
  - both deduplication conventions

Run:  python r2c_family_contrast.py
Out:  r1_out/r2c_family_contrast.tsv   family means, contrast, count stability
      r1_out/r2c_per_tool_adjusted.tsv per-tool gap, CI, raw p, BH q, Holm p
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)

B = 2000
SEED = 20261031
MAJORITY = 25          # > half of 49


def find_table(name):
    for sub in ("results/tables", "paper_handoff/tables", "tables"):
        p = os.path.join(ROOT, sub, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(name)


def auc(y, s):
    P = int(y.sum())
    N = len(y) - P
    if P == 0 or N == 0:
        return np.nan
    r = rankdata(s)
    return (r[y == 1].sum() - P * (P + 1) / 2.0) / (P * N)


def bh(p):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(len(p), np.nan)
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(p[idx])]
    m = len(order)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        val = p[i] * m / (rank + 1)
        prev = min(prev, val)
        q[i] = min(prev, 1.0)
    return q


def holm(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    out = np.full(len(p), np.nan)
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(p[idx])]
    m = len(order)
    running = 0.0
    for rank, i in enumerate(order):
        val = p[i] * (m - rank)
        running = max(running, val)
        out[i] = min(running, 1.0)
    return out


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    base = df[df.y.notna()].copy()
    base["y"] = base.y.astype(int)
    tools = [c for c in base.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]
    short = [c.replace("dbnsfp_", "").replace("_rankscore", "").replace("_converted", "")
             for c in tools]

    # Pre-specified families: exactly the published S3 Table classification.
    cls = pd.read_csv(find_table("driver_miss_rate.tsv"), sep="\t")
    type2 = dict(zip(cls.tool, cls.type2))
    fam = np.array([type2.get(s, "unclassified") for s in short])
    unsup = fam == "unsupervised"
    sup = fam == "supervised"
    print("families from S3 Table: unsupervised %d | supervised %d | unclassified %d"
          % (unsup.sum(), sup.sum(), (~unsup & ~sup).sum()))
    if (~unsup & ~sup).any():
        print("  unclassified:", [s for s, f in zip(short, fam) if f == "unclassified"])

    rows_fam, rows_tool = [], []
    for conv, lab in [("as_published_no_dedup", base),
                      ("pkey_dedup", base.drop_duplicates("pkey"))]:
        sub = lab[lab.role_class.isin(["oncogene", "TSG"])].copy()
        is_onc = (sub.role_class == "oncogene").values
        y = sub.y.values
        S = sub[tools].astype(float).values
        genes = sub.gene.values
        gl = np.array(sorted(pd.unique(genes)))
        idx_by_gene = {g: np.flatnonzero(genes == g) for g in gl}
        print("\n##### %s: rows=%d genes=%d" % (conv, len(sub), len(gl)))

        def gaps_for(rows):
            yy, oo = y[rows], is_onc[rows]
            out = np.full(len(tools), np.nan)
            for j in range(len(tools)):
                col = S[rows, j]
                ok = ~np.isnan(col)
                a = auc(yy[ok & oo], col[ok & oo])
                b = auc(yy[ok & ~oo], col[ok & ~oo])
                if not (np.isnan(a) or np.isnan(b)):
                    out[j] = b - a
            return out

        point = gaps_for(np.arange(len(sub)))
        rng = np.random.RandomState(SEED)
        boot = np.full((B, len(tools)), np.nan)
        for b in range(B):
            pick = rng.choice(gl, len(gl), replace=True)
            boot[b] = gaps_for(np.concatenate([idx_by_gene[g] for g in pick]))
            if (b + 1) % 500 == 0:
                print("  bootstrap %d/%d" % (b + 1, B), flush=True)
        np.save(os.path.join(OUT, "r2c_boot_%s.npy" % conv), boot)

        # ---------- family contrast: ONE pre-specified test ----------
        fu = np.nanmean(boot[:, unsup], axis=1)
        fs = np.nanmean(boot[:, sup], axis=1)
        contrast = fu - fs
        pu, ps = np.nanmean(point[unsup]), np.nanmean(point[sup])
        for nm, arr, pt in [("unsupervised_family_mean_gap", fu, pu),
                            ("supervised_family_mean_gap", fs, ps),
                            ("family_contrast_unsup_minus_sup", contrast, pu - ps)]:
            lo, hi = np.nanpercentile(arr, [2.5, 97.5])
            rows_fam.append({"convention": conv, "quantity": nm, "point": round(float(pt), 4),
                             "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
                             "P_gt_0": round(float(np.nanmean(arr > 0)), 4),
                             "excludes_zero": "YES" if lo > 0 else ("NEG" if hi < 0 else "no")})
        counts = np.nansum(boot > 0, axis=1)
        clo, chi = np.percentile(counts, [2.5, 97.5])
        rows_fam.append({"convention": conv, "quantity": "n_tools_oncogene_worse",
                         "point": int(np.nansum(point > 0)), "ci_lo": float(clo), "ci_hi": float(chi),
                         "P_gt_0": round(float(np.mean(counts >= MAJORITY)), 4),
                         "excludes_zero": "P_gt_0 column = P(count >= %d)" % MAJORITY})

        # ---------- per-tool, with multiplicity made explicit ----------
        praw = np.full(len(tools), np.nan)
        for j in range(len(tools)):
            col = boot[:, j][~np.isnan(boot[:, j])]
            if len(col) < 100:
                continue
            praw[j] = 2 * min(np.mean(col > 0), np.mean(col < 0))
            praw[j] = min(praw[j], 1.0)
        q, h = bh(praw), holm(praw)
        for j, nm in enumerate(short):
            col = boot[:, j][~np.isnan(boot[:, j])]
            if not len(col):
                continue
            lo, hi = np.percentile(col, [2.5, 97.5])
            rows_tool.append({"convention": conv, "tool": nm, "family": fam[j],
                              "gap": round(float(point[j]), 4),
                              "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
                              "p_raw_percentile": round(float(praw[j]), 4),
                              "q_BH": round(float(q[j]), 4), "p_Holm": round(float(h[j]), 4),
                              "ci_excludes_zero": "YES" if lo > 0 else ("NEG" if hi < 0 else "no"),
                              "survives_BH_q05": "YES" if q[j] < 0.05 else "no"})
        nsig = int(((np.array([r["convention"] for r in rows_tool]) == conv) &
                    (np.array([r["q_BH"] for r in rows_tool]) < 0.05)).sum())
        print("  per-tool: CI excludes 0 in %d tools; survives BH q<0.05 in %d; smallest q = %.3f"
              % (sum(1 for r in rows_tool if r["convention"] == conv and r["ci_excludes_zero"] != "no"),
                 nsig, np.nanmin(q)))

    fam_df = pd.DataFrame(rows_fam)
    tool_df = pd.DataFrame(rows_tool)
    fam_df.to_csv(os.path.join(OUT, "r2c_family_contrast.tsv"), sep="\t", index=False)
    tool_df.to_csv(os.path.join(OUT, "r2c_per_tool_adjusted.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 240)
    print("\n================ FAMILY CONTRAST ================")
    print(fam_df.to_string(index=False))
    print("\n[saved] r1_out/r2c_family_contrast.tsv, r2c_per_tool_adjusted.tsv, r2c_boot_*.npy")


if __name__ == "__main__":
    main()
