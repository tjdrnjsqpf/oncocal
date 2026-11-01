"""R2-B — Two-stage (gene-clustered) uncertainty for the 42/49 oncogene-vs-TSG gap.

Problem this fixes
------------------
The manuscript reports "42 of 49 tools ... mean gap 0.030 (sign test p = 3.6e-7)".
The sign test's observation unit is the 49 TOOLS, not variants or genes. The tools score
the same positives and the dbNSFP meta-predictors use each other as features, so they are
heavily correlated: one realisation of the variant sample is amplified 49-fold. That p-value
therefore carries no information about sampling uncertainty in the variants/genes, which is
the uncertainty a reader cares about.

What this script does
---------------------
Keeps the gene as the independent unit (same convention as s2_clustered.py / r1f / r1g) and
resamples GENES with replacement. Within each replicate it recomputes ALL 49 per-tool gaps,
which yields three things the manuscript currently lacks:

  (1) per-tool gap with a gene-clustered 95% CI  -> how many tools' gaps actually exclude 0
  (2) the MEAN gap across 49 tools with a gene-clustered CI (the two-stage estimate)
  (3) the COUNT "oncogene < TSG" with a gene-clustered CI (how stable is "42 of 49")

This converts an over-stated significance claim into an honest, defensible one. Note the
direction-consistency statement ("42 of 49 point in the same direction") survives regardless;
what changes is the uncertainty attached to it.

Run:  python r2b_tool_gap_twostage.py
Out:  r1_out/r2b_tool_gap_ci.tsv       per-tool point + CI
      r1_out/r2b_twostage_summary.tsv  mean-gap and count CIs
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


def auc(y, s):
    """Rank-based AUROC with tie correction. y is 0/1 int array, s is float array."""
    n = len(y)
    if n == 0:
        return np.nan
    P = int(y.sum())
    N = n - P
    if P == 0 or N == 0:
        return np.nan
    r = rankdata(s)
    return (r[y == 1].sum() - P * (P + 1) / 2.0) / (P * N)


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    base = df[df.y.notna()].copy()
    base["y"] = base.y.astype(int)
    # Fig 1A (src/evaluation/part_b.py b1_tool_performance) does NOT deduplicate by pkey;
    # the mechanism analysis does. Report both so the manuscript can be internally consistent.
    CONVENTIONS = [("as_published_no_dedup", base),
                   ("pkey_dedup", base.drop_duplicates("pkey"))]
    tools = [c for c in base.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]
    short = [c.replace("dbnsfp_", "").replace("_rankscore", "").replace("_converted", "")
             for c in tools]
    print("tools=%d" % len(tools))
    all_tab, all_summ = [], []
    for conv_name, lab in CONVENTIONS:
      print("\n########## convention = %s (n=%d, pos=%d)" % (conv_name, len(lab), int(lab.y.sum())))
      # Only oncogene and TSG rows enter this comparison (matches part_b.py b1_tool_performance)
      sub = lab[lab.role_class.isin(["oncogene", "TSG"])].copy()
      is_onc = (sub.role_class == "oncogene").values
      y = sub.y.values
      S = sub[tools].astype(float).values           # n x 49
      genes = sub.gene.values
      gene_list = np.array(sorted(pd.unique(genes)))
      idx_by_gene = {g: np.flatnonzero(genes == g) for g in gene_list}
      print("oncogene+TSG rows=%d (oncogene %d / TSG %d) genes=%d"
          % (len(sub), int(is_onc.sum()), int((~is_onc).sum()), len(gene_list)))

      def gaps_for(rows):
        """Per-tool gap = AUROC(TSG) - AUROC(oncogene) on the given row index."""
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

      allrows = np.arange(len(sub))
      point = gaps_for(allrows)
      n_worse = int(np.nansum(point > 0))
      print("\npoint estimate: oncogene<TSG in %d/%d tools | mean gap %+.4f"
          % (n_worse, int(np.sum(~np.isnan(point))), np.nanmean(point)))

      rng = np.random.RandomState(SEED)
      boot = np.full((B, len(tools)), np.nan)
      mean_gap = np.full(B, np.nan)
      counts = np.full(B, np.nan)
      for b in range(B):
        pick = rng.choice(gene_list, len(gene_list), replace=True)
        rows = np.concatenate([idx_by_gene[g] for g in pick])
        g = gaps_for(rows)
        boot[b] = g
        mean_gap[b] = np.nanmean(g)
        counts[b] = np.nansum(g > 0)
        if (b + 1) % 200 == 0:
            print("  bootstrap %d/%d" % (b + 1, B), flush=True)

      rows_out = []
      for j, nm in enumerate(short):
        col = boot[:, j]
        col = col[~np.isnan(col)]
        if not len(col):
            continue
        lo, hi = np.percentile(col, [2.5, 97.5])
        rows_out.append({"tool": nm, "gap_TSG_minus_oncogene": round(float(point[j]), 4),
                         "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
                         "excludes_zero": "YES" if lo > 0 else ("NEG" if hi < 0 else "no"),
                         "P_gap_gt_0": round(float((col > 0).mean()), 4)})
      tab = pd.DataFrame(rows_out).sort_values("gap_TSG_minus_oncogene", ascending=False)
      tab.insert(0, "convention", conv_name)
      all_tab.append(tab)

      mlo, mhi = np.percentile(mean_gap[~np.isnan(mean_gap)], [2.5, 97.5])
      clo, chi = np.percentile(counts[~np.isnan(counts)], [2.5, 97.5])
      n_excl = int((tab.excludes_zero == "YES").sum())
      n_neg = int((tab.excludes_zero == "NEG").sum())
      summ = pd.DataFrame([
        {"quantity": "mean_gap_across_tools", "point": round(float(np.nanmean(point)), 4),
         "ci_lo": round(float(mlo), 4), "ci_hi": round(float(mhi), 4),
         "excludes_zero": "YES" if mlo > 0 else "no",
         "P_gt_0": round(float((mean_gap[~np.isnan(mean_gap)] > 0).mean()), 4)},
        {"quantity": "n_tools_oncogene_worse", "point": n_worse,
         "ci_lo": round(float(clo), 1), "ci_hi": round(float(chi), 1),
         "excludes_zero": "", "P_gt_0": ""},
        {"quantity": "n_tools_with_own_CI_excluding_0_positive", "point": n_excl,
         "ci_lo": "", "ci_hi": "", "excludes_zero": "", "P_gt_0": ""},
        {"quantity": "n_tools_with_own_CI_excluding_0_negative", "point": n_neg,
         "ci_lo": "", "ci_hi": "", "excludes_zero": "", "P_gt_0": ""},
      ])
      summ.insert(0, "convention", conv_name)
      all_summ.append(summ)

      pd.set_option("display.width", 220)
      print("\n================ TWO-STAGE SUMMARY ================")
      print(summ.to_string(index=False))
      print("\n개별 도구 CI가 0을 제외하는 경우: 양의 방향 %d개 / 음의 방향 %d개 / 나머지 %d개는 0 포함"
          % (n_excl, n_neg, len(tab) - n_excl - n_neg))
      print("\n--- gap 상위 8 ---")
      print(tab.head(8).to_string(index=False))
      print("\n--- gap 하위 8 (oncogene에서 오히려 나은 도구) ---")
      print(tab.tail(8).to_string(index=False))
    pd.concat(all_tab).to_csv(os.path.join(OUT, "r2b_tool_gap_ci.tsv"), sep="\t", index=False)
    pd.concat(all_summ).to_csv(os.path.join(OUT, "r2b_twostage_summary.tsv"), sep="\t", index=False)
    print("\n[saved] r1_out/r2b_tool_gap_ci.tsv, r2b_twostage_summary.tsv (양 convention 모두 포함)")


if __name__ == "__main__":
    main()
