"""R1-A — Label-source decomposition (reviewer #1 comment 2).

Reviewer: "COSMIC recurrence >= 10 is not equivalent to experimentally demonstrated
oncogenicity. Could the author show the main GOF under-scoring result separately for
CIViC-oncogenic, cancerhotspots and recurrence-only variants? A compact table with
effect sizes for each label source would make the benchmark much easier to judge."

Design
------
The label source only affects POSITIVES. Negatives (ClinVar benign / gnomAD common) are
held fixed in every stratum, so each stratum answers: "if we only trusted this one source
of oncogenicity, would the oncogene-vs-TSG AUROC gap still appear?"

Strata (positives):
  all            every positive (current manuscript definition)
  civic           CIViC oncogenic (CC0, curated)
  hotspot         cancerhotspots significant residue
  recurrence      COSMIC sample count >= 10
  civic_only      CIViC and NOT hotspot and NOT recurrence>=10
  hotspot_only    hotspot and NOT CIViC and NOT recurrence>=10
  recurrence_only recurrence>=10 and NOT CIViC and NOT hotspot

Per stratum we recompute, for all 49 dbNSFP rank-score tools:
  oncogene AUROC, TSG AUROC, gap (TSG - oncogene), n tools with oncogene < TSG,
  exact sign test p, mean/median gap.
Only role_class oncogene and TSG enter this comparison (dual-role, fusion-only and
unknown genes cannot be assigned to either arm) -- this matches the published Fig 1A.

Run:  python r1a_label_source.py
Out:  r1_out/r1a_label_source.tsv        per-stratum, per-tool AUROC
      r1_out/r1a_label_source_summary.tsv per-stratum effect sizes (paper table)
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "tables", "r1")
os.makedirs(OUT, exist_ok=True)

RECURRENCE_MIN = 10  # must match build_dataset.py --recurrence-min used for the manuscript


def rank_cols(df):
    return [c for c in df.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]


def short(c):
    return c.replace("dbnsfp_", "").replace("_rankscore", "").replace("_converted", "")


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    tools = rank_cols(df)
    print("dataset=%d rows | rankscore tools=%d" % (len(df), len(tools)))

    pos = df.oncogenic_label == "positive"
    neg = df.oncogenic_label == "negative"
    civic = df.get("civic_oncogenic", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    hot = df.get("is_hotspot", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    rec = df.get("cosmic_sample_count", pd.Series(0, index=df.index)).fillna(0) >= RECURRENCE_MIN

    strata = {
        "all": pos,
        "civic": pos & civic,
        "hotspot": pos & hot,
        "recurrence": pos & rec,
        "civic_only": pos & civic & ~hot & ~rec,
        "hotspot_only": pos & hot & ~civic & ~rec,
        "recurrence_only": pos & rec & ~civic & ~hot,
    }

    rows, summary = [], []
    for nm, mask in strata.items():
        sub = df[mask | neg].copy()
        sub["y"] = sub["oncogenic_label"].map({"positive": 1, "negative": 0}).astype(int)
        n_pos = int(mask.sum())
        for c in tools:
            rec_row = {"stratum": nm, "tool": short(c), "n_pos": n_pos}
            for grp, d in [("oncogene", sub[sub.role_class == "oncogene"]),
                           ("TSG", sub[sub.role_class == "TSG"])]:
                s = d[[c, "y"]].dropna()
                ok = s.y.nunique() == 2 and len(s) >= 50
                rec_row["auroc_%s" % grp] = roc_auc_score(s.y, s[c]) if ok else np.nan
                rec_row["n_%s" % grp] = len(s)
                rec_row["pos_%s" % grp] = int(s.y.sum())
            rows.append(rec_row)

        st = pd.DataFrame([r for r in rows if r["stratum"] == nm]).dropna(
            subset=["auroc_oncogene", "auroc_TSG"])
        gap = st.auroc_TSG - st.auroc_oncogene
        n_worse, n_tools = int((gap > 0).sum()), len(gap)
        p = binomtest(n_worse, n_tools, 0.5, alternative="two-sided").pvalue if n_tools else np.nan
        summary.append({
            "stratum": nm, "n_positive": n_pos, "n_tools": n_tools,
            "n_oncogene_worse": n_worse,
            "frac_oncogene_worse": round(n_worse / n_tools, 3) if n_tools else np.nan,
            "mean_gap_TSG_minus_oncogene": round(gap.mean(), 4) if n_tools else np.nan,
            "median_gap": round(gap.median(), 4) if n_tools else np.nan,
            "gap_q25": round(gap.quantile(.25), 4) if n_tools else np.nan,
            "gap_q75": round(gap.quantile(.75), 4) if n_tools else np.nan,
            "sign_test_p": p,
        })
        print("  %-16s n_pos=%5d | oncogene<TSG in %2d/%2d tools | mean gap %+.4f | sign p=%.3g"
              % (nm, n_pos, n_worse, n_tools, gap.mean() if n_tools else np.nan, p))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "r1a_label_source.tsv"), sep="\t", index=False)
    s = pd.DataFrame(summary)
    s.to_csv(os.path.join(OUT, "r1a_label_source_summary.tsv"), sep="\t", index=False)
    print("\n[saved] r1_out/r1a_label_source{,_summary}.tsv")
    print("\n--- summary (this is the compact table the reviewer asked for) ---")
    pd.set_option("display.width", 220)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
