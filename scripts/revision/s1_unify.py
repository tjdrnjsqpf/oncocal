"""S1 — Recompute mechanism values after unifying to a single definition.
Definition: oncogenic positive **unique variant (dedup pkey)**, missed = AM RAW(am_pathogenicity) < 0.5.
Recompute phyloP/RSA/GERP all from the same set to eliminate mixing of dedup/non-dedup.
Usage: python paper_handoff/revision/s1_unify.py
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from src.features.structure_features import add_structure_features

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
RAW = "am_pathogenicity"; RANK = "dbnsfp_AlphaMissense_rankscore"

pos = df[df.oncogenic_label == "positive"].drop_duplicates("pkey").copy()
pos = add_structure_features(pos)
pos["missed"] = pos[RAW] < 0.5
n = len(pos); nm = int(pos.missed.sum()); nc = n - nm
print("Unified definition: unique oncogenic positive = %d | missed(raw<0.5)=%d (%.1f%%) | caught=%d"
      % (n, nm, 100*nm/n, nc))
print("(note) same set rankscore<0.5 = %d/%d (%.1f%%)"
      % (int((pos[RANK] < 0.5).sum()), pos[RANK].notna().sum(),
         100*(pos[RANK] < 0.5).mean()))

rows = []
for name, col in [("phyloP100way", "dbnsfp_phyloP100way_vertebrate"),
                  ("GERP++_RS", "dbnsfp_GERPpp_RS"),
                  ("RSA", "struct_rsa")]:
    miss = pos.loc[pos.missed, col].dropna()
    caught = pos.loc[~pos.missed, col].dropna()
    u, p = mannwhitneyu(miss, caught, alternative="two-sided")
    # rank-biserial r
    r = 1 - 2*u/(len(miss)*len(caught))
    rows.append({"axis": name, "missed_n": len(miss), "missed_median": round(miss.median(), 3),
                 "caught_n": len(caught), "caught_median": round(caught.median(), 3),
                 "median_diff": round(caught.median()-miss.median(), 3),
                 "rank_biserial_r": round(r, 3), "mannwhitney_p": p})
res = pd.DataFrame(rows)
print("\nMechanism recomputed under unified definition (naive Mann-Whitney, clustering-corrected in S2):")
print(res.to_string(index=False))
res.to_csv(os.path.join(ROOT, "paper_handoff", "revision", "s1_unified_mechanism.tsv"),
           sep="\t", index=False)
print("\n[saved] s1_unified_mechanism.tsv")
