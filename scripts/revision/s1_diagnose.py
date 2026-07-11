"""S1 — Reconciliation diagnostic for AM-missed N.
Directly reproduce, from the dataset, the definitions behind the two manuscript figures (472/2415 vs 1091/3077, total 4168).
Usage: python paper_handoff/revision/s1_diagnose.py
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))

RANK = "dbnsfp_AlphaMissense_rankscore"
RAW = "am_pathogenicity"
PHY = "dbnsfp_phyloP100way_vertebrate"

print("=" * 70)
print("Full dataset: %d rows, unique pkey %d" % (len(df), df.pkey.nunique()))

# oncogenic positive set (two criteria)
pos_all = df[df.oncogenic_label == "positive"]
pos_dedup = pos_all.drop_duplicates("pkey")
print("\noncogenic_label==positive:")
print("  rows(non-dedup)    = %d" % len(pos_all))
print("  unique pkey(dedup) = %d" % len(pos_dedup))
print("  y==1 & label==pos  = %d (rows)" % len(df[(df.oncogenic_label=='positive') & (df.y==1)]))

def summ(name, sub, missed_col, thr=0.5):
    """missed/total by missed_col, and N after handling phyloP/RSA missingness."""
    s = sub[missed_col]
    have = s.notna()
    n_have = have.sum()
    missed = (s < thr) & have
    n_missed = missed.sum()
    rate = n_missed / n_have if n_have else float('nan')
    print("  [%s] %s<%.2f: missed=%d / have-score=%d (%.1f%%) | missing=%d"
          % (name, missed_col, thr, n_missed, n_have, 100*rate, (~have).sum()))
    return n_missed, n_have

print("\n" + "=" * 70)
print("A. RAW score(am_pathogenicity) < 0.5 criterion")
print("-- non-dedup(rows) --")
summ("rows/raw", pos_all, RAW)
print("-- dedup(unique pkey) --")
summ("dedup/raw", pos_dedup, RAW)

print("\nB. RANKSCORE(AlphaMissense_rankscore) < 0.5 criterion")
print("-- non-dedup(rows) --")
summ("rows/rank", pos_all, RANK)
print("-- dedup(unique pkey) --")
summ("dedup/rank", pos_dedup, RANK)

# What if restricted to positives within the training set (y.notna)?
print("\n" + "=" * 70)
print("C. Restricted to training set (y==1, positive)")
lab_pos = df[(df.y == 1)].drop_duplicates("pkey")
print("  training positive unique pkey = %d" % len(lab_pos))
summ("labpos/raw", lab_pos, RAW)
summ("labpos/rank", lab_pos, RANK)

# After handling phyloP missingness (b3 method: non-dedup, raw<0.5, phyloP dropna)
print("\n" + "=" * 70)
print("D. Reproduce b3_why_conservation (non-dedup, raw<0.5, phyloP present)")
pos = pos_all.copy()
pos["am_missed"] = pos[RAW] < 0.5
for cons in [PHY, "dbnsfp_GERPpp_RS"]:
    miss = pos.loc[pos.am_missed, cons].dropna()
    caught = pos.loc[~pos.am_missed, cons].dropna()
    print("  %s: missed n=%d (med %.2f) | caught n=%d (med %.2f)"
          % (cons.replace('dbnsfp_',''), len(miss), miss.median(), len(caught), caught.median()))

print("\nE. dedup version phyloP")
posd = pos_dedup.copy(); posd["am_missed"] = posd[RAW] < 0.5
for cons in [PHY, "dbnsfp_GERPpp_RS"]:
    miss = posd.loc[posd.am_missed, cons].dropna()
    caught = posd.loc[~posd.am_missed, cons].dropna()
    print("  %s: missed n=%d (med %.2f) | caught n=%d (med %.2f)"
          % (cons.replace('dbnsfp_',''), len(miss), miss.median(), len(caught), caught.median()))
