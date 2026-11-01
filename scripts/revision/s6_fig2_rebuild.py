"""S6 — Regenerate Figure 2 (old build 2,415 -> unified 3,922).
Unified definition (S1): unique oncogenic positive (dedup pkey), missed = AM RAW<0.5.
Outputs:
  [A] fig2_master_3922.tsv — shared input for Fig2A scatter (AM raw vs COSMIC recurrence) + Fig2B (phyloP vs RSA)
  [B] recompute corner statistics — the "79% missed vs 14% elsewhere" values for the unified set
Usage: python scripts/revision/s6_fig2_rebuild.py
"""
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np, pandas as pd
from src.features.structure_features import add_structure_features

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
from src.paths import revision_out_dir
OUT = revision_out_dir()
df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))

pos = df[df.oncogenic_label == "positive"].drop_duplicates("pkey").copy()
pos = add_structure_features(pos)
pos["AM_status"] = np.where(pos.am_pathogenicity < 0.5, "AM missed", "AM caught")

master = pos[["pkey", "gene", "role_class", "cosmic_sample_count", "am_pathogenicity",
              "AM_status", "dbnsfp_phyloP100way_vertebrate", "dbnsfp_GERPpp_RS",
              "struct_rsa"]].rename(columns={
    "dbnsfp_phyloP100way_vertebrate": "phyloP100way", "dbnsfp_GERPpp_RS": "GERPpp_RS"})
master.to_csv(os.path.join(OUT, "fig2_master_3922.tsv"), sep="\t", index=False)

n_total = len(master)
n_missed = int((master.AM_status == "AM missed").sum())
print("=" * 64)
print("Unified Fig2 input: N=%d positives | AM-missed=%d (%.1f%%) | caught=%d"
      % (n_total, n_missed, 100*n_missed/n_total, n_total-n_missed))
print("  (old build legend 2,415 -> replaced with %d; missed rate 19.5%%->%.1f%%)" % (n_total, 100*n_missed/n_total))
print("  Fig2A(AM raw vs recurrence): am_pathogenicity/cosmic_sample_count missing %d/%d"
      % (master.am_pathogenicity.isna().sum(), master.cosmic_sample_count.isna().sum()))
print("  Fig2B(phyloP vs RSA): both present = %d" % ((master.phyloP100way.notna() & master.struct_rsa.notna()).sum()))

# ---------- recompute corner statistics ----------
print("\n" + "=" * 64)
print("[corner statistics] low-conservation & exposed corner = (phyloP100way<4 & RSA>0.4)")
d = master.dropna(subset=["phyloP100way", "struct_rsa"]).copy()
d["missed"] = d.AM_status == "AM missed"
d["corner"] = (d.phyloP100way < 4) & (d.struct_rsa > 0.4)
cin = d[d.corner]; cout = d[~d.corner]
p_in = 100*cin.missed.mean(); p_out = 100*cout.missed.mean()
print("  Fig2B target(phyloP/RSA both present) N=%d" % len(d))
print("  inside corner   n=%4d | missed %.1f%%   <- manuscript '79%% missed'" % (len(cin), p_in))
print("  outside corner  n=%4d | missed %.1f%%   <- manuscript '14%% elsewhere'" % (len(cout), p_out))
# reverse direction (note): fraction of missed variants located in the corner
print("  (note) fraction of all missed located in corner = %.1f%%"
      % (100*cin.missed.sum()/d.missed.sum()))

# threshold sensitivity (manuscript definition may differ, so list several cuts)
print("\n[sensitivity] by corner definition (missed_in%% / missed_out%% / n_in):")
for pc, rc in [(4,0.4),(3,0.4),(4,0.5),(2,0.5),(4,0.35)]:
    cc = (d.phyloP100way<pc)&(d.struct_rsa>rc)
    print("  phyloP<%.1f & RSA>%.2f : %.1f%% / %.1f%% / n_in=%d"
          % (pc, rc, 100*d[cc].missed.mean(), 100*d[~cc].missed.mean(), cc.sum()))

pd.DataFrame([{"corner_def":"phyloP<4 & RSA>0.4","N_fig2B":len(d),
               "n_corner":len(cin),"missed_in_pct":round(p_in,1),
               "n_elsewhere":len(cout),"missed_elsewhere_pct":round(p_out,1)}]
             ).to_csv(os.path.join(OUT,"s6_corner_stats.tsv"),sep="\t",index=False)
print("\n[saved] fig2_master_3922.tsv, s6_corner_stats.tsv")
