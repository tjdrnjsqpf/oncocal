"""Build the unified dataset — join scores, labels, and recurrence to the panel missense universe by pkey.

Architecture (2026-06-16, after COSMIC arrival):
  variant universe = AlphaMissense panel missense (data/interim/alphamissense_panel.tsv)
    → 'all possible missense' within the panel + flagship scores. Expand multi-tool score columns when dbNSFP arrives.
  Left-join onto it:
    - CIViC      : civic_oncogenic, civic_mechanism
    - hotspots   : is_hotspot, hotspot_qvalue
    - ClinVar    : clinvar_benign
    - COSMIC rec : cosmic_sample_count (MutantCensus, for label definition)
    - panel      : role_class, mech_proxy (gene-level)

Labeling rule (reflecting civic-label limitations + JAK2 V617F rescue):
  positive: is_hotspot | civic_oncogenic | (cosmic_sample_count >= --recurrence-min)
            └ JAK2 V617F is missing from hotspots but is rescued by COSMIC recurrence (42,930).
  negative: clinvar_benign AND no positive signal
  conflict: positive and negative simultaneously → excluded from training
  unknown : no signal → excluded from training (never treat as negative — key pitfall)

Mechanism: prefer civic_mechanism (per-variant); otherwise panel mech_proxy (per-gene, weak prior).

⚠️ circularity (§6): cosmic_sample_count is used for **label definition** → excluded from the default features.
   When used as a feature, only the ablation arm uses GenomeScreens (unbiased) recurrence.
⚠️ leakage: gene/position group-split is generated separately in make_group_split.py.

Output: data/processed/dataset.parquet (full universe) + labels_summary.tsv (labeled subset)
Usage: python -m src.data.build_dataset --recurrence-min 10
"""
from __future__ import annotations
import argparse
import os

import pandas as pd

INTERIM = "data/interim"
PROCESSED = "data/processed"


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _load(name: str, root: str):
    p = os.path.join(root, INTERIM, name)
    if not os.path.exists(p):
        print("[build] (not found, skip): %s" % p)
        return None
    return pd.read_csv(p, sep="\t")


def build(root: str, recurrence_min: int, gnomad_common_af: float = 1e-3) -> pd.DataFrame:
    am = _load("alphamissense_panel.tsv", root)
    if am is None:
        raise SystemExit("alphamissense_panel.tsv not found — run alphamissense_join first.")
    base = am.copy()

    civic = _load("civic_labels.tsv", root)
    hot = _load("cancerhotspots_labels.tsv", root)
    clin = _load("clinvar_labels.tsv", root)
    rec = _load("cosmic_recurrence__mutantcensus.tsv", root)
    panel = _load("cancer_gene_panel.tsv", root)
    dbnsfp = _load("dbnsfp_panel.tsv", root)  # auto-join on arrival (skip if absent)

    if civic is not None:
        base = base.merge(civic[["pkey", "civic_oncogenic", "civic_mechanism", "civic_onc_support"]],
                          on="pkey", how="left")
    if hot is not None:
        base = base.merge(hot[["pkey", "is_hotspot", "hotspot_qvalue", "hotspot_alt_count"]],
                          on="pkey", how="left")
    if clin is not None:
        base = base.merge(clin[["pkey", "clinvar_benign", "clinvar_significance"]],
                          on="pkey", how="left")
    if rec is not None:
        # ⚠️ recurrence is per GRCh38 genome coordinate (avoids isoform numbering issues) → join with AM (GRCh38) by coordinate.
        gk = ["grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt"]
        for c in ["grch38_chr", "grch38_ref", "grch38_alt"]:
            base[c] = base[c].astype(str)
            rec[c] = rec[c].astype(str)
        base["grch38_pos"] = base["grch38_pos"].astype(int)
        rec["grch38_pos"] = rec["grch38_pos"].astype(int)
        base = base.merge(rec[gk + ["cosmic_sample_count"]], on=gk, how="left")
    if panel is not None:
        base = base.merge(panel[["gene", "role_class", "mech_proxy", "tier"]], on="gene", how="left")
    if dbnsfp is not None:
        # Multi-tool scores/conservation/gnomAD AF — GRCh38 coordinate join (adds 40+ tool columns, no structural change)
        gk = ["grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt"]
        for c in ["grch38_chr", "grch38_ref", "grch38_alt"]:
            base[c] = base[c].astype(str); dbnsfp[c] = dbnsfp[c].astype(str)
        base["grch38_pos"] = base["grch38_pos"].astype(int)
        dbnsfp["grch38_pos"] = dbnsfp["grch38_pos"].astype(int)
        drop = [c for c in ["aaref", "aaalt", "dbnsfp_gene"] if c in dbnsfp.columns]
        base = base.merge(dbnsfp.drop(columns=drop), on=gk, how="left")
        print("[build] joined %d dbNSFP score columns" % sum(c.startswith("dbnsfp_") for c in base.columns))

    for c, d in [("civic_oncogenic", False), ("is_hotspot", False), ("clinvar_benign", False)]:
        base[c] = base.get(c, d)
        base[c] = base[c].fillna(False).astype(bool)
    base["cosmic_sample_count"] = base.get("cosmic_sample_count", 0)
    base["cosmic_sample_count"] = base["cosmic_sample_count"].fillna(0).astype(int)

    # ---- Labels ----
    base["recurrent_driver"] = base["cosmic_sample_count"] >= recurrence_min
    base["pos_signal"] = base["is_hotspot"] | base["civic_oncogenic"] | base["recurrent_driver"]
    # gnomAD common variant = benign negative supplement (only when dbNSFP present)
    af_col = "dbnsfp_gnomAD4.1_joint_AF"
    if af_col in base.columns:
        base["gnomad_common"] = pd.to_numeric(base[af_col], errors="coerce").fillna(0) >= gnomad_common_af
    else:
        base["gnomad_common"] = False
    base["neg_signal"] = base["clinvar_benign"] | base["gnomad_common"]
    base["label_conflict"] = base["pos_signal"] & base["neg_signal"]

    def _lab(r):
        if r["label_conflict"]:
            return "conflict"
        if r["pos_signal"]:
            return "positive"
        if r["neg_signal"]:
            return "negative"
        return "unknown"
    base["oncogenic_label"] = base.apply(_lab, axis=1)
    base["y"] = base["oncogenic_label"].map({"positive": 1, "negative": 0})

    # ---- Mechanism ----
    base["mechanism"] = base.get("civic_mechanism")
    base["mechanism"] = base["mechanism"].fillna(base.get("mech_proxy"))
    return base


def summarize(df: pd.DataFrame, recurrence_min: int) -> None:
    print("[build] universe (panel missense w/ AM): %d  (genes=%d)" % (len(df), df["gene"].nunique()))
    print(df["oncogenic_label"].value_counts().to_string())
    tr = df[df["y"].notna()]
    print("\ntrainable pos=%d / neg=%d  (recurrence_min=%d)"
          % (int((tr["y"] == 1).sum()), int((tr["y"] == 0).sum()), recurrence_min))
    print("  positive source: hotspot=%d, CIViC=%d, recurrent(>=%d)=%d"
          % (int(df["is_hotspot"].sum()), int(df["civic_oncogenic"].sum()),
             recurrence_min, int(df["recurrent_driver"].sum())))
    # Key driver check: under-scoring preview
    canon = ["JAK2:p.V617F", "BRAF:p.V600E", "EGFR:p.L858R", "KRAS:p.G12D", "IDH1:p.R132H"]
    sub = df[df["pkey"].isin(canon)][["pkey", "am_pathogenicity", "am_class",
                                      "cosmic_sample_count", "is_hotspot", "oncogenic_label"]]
    if len(sub):
        print("\n[under-scoring preview] canonical driver — AM score vs recurrence:")
        print(sub.to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description="Build unified panel dataset (AM universe + labels + recurrence).")
    root = _repo_root()
    ap.add_argument("--recurrence-min", type=int, default=10,
                    help="if COSMIC distinct-sample count >= this value, candidate recurrent driver (positive)")
    ap.add_argument("--gnomad-common-af", type=float, default=1e-3,
                    help="if gnomAD joint AF >= this value, common=benign negative (when dbNSFP present)")
    ap.add_argument("--out-parquet", default=os.path.join(root, PROCESSED, "dataset.parquet"))
    ap.add_argument("--out-labeled", default=os.path.join(root, PROCESSED, "labels_summary.tsv"))
    args = ap.parse_args()

    df = build(root, args.recurrence_min, args.gnomad_common_af)
    os.makedirs(os.path.join(root, PROCESSED), exist_ok=True)
    df.to_parquet(args.out_parquet, index=False)
    df[df["y"].notna()].to_csv(args.out_labeled, sep="\t", index=False)
    summarize(df, args.recurrence_min)
    print("\n[build] wrote -> %s (full), %s (labeled subset)" % (args.out_parquet, args.out_labeled))


if __name__ == "__main__":
    main()
