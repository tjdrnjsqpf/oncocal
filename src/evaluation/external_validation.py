"""External validation — evaluate the model & tools on independent DMS (ProteinGym) functional data (filling the biggest reviewer gap).

Method: train the full model (tools+mech+struct) on our labels → predict on ProteinGym DMS variants →
  compare correlation (Spearman) with the experimental fitness score (DMS_score) against AlphaMissense and the best tool.
  - DMS is independent of our training labels (a different experimental axis).
  - 6 assays: oncogene (SRC activity, YAP1) + TSG-LOF (PTEN×2, BRCA1, MSH2).
  - ⚠️ For TSG-LOF assays, oncogenicity ∝ -fitness → compare by |Spearman| (report the sign). SRC=activity assay (different GOF direction).
  - The model may have seen some of these variants/genes during training (not gene-group) → the 'relative' improvement over tools is the key point.

Usage: python -m src.evaluation.external_validation
"""
from __future__ import annotations
import glob
import os

import numpy as np
import pandas as pd
import lightgbm as lgb
from scipy.stats import spearmanr

from ..features.mechanism_features import add_mechanism_features, MECH_COLS
from ..features.structure_features import add_structure_features, STRUCT_COLS
from ..data.variant_keys import parse_protein_change, protein_key

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBL = os.path.join(ROOT, "results", "tables")
ADIR = os.path.join(ROOT, "data", "external", "proteingym", "assays")
PARAMS = dict(objective="binary", learning_rate=0.03, num_leaves=31, min_child_samples=40,
              subsample=0.8, subsample_freq=1, colsample_bytree=0.7, reg_lambda=1.0,
              n_estimators=400, verbose=-1)


AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def assay_gene(fn):
    return os.path.basename(fn).split("_")[0].upper()


def canonical_seq(uniprot):
    """Canonical 1-letter sequence from the AF pdb (position=resnum)."""
    p = glob.glob(os.path.join(ROOT, "data", "external", "alphafold", "AF-%s-F1-model_v*.pdb" % uniprot))
    if not p:
        return None
    res = {}
    for line in open(p[0]):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res[int(line[22:26])] = AA3.get(line[17:20].strip(), "X")
    if not res:
        return None
    return "".join(res.get(i, "X") for i in range(1, max(res) + 1))


def assay_offset(target_seq, canon):
    """Where the assay construct (target_seq) starts within canonical → offset to add to the local pos.
    Full-length (0), or substring search using a fragment as seed. Returns None on failure."""
    if not target_seq or not canon:
        return 0
    if target_seq[:40] in canon:
        return canon.index(target_seq[:40])
    if target_seq[10:50] in canon:
        return canon.index(target_seq[10:50]) - 10
    return None


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet")).drop_duplicates("pkey")
    df = add_structure_features(add_mechanism_features(df))
    tool_feats = [c for c in df.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]
    feats = tool_feats + MECH_COLS + STRUCT_COLS

    # Train the full model (labeled portion)
    lab = df[df.y.notna()]
    m = lgb.LGBMClassifier(**PARAMS).fit(lab[feats].astype(float).values, lab.y.astype(int).values)
    df["model_pred"] = m.predict_proba(df[feats].astype(float).values)[:, 1]
    feat_index = df.set_index("pkey")
    AM = "dbnsfp_AlphaMissense_rankscore"
    REVEL = "dbnsfp_REVEL_rankscore"

    # Reference (target_seq) + gene->uniprot (for construct numbering offset correction)
    ref = pd.read_csv(os.path.join(ROOT, "data", "external", "proteingym", "DMS_substitutions.csv"))
    id2seq = dict(zip(ref.DMS_id, ref.target_seq))
    gu = pd.read_csv(os.path.join(ROOT, "data", "interim", "gene_uniprot.tsv"), sep="\t")
    g2u = dict(zip(gu.gene.str.upper(), gu.uniprot.fillna("")))

    rows = []
    for fn in sorted(glob.glob(os.path.join(ADIR, "*.csv"))):
        gene = assay_gene(fn)
        dms_id = os.path.basename(fn).replace(".csv", "")
        canon = canonical_seq(g2u.get(gene, ""))
        off = assay_offset(id2seq.get(dms_id), canon)
        dms = pd.read_csv(fn)
        # Single variants only + construct->canonical position correction + ref AA consistency check
        pk, keep = [], []
        for mut in dms["mutant"].astype(str):
            pc = parse_protein_change(mut)
            ok = False
            if pc and pc[0] != pc[2] and off is not None and canon:
                cpos = pc[1] + off
                if 1 <= cpos <= len(canon) and canon[cpos - 1] == pc[0]:  # ref AA matches canonical
                    pk.append(protein_key(gene, pc[0], cpos, pc[2])); ok = True
            if not ok:
                pk.append(None)
            keep.append(ok)
        dms = dms.assign(pkey=pk)[keep].dropna(subset=["pkey"])
        merged = dms.merge(feat_index[["model_pred", AM, REVEL]], left_on="pkey", right_index=True, how="inner")
        if len(merged) < 30:
            rows.append({"assay": os.path.basename(fn).replace(".csv", ""), "gene": gene,
                         "n_dms": len(dms), "n_covered": len(merged), "note": "insufficient coverage"})
            continue
        s_model = spearmanr(merged.model_pred, merged.DMS_score, nan_policy="omit").correlation
        s_am = spearmanr(merged[AM], merged.DMS_score, nan_policy="omit").correlation
        s_rev = spearmanr(merged[REVEL], merged.DMS_score, nan_policy="omit").correlation
        rows.append({"assay": os.path.basename(fn).replace(".csv", ""), "gene": gene,
                     "n_dms": len(dms), "n_covered": len(merged),
                     "spearman_MODEL": round(s_model, 3), "spearman_AlphaMissense": round(s_am, 3),
                     "spearman_REVEL": round(s_rev, 3)})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(TBL, "external_validation_dms.tsv"), sep="\t", index=False)
    print(res.to_string(index=False))

    ok = res.dropna(subset=["spearman_MODEL"]) if "spearman_MODEL" in res else res
    if "spearman_MODEL" in res.columns and len(ok):
        print("\n|Spearman| mean: MODEL=%.3f, AlphaMissense=%.3f, REVEL=%.3f"
              % (ok.spearman_MODEL.abs().mean(), ok.spearman_AlphaMissense.abs().mean(),
                 ok.spearman_REVEL.abs().mean()))
        print("Assays where MODEL has higher |Spearman|: %d/%d"
              % (int((ok.spearman_MODEL.abs() > ok.spearman_AlphaMissense.abs()).sum()), len(ok)))
    print("\n[done] -> results/tables/external_validation_dms.tsv")


if __name__ == "__main__":
    main()
