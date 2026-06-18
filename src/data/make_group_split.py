"""Build a gene/position-level group-split index (block data leakage / circularity).

Key pitfall (research_plan §6, gap_analysis, Raimondi 2021):
  - dbNSFP ensemble scores (REVEL etc.) use other tools as features, and tools tend to recognize
    the 'driver gene' rather than 'function'. If variants of the same gene mix into train/test, performance is overestimated.
  - Therefore **gene-level group split** is the default (test genes are unseen during training).
  - Position-level (gene:pos) is a less strict auxiliary split.

Method: StratifiedGroupKFold (stratify by y, split by group) → assign fold indices +
      mark fold 0 as the held-out test.

Output: data/processed/split_<level>.tsv  (pkey, group, y, fold, is_test)
Usage: python -m src.data.make_group_split --level gene --folds 5
"""
from __future__ import annotations
import argparse
import os

import sys

import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

try:
    from .variant_keys import parse_protein_change
except ImportError:
    sys.path.insert(0, os.path.dirname(__file__))
    from variant_keys import parse_protein_change

PROCESSED = "data/processed"


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _pos(change: str):
    pc = parse_protein_change(str(change))
    return pc[1] if pc else -1


def make_split(df: pd.DataFrame, level: str, folds: int, seed: int) -> pd.DataFrame:
    df = df[df["y"].notna()].copy()
    # ⚠️ Modeling unit = unique pkey (protein variant). If multiple codons (genomic variants) producing
    #    the same AA enter as duplicate pkeys, the same variant splits across train/test → leakage → dedup by pkey.
    df = df.drop_duplicates("pkey").reset_index(drop=True)
    df["y"] = df["y"].astype(int)
    if "pos" not in df.columns:
        df["pos"] = df["protein_change"].map(_pos)
    if level == "gene":
        df["group"] = df["gene"]
    elif level == "position":
        df["group"] = df["gene"].astype(str) + ":" + df["pos"].astype(str)
    else:
        raise ValueError("level must be gene|position")

    sgkf = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed)
    df["fold"] = -1
    for fold, (_, te) in enumerate(sgkf.split(df, df["y"], df["group"])):
        df.iloc[te, df.columns.get_loc("fold")] = fold
    df["is_test"] = df["fold"] == 0   # fold 0 = held-out test
    return df[["pkey", "gene", "pos", "group", "y", "fold", "is_test"]]


def main():
    ap = argparse.ArgumentParser(description="Gene/position group-split index (anti-leakage).")
    root = _repo_root()
    ap.add_argument("--labeled", default=os.path.join(root, PROCESSED, "labels_summary.tsv"))
    ap.add_argument("--level", choices=["gene", "position"], default="gene")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.labeled, sep="\t")
    split = make_split(df, args.level, args.folds, args.seed)
    out = args.out or os.path.join(root, PROCESSED, "split_%s.tsv" % args.level)
    split.to_csv(out, sep="\t", index=False)

    print("[split] level=%s folds=%d  labeled=%d  groups=%d"
          % (args.level, args.folds, len(split), split["group"].nunique()))
    print("  per fold (n, pos rate):")
    for fold in sorted(split["fold"].unique()):
        s = split[split["fold"] == fold]
        print("    fold %d: n=%6d  pos=%.3f  genes=%d"
              % (fold, len(s), s["y"].mean(), s["gene"].nunique()))
    # Leakage check: the same group must not appear in multiple folds
    leak = split.groupby("group")["fold"].nunique()
    print("  number of groups spanning 2+ folds (=leakage, should be 0): %d" % int((leak > 1).sum()))
    print("[split] wrote -> %s" % out)


if __name__ == "__main__":
    main()
