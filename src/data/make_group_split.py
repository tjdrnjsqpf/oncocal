"""Gene/position 단위 group-split 인덱스 생성 (data leakage / circularity 차단).

핵심 함정 (research_plan §6, gap_analysis, Raimondi 2021):
  - dbNSFP ensemble 점수(REVEL 등)는 다른 도구를 feature 로 쓰고, 도구들은 '기능'이 아니라
    '드라이버 유전자'를 인식하는 경향. 같은 유전자 변이가 train/test 에 섞이면 성능 과대평가.
  - 따라서 **gene 단위 group split** 이 기본(test 유전자는 학습에서 미관측).
  - position 단위(gene:pos)는 덜 엄격한 보조 split.

방법: StratifiedGroupKFold (y 로 stratify, group 으로 분리) → fold 인덱스 부여 +
      fold 0 을 held-out test 로 표시.

출력: data/processed/split_<level>.tsv  (pkey, group, y, fold, is_test)
사용: python -m src.data.make_group_split --level gene --folds 5
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
    # ⚠️ 모델링 단위 = unique pkey(단백질 변이). 같은 AA를 만드는 여러 코돈(유전체 변이)이
    #    중복 pkey로 들어오면 같은 변이가 train/test에 갈려 leakage → pkey로 dedup.
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
    print("  fold별 (n, pos율):")
    for fold in sorted(split["fold"].unique()):
        s = split[split["fold"] == fold]
        print("    fold %d: n=%6d  pos=%.3f  genes=%d"
              % (fold, len(s), s["y"].mean(), s["gene"].nunique()))
    # 누수 점검: 같은 group 이 여러 fold 에 있으면 안 됨
    leak = split.groupby("group")["fold"].nunique()
    print("  group이 2개+ fold에 걸친 수(=누수, 0이어야 함): %d" % int((leak > 1).sum()))
    print("[split] wrote -> %s" % out)


if __name__ == "__main__":
    main()
