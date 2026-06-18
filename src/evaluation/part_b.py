"""Part B — Discordance & 과소평가 분석 (다도구, dbNSFP rankscore).

산출 (results/tables/, results/figures/):
  B1. 도구별 성능표: oncogenicity AUROC/AUPRC (전체·oncogene·TSG)  -> tool_performance.tsv
  B2. driver 과소평가율: 고recurrence oncogenic positive 를 각 도구가 rankscore<0.5 로 놓치는 비율
      (비지도 구조/진화 도구 vs 지도 앙상블 대비)                  -> driver_miss_rate.tsv
  B3. '왜' 검증: oncogenic positive 중 AlphaMissense 놓침 vs 잡힘의 보존도(GERP/phyloP) 차이
      (Mann-Whitney) + 보존도가 과소평가를 설명하는 정도            -> 콘솔 + why_conservation.tsv
  B5. canonical driver 패널 다도구 표                              -> canonical_driver_panel.tsv

사용: python -m src.evaluation.part_b
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBL = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")

# 도구 분류(거칠지만 유용): 비지도 진화/구조/LM vs 지도 앙상블/supervised
UNSUP = ["AlphaMissense", "ESM1b_converted", "popEVE_converted", "PrimateAI", "SIFT_converted",
         "PROVEAN_converted", "GERPpp_RS", "phyloP100way_vertebrate", "phastCons100way_vertebrate"]
SUP = ["REVEL", "MetaRNN", "MetaSVM", "MetaLR", "VEST4", "M-CAP", "ClinPred", "BayesDel_addAF",
       "CADD_raw", "DEOGEN2", "MVP", "MPC", "LIST-S2", "VARITY_ER", "MutPred2"]


def _rank_cols(df):
    return [c for c in df.columns if c.startswith("dbnsfp_") and c.endswith("_rankscore")]


def _short(c):
    return c.replace("dbnsfp_", "").replace("_rankscore", "").replace("_converted", "")


def b1_tool_performance(lab):
    rows = []
    for c in _rank_cols(lab):
        for grp, d in [("all", lab), ("oncogene", lab[lab.role_class == "oncogene"]),
                       ("TSG", lab[lab.role_class == "TSG"])]:
            s = d[[c, "y"]].dropna()
            if s["y"].nunique() < 2 or len(s) < 50:
                continue
            rows.append({"tool": _short(c), "group": grp, "n": len(s),
                         "pos": int(s["y"].sum()),
                         "AUROC": roc_auc_score(s["y"], s[c]),
                         "AUPRC": average_precision_score(s["y"], s[c])})
    perf = pd.DataFrame(rows)
    wide = perf[perf.group == "all"].set_index("tool")[["AUROC", "AUPRC", "n", "pos"]].sort_values("AUROC", ascending=False)
    return perf, wide


def b2_driver_miss_rate(df, rec_min=50):
    drivers = df[(df.oncogenic_label == "positive") & (df.cosmic_sample_count >= rec_min)]
    rows = []
    for c in _rank_cols(df):
        s = drivers[c].dropna()
        if len(s) < 20:
            continue
        tool = _short(c)
        rows.append({"tool": tool, "n_drivers": len(s),
                     "miss_rate(<0.5)": float((s < 0.5).mean()),
                     "median_rankscore": float(s.median()),
                     "type": "unsup" if tool in UNSUP else ("sup" if tool in SUP else "other")})
    mr = pd.DataFrame(rows).sort_values("miss_rate(<0.5)", ascending=False)
    return mr, drivers


def b3_why_conservation(df):
    pos = df[df.oncogenic_label == "positive"].copy()
    pos["am_missed"] = pos["am_pathogenicity"] < 0.5  # 브리지 raw 점수 기준
    out = {}
    for cons in ["dbnsfp_GERPpp_RS", "dbnsfp_phyloP100way_vertebrate"]:
        miss = pos.loc[pos.am_missed, cons].dropna()
        caught = pos.loc[~pos.am_missed, cons].dropna()
        if len(miss) < 10 or len(caught) < 10:
            continue
        u, p = mannwhitneyu(miss, caught, alternative="two-sided")
        out[cons] = {"missed_median": float(miss.median()), "caught_median": float(caught.median()),
                     "missed_n": len(miss), "caught_n": len(caught), "mannwhitney_p": float(p)}
    return pd.DataFrame(out).T, pos


def b5_canonical_panel(df):
    panel = ["JAK2:p.V617F", "MPL:p.W515L", "MPL:p.W515K", "FGFR3:p.Y373C", "FGFR3:p.S249C",
             "STAT3:p.Y640F", "KIT:p.D816V", "FLT3:p.D835Y", "GNAQ:p.Q209L", "GNA11:p.Q209L",
             "BRAF:p.V600E", "EGFR:p.L858R", "KRAS:p.G12D", "IDH1:p.R132H", "PIK3CA:p.H1047R",
             "NRAS:p.Q61K", "KIT:p.V559D", "ERBB2:p.R678Q"]
    cols = ["pkey", "role_class", "cosmic_sample_count", "am_pathogenicity",
            "dbnsfp_AlphaMissense_rankscore", "dbnsfp_REVEL_rankscore", "dbnsfp_CADD_raw_rankscore",
            "dbnsfp_ESM1b_converted_rankscore", "dbnsfp_popEVE_converted_rankscore",
            "dbnsfp_MetaRNN_rankscore", "dbnsfp_GERPpp_RS", "dbnsfp_phyloP100way_vertebrate"]
    sub = df[df.pkey.isin(panel)][cols].copy()
    sub.columns = [_short(c) if c.startswith("dbnsfp_") else c for c in sub.columns]
    return sub.sort_values("cosmic_sample_count", ascending=False)


def main():
    os.makedirs(TBL, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))
    lab = df[df.y.notna()].copy(); lab["y"] = lab["y"].astype(int)
    print("labeled=%d (pos=%d, neg=%d) | rankscore 도구=%d"
          % (len(lab), int(lab.y.sum()), int((lab.y == 0).sum()), len(_rank_cols(lab))))

    perf, wide = b1_tool_performance(lab)
    perf.to_csv(os.path.join(TBL, "tool_performance.tsv"), sep="\t", index=False)
    print("\n[B1] oncogenicity AUROC 상위/하위 (전체):")
    print(wide.head(8).round(3).to_string()); print("  ...")
    print(wide.tail(5).round(3).to_string())
    og = perf[perf.group == "oncogene"].set_index("tool")["AUROC"]
    ts = perf[perf.group == "TSG"].set_index("tool")["AUROC"]
    comp = pd.concat([og.rename("oncogene"), ts.rename("TSG")], axis=1).dropna()
    comp["TSG-oncogene"] = comp["TSG"] - comp["oncogene"]
    print("\n[B1] oncogene vs TSG AUROC (양수=oncogene를 더 못잡음): 평균 차 %.3f, %d/%d 도구가 oncogene<TSG"
          % (comp["TSG-oncogene"].mean(), int((comp["TSG-oncogene"] > 0).sum()), len(comp)))

    mr, drivers = b2_driver_miss_rate(df)
    mr.to_csv(os.path.join(TBL, "driver_miss_rate.tsv"), sep="\t", index=False)
    print("\n[B2] 고recurrence(>=50) driver %d개 — 도구별 과소평가율(rankscore<0.5):" % len(drivers))
    print("  과소평가 최다 5:"); print(mr.head(5)[["tool", "miss_rate(<0.5)", "type"]].to_string(index=False))
    print("  과소평가 최소 5:"); print(mr.tail(5)[["tool", "miss_rate(<0.5)", "type"]].to_string(index=False))
    for t in ["unsup", "sup"]:
        m = mr[mr.type == t]["miss_rate(<0.5)"]
        if len(m): print("  평균 과소평가율 [%s]: %.3f (n도구=%d)" % (t, m.mean(), len(m)))

    why, pos = b3_why_conservation(df)
    why.to_csv(os.path.join(TBL, "why_conservation.tsv"), sep="\t")
    print("\n[B3] '왜' — AlphaMissense 놓친 vs 잡은 oncogenic positive 의 보존도:")
    print(why.round(3).to_string())

    panel = b5_canonical_panel(df)
    panel.to_csv(os.path.join(TBL, "canonical_driver_panel.tsv"), sep="\t", index=False)
    print("\n[B5] canonical driver 패널 -> results/tables/canonical_driver_panel.tsv (%d 변이)" % len(panel))
    print("[done] 표 4종 -> results/tables/")


if __name__ == "__main__":
    main()
