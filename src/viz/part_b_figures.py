"""Three core Part B figures (results/figures/, PNG+SVG+PDF).

  fig1_conservation_why  — conservation (phyloP100) of AM-missed vs caught oncogenic drivers = the 'why'
  fig2_tool_miss_rate    — per-tool under-scoring rate of high-recurrence drivers (unsupervised vs supervised)
  fig3_oncogene_vs_tsg   — per-tool oncogene vs TSG AUROC

For the paper, figure text is in English (removes Korean font dependency). Each figure saved as PNG/SVG/PDF.
Usage: python -m src.viz.part_b_figures
"""
from __future__ import annotations
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TBL = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.titlesize": 11, "axes.labelsize": 11,
    "svg.fonttype": "none", "pdf.fonttype": 42,  # editable text
})


def _save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(FIG, "%s.%s" % (name, ext)),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet"))

    # --- fig1: conservation = why ---
    pos = df[df.oncogenic_label == "positive"].copy()
    pos["AM_missed"] = pos["am_pathogenicity"] < 0.5
    miss = pos.loc[pos.AM_missed, "dbnsfp_phyloP100way_vertebrate"].dropna()
    caught = pos.loc[~pos.AM_missed, "dbnsfp_phyloP100way_vertebrate"].dropna()
    import numpy as np
    rng = np.random.default_rng(42)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    # Show individual data points (jitter strip) behind the box
    for xpos, data, c in [(1, caught, "#4C9F70"), (2, miss, "#C0504D")]:
        jx = xpos + rng.normal(0, 0.07, size=len(data))
        ax.scatter(jx, data, s=5, color=c, alpha=0.18, edgecolor="none", zorder=1)
    bp = ax.boxplot([caught, miss], positions=[1, 2],
                    tick_labels=["AM caught\n(n=%d)" % len(caught), "AM missed\n(n=%d)" % len(miss)],
                    patch_artist=True, showfliers=False, widths=0.55, zorder=3)
    for p, c in zip(bp["boxes"], ["#4C9F70", "#C0504D"]):
        p.set_facecolor(c); p.set_alpha(.55); p.set_edgecolor("black")
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(1.8)
    for w in bp["whiskers"] + bp["caps"]:
        w.set_color("black")
    ax.set_ylabel("phyloP100way conservation")
    ax.set_title("Underestimated drivers sit at\nlow-conservation residues (Mann-Whitney p<1e-300)")
    ax.axhline(0, color="grey", lw=.5, ls="--")
    _save(fig, "fig1_conservation_why")

    # --- fig2: per-tool driver underestimation rate ---
    mr = pd.read_csv(os.path.join(TBL, "driver_miss_rate.tsv"), sep="\t").sort_values("miss_rate(<0.5)")
    color = {"unsup": "#C0504D", "sup": "#4C9F70", "other": "#B0B0B0"}
    fig, ax = plt.subplots(figsize=(6.2, 9))
    ax.barh(mr["tool"], mr["miss_rate(<0.5)"], color=[color[t] for t in mr["type"]])
    ax.set_xlabel("Underestimation rate of high-recurrence drivers (rankscore < 0.5)")
    ax.set_title("Per-tool underestimation of cancer drivers")
    ax.legend(handles=[Patch(color="#C0504D", label="unsupervised (structure/evolution)"),
                       Patch(color="#4C9F70", label="supervised ensemble"),
                       Patch(color="#B0B0B0", label="other")], fontsize=9, loc="lower right")
    ax.margins(y=0.01)
    _save(fig, "fig2_tool_miss_rate")

    # --- fig3: oncogene vs TSG AUROC ---
    perf = pd.read_csv(os.path.join(TBL, "tool_performance.tsv"), sep="\t")
    og = perf[perf.group == "oncogene"].set_index("tool")["AUROC"]
    ts = perf[perf.group == "TSG"].set_index("tool")["AUROC"]
    comp = pd.concat([og.rename("og"), ts.rename("ts")], axis=1).dropna()
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(comp["ts"], comp["og"], s=20, color="#34568B", alpha=.8, edgecolor="white", linewidth=.3)
    lo, hi = comp.min().min() - .01, comp.max().max() + .01
    ax.plot([lo, hi], [lo, hi], "k--", lw=.8, label="equal")
    ax.set_xlabel("TSG (LOF) AUROC")
    ax.set_ylabel("Oncogene (GOF) AUROC")
    ax.set_title("Most tools score oncogenes (GOF) worse than TSGs\n(below diagonal: %d/%d tools)"
                 % (int((comp.og < comp.ts).sum()), len(comp)))
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=9)
    _save(fig, "fig3_oncogene_vs_tsg")

    print("[figures] saved PNG+SVG+PDF -> results/figures/")
    for f in sorted(os.listdir(FIG)):
        print("   ", f)


if __name__ == "__main__":
    main()
