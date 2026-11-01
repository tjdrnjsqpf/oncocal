"""R2-A — Symmetric rescue accounting at matched operating points (released seed).

What this replaces
------------------
src/models/fig4_multiseed.py L141-150 counts rescue in one direction only:

    am_missed = positives with AlphaMissense raw < 0.5
    rescue    = fraction of those with OncoCal >= 0.5

It never counts positives that AlphaMissense caught and OncoCal loses. At that operating
point the two models sit at very different places on their own scales, so the comparison
is uncalibrated: AlphaMissense runs at sensitivity 0.79 / specificity 0.90 while OncoCal
runs at 0.75 / 0.97. Counting both directions there gives a NET LOSS. Matching the
operating points reverses that and gives a larger, honest rescue figure.

This script works from `tables/oof_predictions.tsv` — the single representative seed that
is released as S9 Table — so every number it prints can be recomputed by a reader from the
released data alone. (A 10-seed variant was prototyped but each seed costs ~2.5 h on a
shared machine for a refinement that does not change the conclusion; the released seed is
also the one a reviewer can actually check.)

Run:  python scripts/revision_r1/r2a_rescue_repseed.py
Out:  tables/r1/r2a_rescue_repseed.tsv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.paths import find_table  # noqa: E402

OUT_DIR = os.path.join(ROOT, "tables", "r1")
AM_CUT = 0.5                       # the cut the published rescue rule uses
CORE = ["oncogene", "TSG", "both"]
MODELS = ["p_model_tools", "p_model_full"]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    oof = pd.read_csv(find_table("oof_predictions.tsv"), sep="\t")
    rows = []
    for set_nm, sub in [("all_roles", oof),
                        ("core_roles", oof[oof.role_class.isin(CORE)])]:
        pos, neg = sub.y_true == 1, sub.y_true == 0
        am = sub.AlphaMissense.fillna(0.5)
        am_spec = float((am[neg] < AM_CUT).mean())
        am_sens = float((am[pos] >= AM_CUT).mean())
        missed = am[pos] < AM_CUT
        for col in MODELS:
            m = sub[col]
            # threshold at which OncoCal's specificity equals AlphaMissense's at its own cut
            cut_matched = float(np.quantile(m[neg], am_spec))
            for scheme, cut in [("published_0.5", AM_CUT), ("spec_matched", cut_matched)]:
                gained = int((missed & (m[pos] >= cut)).sum())
                lost = int((~missed & (m[pos] < cut)).sum())
                rows.append(dict(
                    role_set=set_nm, model=col, scheme=scheme, model_cut=round(cut, 4),
                    n_am_missed=int(missed.sum()), gained=gained, lost=lost,
                    net=gained - lost,
                    rescue_pct=round(100 * gained / max(int(missed.sum()), 1), 1),
                    model_sens=round(float((m[pos] >= cut).mean()), 4),
                    model_spec=round(float((m[neg] < cut).mean()), 4),
                    am_sens=round(am_sens, 4), am_spec=round(am_spec, 4)))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT_DIR, "r2a_rescue_repseed.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(res.to_string(index=False))

    # JAK2 V617F expressed relative to each model's own threshold: the familiar
    # "0.334 -> 0.57" comparison is between two different scales.
    j = oof[oof.pkey == "JAK2:p.V617F"]
    if len(j):
        j = j.iloc[0]
        core = oof[oof.role_class.isin(CORE)]
        cneg = core.y_true == 0
        am_spec = float((core.AlphaMissense.fillna(0.5)[cneg] < AM_CUT).mean())
        print("\nJAK2 V617F relative to each model's own threshold (core roles, matched specificity):")
        print("  AlphaMissense %.4f / %.3f = %.2fx" % (j.AlphaMissense, AM_CUT, j.AlphaMissense / AM_CUT))
        for col in MODELS:
            cut = float(np.quantile(core[col][cneg], am_spec))
            print("  %-14s %.4f / %.4f = %.2fx" % (col, j[col], cut, j[col] / cut))
        print("  p_model_mech  %.4f (below 0.5: not recovered under the mechanism feature set)"
              % j.p_model_mech)
    print("\n[saved] tables/r1/r2a_rescue_repseed.tsv")


if __name__ == "__main__":
    main()
