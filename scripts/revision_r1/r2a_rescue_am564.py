"""R2-A at AlphaMissense's published boundary (0.564) — the figures Results section 3 now reports.

Why this exists
---------------
`r2a_rescue_repseed.py` hardcodes `AM_CUT = 0.5`, the cut the original submission used. The
revision adopts AlphaMissense's published boundaries everywhere else (Reviewer 2, point 4), so
the rescue paragraph was the one place where the old cut survived — and it mattered: the
Abstract defines the missed set as the 811 positives below 0.564, while the rescue rate was
computed against the 746 below 0.5. A reader computing 315/811 would get 38.8%, not the
quoted 42.2%.

Recomputing at 0.564 leaves the conclusion intact and the denominator consistent:

                   missed   matched cut   gained   lost   net   rescue%   sens vs AM @ spec
    AM_CUT 0.5        746        0.1956      315     88   227     42.2%   0.8512 / 0.7861 @ 0.8980
    AM_CUT 0.564      811        0.2288      336     86   250     41.4%   0.8391 / 0.7674 @ 0.9122

The identity `net = (sens_model - sens_AM) * n_pos` holds exactly at both cuts (227 and 250 on
3,487 core-role positives), which is why the net is reported as the sensitivity difference at
fixed specificity rather than as a count that happens to come out positive.

Method is identical to `r2a_rescue_repseed.py` — same released seed, same specificity matching,
same columns — with the cut taken from argv so both can be produced from one code path.

Run:  python r2a_rescue_am564.py [cut ...]        # default: 0.5 0.564
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CORE = ["oncogene", "TSG", "both"]
MODELS = ["p_model_tools", "p_model_mech", "p_model_full"]
MAIN = "p_model_full"          # the feature set the Results quote


def load_oof():
    """The released seed: tables/oof_predictions.tsv, published as S9 Table."""
    for p in (os.path.abspath(os.path.join(HERE, "..", "..", "tables",
                                           "oof_predictions.tsv")),
              os.path.join(HERE, "oof_predictions.tsv")):
        if os.path.exists(p):
            return pd.read_csv(p, sep="\t")
    # fall back to the submitted S9 Table, so this runs from the submission alone
    import openpyxl
    xl = os.path.join(HERE, "..", "..", "..", "4_Supporting_Information", "S9_Table.xlsx")
    rows = list(openpyxl.load_workbook(xl, read_only=True, data_only=True)
                ["S9_oof_predictions"].iter_rows(values_only=True))
    df = pd.DataFrame([r for r in rows[3:] if r[0]], columns=list(rows[2]))
    df["y_true"] = df.y_true.astype(int)
    for c in MODELS + ["AlphaMissense", "VARITY_R"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def rescue(sub, am_cut, col):
    pos, neg = sub.y_true == 1, sub.y_true == 0
    am = sub.AlphaMissense.fillna(am_cut)
    am_spec = float((am[neg] < am_cut).mean())
    am_sens = float((am[pos] >= am_cut).mean())
    missed = am[pos] < am_cut
    m = sub[col]
    cut = float(np.quantile(m[neg], am_spec))
    gained = int((missed & (m[pos] >= cut)).sum())
    lost = int((~missed & (m[pos] < cut)).sum())
    sens = float((m[pos] >= cut).mean())
    return dict(am_cut=am_cut, model=col, model_cut=round(cut, 4),
                n_missed=int(missed.sum()), gained=gained, lost=lost, net=gained - lost,
                rescue_pct=round(100 * gained / max(int(missed.sum()), 1), 1),
                model_sens=round(sens, 4), am_sens=round(am_sens, 4),
                am_spec=round(am_spec, 4), n_pos=int(pos.sum()),
                sens_diff_x_npos=round((sens - am_sens) * int(pos.sum()), 1))


def gene_bootstrap(core, am_cut, col, B=2000, seed=20261102):
    """C9 — a gene-clustered interval for the rescue quantities.

    The counts were the one place a headline number carried neither an interval nor gene
    clustering. Genes are resampled with replacement and the whole operating-point
    construction is redone inside each replicate: AlphaMissense's specificity at `am_cut`,
    the matched model threshold, then the sensitivity difference. Counts scale with the
    resample, so the quantity bootstrapped is the sensitivity difference at matched
    specificity (which equals net / n_pos exactly) and the rescue rate.
    """
    genes = core.gene.unique()
    where = {g: np.where(core.gene.values == g)[0] for g in genes}
    rng = np.random.default_rng(seed)
    keep = {"d_sens": [], "rate": []}
    for _ in range(B):
        pick = rng.choice(genes, size=len(genes), replace=True)
        rr = core.iloc[np.concatenate([where[g] for g in pick])]
        if int((rr.y_true == 1).sum()) < 10 or int((rr.y_true == 0).sum()) < 10:
            continue
        r = rescue(rr, am_cut, col)
        keep["d_sens"].append(r["model_sens"] - r["am_sens"])
        keep["rate"].append(r["rescue_pct"] / 100)
    ci = {k: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
          for k, v in keep.items()}
    return ci, len(keep["d_sens"])


def main(cuts):
    oof = load_oof()
    core = oof[oof.role_class.isin(CORE)].copy()
    core["gene"] = core.pkey.str.split(":").str[0]
    res = pd.DataFrame([rescue(core, c, m) for c in cuts for m in MODELS])
    pd.set_option("display.width", 250)
    print(res.to_string(index=False))

    for c in cuts:
        ci, n = gene_bootstrap(core, c, MAIN)
        pt = rescue(core, c, MAIN)
        print(f"\ngene-clustered bootstrap at AM_CUT {c} ({MAIN}, {n} replicates over "
              f"{core.gene.nunique()} genes):")
        d = pt["model_sens"] - pt["am_sens"]
        print(f"  sensitivity difference at matched specificity  {d:+.4f}  "
              f"95% CI {ci['d_sens'][0]:+.4f} to {ci['d_sens'][1]:+.4f}"
              f"   {'excludes zero' if ci['d_sens'][0] > 0 else 'includes zero'}")
        print(f"  rescue rate among AlphaMissense-missed        {pt['rescue_pct']:5.1f}%  "
              f"95% CI {100*ci['rate'][0]:.1f}% to {100*ci['rate'][1]:.1f}%")
        print(f"  net as a count on {pt['n_pos']} positives          {d*pt['n_pos']:.0f}  "
              f"95% CI {ci['d_sens'][0]*pt['n_pos']:.0f} to {ci['d_sens'][1]*pt['n_pos']:.0f}")

    j = oof[oof.pkey == "JAK2:p.V617F"]
    if len(j):
        j = j.iloc[0]
        neg = core.y_true == 0
        for c in cuts:
            spec = float((core.AlphaMissense.fillna(c)[neg] < c).mean())
            print(f"\nJAK2 V617F at AM_CUT {c} (matched specificity {spec:.4f}):")
            print(f"  AlphaMissense  {j.AlphaMissense:.3f} / {c} = {j.AlphaMissense / c:.2f}x")
            for col in MODELS:
                cut = float(np.quantile(core[col][neg], spec))
                print(f"  {col:15s}{j[col]:.3f} / {cut:.4f} = {j[col] / cut:.2f}x"
                      f"  {'recovered' if j[col] >= cut else 'not recovered'}")

    out = os.path.join(HERE, "..", "..", "tables", "r1", "r2a_rescue_am564.tsv")
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    print(f"\n[saved] {out}")


if __name__ == "__main__":
    main([float(x) for x in sys.argv[1:]] or [0.5, 0.564])
