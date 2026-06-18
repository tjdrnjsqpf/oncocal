"""cancerhotspots.org v2 -> residue-hotspot-based auxiliary positive labels.

Design rationale (phase1_execution.md C, civic-label limitations):
  - CIViC 'Oncogenic' evidence is scarce (~46 missense positives) + canonical drivers unflagged.
    cancerhotspots' statistically significant recurrent hotspot residues are the main pillar of positives.
  - SNV-hotspots sheet: 1 row = (residue, specific alt AA). Columns:
      Reference_Amino_Acid='Q:422' (ref residue:total residue count), Variant_Amino_Acid='R:204' (alt:count for that alt),
      Mutation_Count=total variants at the residue, qvalue/log10_pvalue=hotspot significance.
  - ⚠️ Amino_Acid_Position is string dtype → needs int conversion.
  - ⚠️ Pitfall: hotspot recurrence can enter both feature/label, risking circularity.
    Here we build only 'residue-level hotspot labels' and do not put the recurrence count directly into the label
    (counts kept for analysis/auxiliary use only). Check via ablation in build_dataset.

Output: data/interim/cancerhotspots_labels.tsv (per pkey, alt-specific missense)
Usage: python -m src.data.cancerhotspots_labels
"""
from __future__ import annotations
import argparse
import os
import sys

import pandas as pd

try:
    from .variant_keys import AA1, protein_key
except ImportError:
    sys.path.insert(0, os.path.dirname(__file__))
    from variant_keys import AA1, protein_key


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _split_aa_count(s):
    """'Q:422' -> ('Q', 422). Malformed input returns (None, None)."""
    if not isinstance(s, str) or ":" not in s:
        return None, None
    aa, _, cnt = s.partition(":")
    aa = aa.strip()
    try:
        return aa, int(cnt)
    except ValueError:
        return aa, None


def build_hotspot_labels(xls_path: str, qvalue_max: float = 0.05) -> pd.DataFrame:
    df = pd.read_excel(xls_path, sheet_name="SNV-hotspots")
    rows = []
    for _, r in df.iterrows():
        gene = str(r["Hugo_Symbol"]).strip()
        try:
            pos = int(str(r["Amino_Acid_Position"]).strip())
        except (ValueError, TypeError):
            continue
        ref_aa, residue_count = _split_aa_count(r["Reference_Amino_Acid"])
        alt_aa, alt_count = _split_aa_count(r["Variant_Amino_Acid"])
        # Single missense only: ref/alt standard 20 AA, alt != ref, nonsense(*) excluded
        if not ref_aa or not alt_aa:
            continue
        if ref_aa not in AA1 or alt_aa not in AA1 or alt_aa == ref_aa:
            continue
        qv = pd.to_numeric(r.get("qvalue"), errors="coerce")
        rows.append({
            "pkey": protein_key(gene, ref_aa, pos, alt_aa),
            "gene": gene, "ref_aa": ref_aa, "pos": pos, "alt_aa": alt_aa,
            "protein_change": "%s%d%s" % (ref_aa, pos, alt_aa),
            "hotspot_residue_count": residue_count,   # total variants at the residue (auxiliary, do not use directly as label)
            "hotspot_alt_count": alt_count,            # count of this alt variant
            "hotspot_qvalue": qv,
            "hotspot_log10p": pd.to_numeric(r.get("log10_pvalue"), errors="coerce"),
            "is_hotspot": bool(pd.notna(qv) and qv <= qvalue_max),
        })
    out = pd.DataFrame(rows).drop_duplicates("pkey").reset_index(drop=True)
    out = out.sort_values("hotspot_alt_count", ascending=False).reset_index(drop=True)
    return out


def main():
    ap = argparse.ArgumentParser(description="cancerhotspots v2 -> residue hotspot positive labels.")
    root = _repo_root()
    ap.add_argument("--xls", default=os.path.join(root, "data", "external", "cancerhotspots", "hotspots_v2.xls"))
    ap.add_argument("--out", default=os.path.join(root, "data", "interim", "cancerhotspots_labels.tsv"))
    ap.add_argument("--qvalue-max", type=float, default=0.05)
    args = ap.parse_args()

    df = build_hotspot_labels(args.xls, args.qvalue_max)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print("[cancerhotspots] missense alt-specific rows: %d" % len(df))
    print("  is_hotspot (q<=%.3g): %d" % (args.qvalue_max, int(df["is_hotspot"].sum())))
    print("  unique genes: %d" % df["gene"].nunique())
    print("[cancerhotspots] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
