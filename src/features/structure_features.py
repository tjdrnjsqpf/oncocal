"""AlphaFold structure -> per-residue structure features (RSA/burial, pLDDT). Structure-based mechanism signal.

Hypothesis (research_plan §2, Part B residual): 'non-disruptive activation' (low ΔΔG) tends to occur at surface/functional sites.
  - **RSA** (relative solvent accessibility, freesasa relativeTotal 0–1): buried (low)=disruptive→tools catch it / surface (high)=non-disruptive→GOF possible.
  - **pLDDT** (AF confidence, CA B-factor): if low, disordered → RSA·ΔΔG less trustworthy (both a feature and a quality flag).
  - residue number = UniProt canonical → maps directly to the variant protein position (AM also uses UniProt numbering).

Output: data/interim/structure_features.tsv (uniprot, pos, aa3, plddt, rsa, buried)
Usage: python -m src.features.structure_features
"""
from __future__ import annotations
import glob
import os
import re
import warnings

import pandas as pd

warnings.filterwarnings("ignore")
import freesasa
freesasa.setVerbosity(freesasa.silent)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
AFDIR = os.path.join(ROOT, "data", "external", "alphafold")
FNAME = re.compile(r"AF-([A-Z0-9]+)-F\d+-model_v\d+\.pdb$")


def plddt_per_residue(pdb: str) -> dict:
    out = {}
    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    out[int(line[22:26])] = float(line[60:66])
                except ValueError:
                    pass
    return out


def process(pdb: str):
    uni = FNAME.search(os.path.basename(pdb)).group(1)
    plddt = plddt_per_residue(pdb)
    try:
        st = freesasa.Structure(pdb)
        ra = freesasa.calc(st).residueAreas()
    except Exception:
        return []
    rows = []
    for chain in ra:
        for posn, area in ra[chain].items():
            try:
                pos = int(posn)
            except ValueError:
                continue
            rsa = area.relativeTotal
            rows.append({
                "uniprot": uni, "pos": pos, "aa3": area.residueType,
                "plddt": plddt.get(pos),
                "rsa": round(rsa, 4) if rsa == rsa else None,
                "buried": int(rsa < 0.25) if rsa == rsa else None,
            })
        break  # single chain (AF monomer)
    return rows


STRUCT_COLS = ["struct_plddt", "struct_rsa", "struct_buried"]


def add_structure_features(df: pd.DataFrame) -> pd.DataFrame:
    """Join structure features to the variant df (gene, protein_change): gene->uniprot, protein pos -> (uniprot,pos)."""
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "data"))
    from variant_keys import parse_protein_change
    g2u = pd.read_csv(os.path.join(ROOT, "data", "interim", "gene_uniprot.tsv"), sep="\t")
    g2u = dict(zip(g2u["gene"].str.upper(), g2u["uniprot"].fillna("")))
    st = pd.read_csv(os.path.join(ROOT, "data", "interim", "structure_features.tsv"), sep="\t")
    st = st.rename(columns={"plddt": "struct_plddt", "rsa": "struct_rsa", "buried": "struct_buried"})
    st = st[["uniprot", "pos", "struct_plddt", "struct_rsa", "struct_buried"]]

    df = df.copy()
    df["_uni"] = df["gene"].astype(str).str.upper().map(g2u)
    df["_pos"] = df["protein_change"].astype(str).map(
        lambda s: (parse_protein_change(s) or (None, None, None))[1])
    out = df.merge(st.rename(columns={"uniprot": "_uni", "pos": "_pos"}),
                   on=["_uni", "_pos"], how="left").drop(columns=["_uni", "_pos"])
    return out


def main():
    pdbs = sorted(glob.glob(os.path.join(AFDIR, "AF-*-F1-model_v*.pdb")))
    print("[struct_feat] processing %d AF structures..." % len(pdbs))
    allrows = []
    for i, pdb in enumerate(pdbs):
        allrows.extend(process(pdb))
        if (i + 1) % 100 == 0:
            print("  ...%d/%d (%d residues)" % (i + 1, len(pdbs), len(allrows)), flush=True)
    df = pd.DataFrame(allrows)
    out = os.path.join(ROOT, "data", "interim", "structure_features.tsv")
    df.to_csv(out, sep="\t", index=False)
    print("[struct_feat] residues=%d, uniprots=%d, buried rate=%.2f, pLDDT median=%.1f"
          % (len(df), df["uniprot"].nunique(), df["buried"].dropna().mean(),
             df["plddt"].dropna().median()))
    print("[struct_feat] -> %s" % out)


if __name__ == "__main__":
    main()
