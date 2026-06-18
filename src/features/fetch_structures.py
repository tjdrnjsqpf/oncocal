"""Panel gene -> UniProt mapping + AlphaFold structure (v6 pdb) download.

Common prerequisite for structure-based mechanism features (ΔΔG, RSA/burial, pLDDT, 3D clustering).
  - gene->uniprot: CGC panel uniprot (720) + missing entries (BRAF etc., 48) supplemented via the dbNSFP gene table (Uniprot_acc).
  - Latest AlphaFold DB = model_v6. Confirm the exact URL via the API, then download the pdb.
  - For proteins >2700 residues, only F1 (fragment split) — limitation noted.

Output: data/interim/gene_uniprot.tsv, data/external/alphafold/AF-<uniprot>-F1-model_v6.pdb
Usage: python -m src.features.fetch_structures
"""
from __future__ import annotations
import gzip
import os
import time
import urllib.request

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
AFDIR = os.path.join(ROOT, "data", "external", "alphafold")


def build_gene_uniprot() -> pd.DataFrame:
    panel = pd.read_csv(os.path.join(ROOT, "data", "interim", "cancer_gene_panel.tsv"), sep="\t")
    # dbNSFP gene table: Gene_name -> Uniprot_acc (supplement missing entries)
    g2u = {}
    with gzip.open(os.path.join(ROOT, "data", "external", "dbnsfp", "dbNSFP5.3_gene.gz"), "rt",
                   errors="replace") as f:
        head = f.readline().rstrip("\n").split("\t")
        gi, ui = head.index("Gene_name"), head.index("Uniprot_acc")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= ui:
                continue
            u = (p[ui] or "").split(";")[0].strip()
            if p[gi] and u and u != ".":
                g2u[p[gi].upper()] = u
    rows = []
    for _, r in panel.iterrows():
        gene = str(r["gene"]).upper()
        uni = r["uniprot"] if isinstance(r["uniprot"], str) and len(str(r["uniprot"])) > 1 else ""
        if not uni:
            uni = g2u.get(gene, "")
        rows.append({"gene": gene, "uniprot": uni})
    df = pd.DataFrame(rows)
    return df


def fetch(df: pd.DataFrame):
    os.makedirs(AFDIR, exist_ok=True)
    ok = fail = skip = 0
    fails = []
    for i, u in enumerate(sorted(set(df["uniprot"]) - {""})):
        out = os.path.join(AFDIR, "AF-%s-F1-model_v6.pdb" % u)
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            skip += 1; continue
        url = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v6.pdb" % u
        try:
            urllib.request.urlretrieve(url, out)
            ok += 1
        except Exception:
            fail += 1; fails.append(u)
            if os.path.exists(out):
                os.remove(out)
        if (i + 1) % 50 == 0:
            print("  ...%d uniprots, ok=%d skip=%d fail=%d" % (i + 1, ok, skip, fail), flush=True)
        time.sleep(0.05)
    return ok, skip, fail, fails


def main():
    df = build_gene_uniprot()
    df.to_csv(os.path.join(ROOT, "data", "interim", "gene_uniprot.tsv"), sep="\t", index=False)
    nmap = int((df["uniprot"] != "").sum())
    print("[struct] gene->uniprot: %d/%d mapped" % (nmap, len(df)))
    ok, skip, fail, fails = fetch(df)
    print("[struct] AlphaFold pdb: ok=%d skip=%d fail=%d" % (ok, skip, fail))
    if fails:
        print("  failed uniprot (first 15):", fails[:15])
    print("[struct] -> data/external/alphafold/ , data/interim/gene_uniprot.tsv")


if __name__ == "__main__":
    main()
