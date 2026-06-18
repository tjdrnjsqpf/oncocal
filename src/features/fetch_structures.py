"""패널 유전자 -> UniProt 매핑 + AlphaFold 구조(v6 pdb) 다운로드.

구조 기반 기전 feature(ΔΔG·RSA/burial·pLDDT·3D clustering)의 공통 전제.
  - gene->uniprot: CGC 패널 uniprot(720) + 결손(BRAF 등 48개)은 dbNSFP gene table(Uniprot_acc)로 보강.
  - AlphaFold DB 최신 = model_v6. API로 정확 URL 확인 후 pdb 다운로드.
  - >2700 잔기 단백질은 F1만(조각 분할) — 한계 기록.

출력: data/interim/gene_uniprot.tsv, data/external/alphafold/AF-<uniprot>-F1-model_v6.pdb
사용: python -m src.features.fetch_structures
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
    # dbNSFP gene table: Gene_name -> Uniprot_acc (결손 보강)
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
    print("[struct] gene->uniprot: %d/%d 매핑" % (nmap, len(df)))
    ok, skip, fail, fails = fetch(df)
    print("[struct] AlphaFold pdb: ok=%d skip=%d fail=%d" % (ok, skip, fail))
    if fails:
        print("  실패 uniprot(앞 15):", fails[:15])
    print("[struct] -> data/external/alphafold/ , data/interim/gene_uniprot.tsv")


if __name__ == "__main__":
    main()
