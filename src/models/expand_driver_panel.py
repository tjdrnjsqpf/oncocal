"""canonical_driver_panel.tsv 재구성+확장 (idempotent, 행 추가만, 새 학습 없음).
oncogene ~30 + TSG ~10 + both(TP53/ESR1 등) 로 확대. 선정 = 라이선스-청정 positive(hotspot|CIViC|COSMIC rec≥10).
점수는 dataset.parquet/dbNSFP에서 그대로. 없는 값 빈칸. 유전자당 ≤3 변이(다양성). 선정 근거 출력.
사용: python -m src.models.expand_driver_panel
"""
from __future__ import annotations
import os
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PANEL = os.path.join(ROOT, "paper_handoff", "tables", "canonical_driver_panel.tsv")

COLMAP = {
    "dbnsfp_AlphaMissense_rankscore": "AlphaMissense", "dbnsfp_REVEL_rankscore": "REVEL",
    "dbnsfp_CADD_raw_rankscore": "CADD_raw", "dbnsfp_ESM1b_converted_rankscore": "ESM1b",
    "dbnsfp_popEVE_converted_rankscore": "popEVE", "dbnsfp_MetaRNN_rankscore": "MetaRNN",
    "dbnsfp_GERPpp_RS": "GERPpp_RS", "dbnsfp_phyloP100way_vertebrate": "phyloP100way_vertebrate",
}
OUTCOLS = ["pkey", "role_class", "cosmic_sample_count", "am_pathogenicity"] + list(COLMAP.values())

# 원본 17 (재현용 base) — part_b b5 canonical 패널
BASE = ["JAK2:p.V617F", "BRAF:p.V600E", "KRAS:p.G12D", "IDH1:p.R132H", "PIK3CA:p.H1047R",
        "EGFR:p.L858R", "NRAS:p.Q61K", "FGFR3:p.S249C", "KIT:p.D816V", "FGFR3:p.Y373C",
        "FLT3:p.D835Y", "GNA11:p.Q209L", "GNAQ:p.Q209L", "MPL:p.W515L", "ERBB2:p.R678Q",
        "STAT3:p.Y640F", "KIT:p.V559D"]
CURATED = [
    "KRAS:p.G12V", "KRAS:p.G12C", "KRAS:p.G13D", "KRAS:p.Q61H", "KRAS:p.G12R",
    "NRAS:p.Q61R", "NRAS:p.G12D", "HRAS:p.G12V", "HRAS:p.Q61R",
    "BRAF:p.G469A", "BRAF:p.K601E", "BRAF:p.G466V", "PIK3CA:p.E545K", "PIK3CA:p.E542K",
    "PIK3CA:p.H1047L", "AKT1:p.E17K", "CTNNB1:p.S45F", "CTNNB1:p.T41A", "CTNNB1:p.S33Y",
    "ESR1:p.Y537S", "ESR1:p.D538G", "RAC1:p.P29S", "SF3B1:p.K700E", "IDH1:p.R132C",
    "IDH2:p.R140Q", "IDH2:p.R172K", "GNAS:p.R201C", "GNAS:p.R201H", "EGFR:p.G719S",
    "EGFR:p.T790M", "EGFR:p.L861Q", "ERBB2:p.S310F", "ERBB2:p.V842I", "MAP2K1:p.K57N",
    "MYD88:p.L265P", "FGFR3:p.R248C", "U2AF1:p.S34F",
    "TP53:p.R175H", "TP53:p.R248Q", "TP53:p.R273H", "TP53:p.R282W", "TP53:p.Y220C",
    "PTEN:p.R130Q", "PTEN:p.R130G", "PTEN:p.C124S", "FBXW7:p.R465C", "FBXW7:p.R505C",
    "SMAD4:p.R361H", "SMAD4:p.R361C", "VHL:p.L178P", "KEAP1:p.R320Q", "CDKN2A:p.D108Y",
    "STK11:p.F354L",
]


def main():
    df = pd.read_parquet(os.path.join(ROOT, "data", "processed", "dataset.parquet")).drop_duplicates("pkey")
    pos = df[df["oncogenic_label"] == "positive"].copy()
    posset = set(pos["pkey"])
    GENE_CAP = 3
    gene_n, chosen, why = {}, [], []

    def add(pk, src):
        if pk in chosen or pk not in posset:
            return False
        g = pk.split(":")[0]
        if gene_n.get(g, 0) >= GENE_CAP:
            return False
        r = pos[pos.pkey == pk].iloc[0]
        chosen.append(pk); gene_n[g] = gene_n.get(g, 0) + 1
        why.append((pk, r.role_class, int(r.cosmic_sample_count), bool(r.is_hotspot),
                    bool(r.civic_oncogenic), src))
        return True

    def role_total(role):
        return sum(1 for w in why if w[1] == role)

    for pk in BASE:
        add(pk, "base17")
    for pk in CURATED:
        add(pk, "curated")
    for role, tgt in [("oncogene", 30), ("TSG", 10)]:
        for _, r in pos[pos.role_class == role].sort_values("cosmic_sample_count", ascending=False).iterrows():
            if role_total(role) >= tgt:
                break
            add(r.pkey, "recurrence-autofill")

    sub = df[df.pkey.isin(chosen)].rename(columns=COLMAP)
    out = sub[[c for c in OUTCOLS if c in sub.columns]].copy()
    out = out.sort_values(["role_class", "cosmic_sample_count"], ascending=[True, False]).reset_index(drop=True)
    out.to_csv(PANEL, sep="\t", index=False)

    print("총 %d행 (base %d + 확장)" % (len(out), len(BASE)))
    print("role_class 분포:", out.role_class.value_counts().to_dict())
    print("  (주의: TP53·ESR1은 CGC에서 role='both' — mutant p53 GOF 때문. heatmap에서 oncogene/both/TSG 3그룹 권장)")
    print("\n변이 선정 근거(role/recurrence/signal):")
    for pk, role, n, hot, civ, src in sorted(why, key=lambda x: (x[1], -x[2])):
        sig = ",".join([s for s, b in [("hotspot", hot), ("CIViC", civ)] if b]) or "rec"
        print("  %-18s %-9s rec=%-6d [%s] (%s)" % (pk, role, n, sig, src))
    print("\n[done] -> paper_handoff/tables/canonical_driver_panel.tsv")


if __name__ == "__main__":
    main()
