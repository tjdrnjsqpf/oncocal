"""CIViC -> oncogenicity + 기전(GOF/LOF) 라벨 추출 (1차 학습 라벨, CC0).

설계 근거 (data_sources.md, phase1_execution.md C):
  - 라벨 backbone = CIViC(CC0). OncoKB는 학습 금지 → 미사용.
  - CIViC 좌표는 GRCh37 → join 키는 **gene + protein change** (variant_keys.protein_key).
  - oncogenicity = evidence_type == 'Oncogenic'  (significance='Oncogenicity')
        evidence_direction: 'Supports' = oncogenic 지지, 'Does Not Support' = 반증.
  - 기전(mechanism) = evidence_type == 'Functional' 의 significance:
        'Gain of Function' -> GOF,  'Loss of Function' -> LOF,  'Dominant Negative' -> DN.
  - variant_origin 으로 Somatic 위주 필터(germline-only 근거는 제외 가능).

⚠️ 핵심 함정(research_plan §6): CIViC oncogenic missense 단일변이는 수십 개 수준으로 작다.
   → 이 라벨은 '고신뢰 positive 코어 + 기전 소스'이며, positive 백본은 cancerhotspots/COSMIC와
     build_dataset.py 에서 결합한다. unknown 을 negative 로 쓰지 않는다.

출력: data/interim/civic_labels.tsv  (pkey 단위, 1행=1 protein-level 변이)
사용: python -m src.data.civic_labels  (repo 루트에서)
"""
from __future__ import annotations
import argparse
import os
import sys
from collections import defaultdict

import pandas as pd

# 패키지/스크립트 양쪽 실행 지원
try:
    from .variant_keys import parse_protein_change, is_missense_change, protein_key
except ImportError:  # python src/data/civic_labels.py 직접 실행 시
    sys.path.insert(0, os.path.dirname(__file__))
    from variant_keys import parse_protein_change, is_missense_change, protein_key

MECH_MAP = {
    "Gain of Function": "GOF",
    "Loss of Function": "LOF",
    "Dominant Negative": "DN",
}
SOMATIC_ORIGINS = {"Somatic", "Mixed", "Combined", "Unknown", "N/A"}  # germline-only 제외


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _resolve_pchange(row) -> str | None:
    """변이의 protein change 문자열 결정: variant 컬럼 우선, 실패 시 hgvs_descriptions."""
    for cand in (row.get("variant"), row.get("hgvs_descriptions")):
        if isinstance(cand, str) and cand.strip():
            if parse_protein_change(cand) is not None:
                return cand
    return None


def build_civic_labels(civic_dir: str) -> pd.DataFrame:
    vs = pd.read_csv(os.path.join(civic_dir, "nightly-VariantSummaries.tsv"),
                     sep="\t", dtype=str).fillna("")
    ev = pd.read_csv(os.path.join(civic_dir, "nightly-ClinicalEvidenceSummaries.tsv"),
                     sep="\t", dtype=str).fillna("")

    # 1) single-variant molecular profile 기준 variant <-> MP 매핑
    vs = vs[vs["single_variant_molecular_profile_id"] != ""].copy()
    vs["mp"] = vs["single_variant_molecular_profile_id"]

    # 2) missense + 단일 잔기 치환만 남기고 pkey 부여
    keep = []
    for _, r in vs.iterrows():
        pch = _resolve_pchange(r)
        if pch is None:
            continue
        pc = parse_protein_change(pch)
        if pc is None:
            continue
        ref, pos, alt = pc
        if not is_missense_change(ref, alt):
            continue
        keep.append({
            "mp": r["mp"], "civic_variant_id": r["variant_id"],
            "gene": r["gene"], "ref_aa": ref, "pos": pos, "alt_aa": alt,
            "protein_change": "%s%d%s" % (ref, pos, alt),
            "pkey": protein_key(r["gene"], ref, pos, alt),
            "chromosome": r["chromosome"], "start": r["start"],
            "reference_build": r["reference_build"],
            "representative_transcript": r["representative_transcript"],
            "civic_url": r["variant_civic_url"],
        })
    vmap = pd.DataFrame(keep)
    if vmap.empty:
        raise SystemExit("CIViC: missense 단일변이를 찾지 못함 — 입력 파일 확인 필요.")

    # 3) 근거(evidence) 집계: MP 단위로 oncogenicity / mechanism 카운트
    onc_sup = defaultdict(int); onc_ref = defaultdict(int)
    mech_sup = defaultdict(lambda: defaultdict(int))  # mp -> {GOF/LOF/DN: n}
    ev_types = defaultdict(set)
    for _, e in ev.iterrows():
        mp = e["molecular_profile_id"]
        if e["variant_origin"] and e["variant_origin"] not in SOMATIC_ORIGINS:
            # 'Rare Germline'/'Common Germline' 만 있는 근거는 oncogenicity/기전에서 제외
            germline_only = True
        else:
            germline_only = False
        et, direction, sig = e["evidence_type"], e["evidence_direction"], e["significance"]
        ev_types[mp].add(et)
        if et == "Oncogenic" and not germline_only:
            if direction == "Supports":
                onc_sup[mp] += 1
            elif direction == "Does Not Support":
                onc_ref[mp] += 1
        elif et == "Functional" and direction == "Supports" and sig in MECH_MAP:
            mech_sup[mp][MECH_MAP[sig]] += 1

    # 4) MP 집계를 variant pkey에 결합
    rows = []
    for _, v in vmap.iterrows():
        mp = v["mp"]
        ns, nr = onc_sup.get(mp, 0), onc_ref.get(mp, 0)
        mech_counts = mech_sup.get(mp, {})
        # 기전 결정: 최다 득표(동률이면 우선순위 GOF>LOF>DN). DN은 LOF 계열로도 해석되나 별도 보존.
        mech = None
        if mech_counts:
            mech = sorted(mech_counts.items(), key=lambda kv: (-kv[1], ["GOF", "LOF", "DN"].index(kv[0])))[0][0]
        rows.append({
            **{k: v[k] for k in ["pkey", "gene", "ref_aa", "pos", "alt_aa",
                                 "protein_change", "chromosome", "start",
                                 "reference_build", "representative_transcript",
                                 "civic_variant_id", "civic_url"]},
            "civic_onc_support": ns,
            "civic_onc_refute": nr,
            "civic_oncogenic": bool(ns > 0 and ns >= nr),  # 순지지 시 oncogenic
            "civic_mechanism": mech,
            "civic_mech_GOF": mech_counts.get("GOF", 0),
            "civic_mech_LOF": mech_counts.get("LOF", 0),
            "civic_mech_DN": mech_counts.get("DN", 0),
            "civic_evidence_types": "|".join(sorted(ev_types.get(mp, set()))),
        })
    out = pd.DataFrame(rows)

    # 동일 pkey가 여러 transcript/variant_id로 중복될 수 있음 → 카운트 합산해 1행으로
    agg = {
        "gene": "first", "ref_aa": "first", "pos": "first", "alt_aa": "first",
        "protein_change": "first", "chromosome": "first", "start": "first",
        "reference_build": "first", "representative_transcript": "first",
        "civic_variant_id": lambda s: "|".join(sorted(set(s))),
        "civic_url": "first",
        "civic_onc_support": "sum", "civic_onc_refute": "sum",
        "civic_oncogenic": "any",
        "civic_mech_GOF": "sum", "civic_mech_LOF": "sum", "civic_mech_DN": "sum",
        "civic_evidence_types": lambda s: "|".join(sorted(set("|".join(s).split("|")) - {""})),
    }
    out = out.groupby("pkey", as_index=False).agg(agg)
    # 합산 후 기전 재결정
    def _mech(r):
        c = {"GOF": r["civic_mech_GOF"], "LOF": r["civic_mech_LOF"], "DN": r["civic_mech_DN"]}
        c = {k: v for k, v in c.items() if v > 0}
        if not c:
            return None
        return sorted(c.items(), key=lambda kv: (-kv[1], ["GOF", "LOF", "DN"].index(kv[0])))[0][0]
    out["civic_mechanism"] = out.apply(_mech, axis=1)
    out = out.sort_values(["civic_oncogenic", "civic_onc_support"], ascending=False).reset_index(drop=True)
    return out


def main():
    ap = argparse.ArgumentParser(description="Extract CIViC oncogenicity/mechanism labels (gene+protein keyed).")
    root = _repo_root()
    ap.add_argument("--civic-dir", default=os.path.join(root, "data", "external", "civic"))
    ap.add_argument("--out", default=os.path.join(root, "data", "interim", "civic_labels.tsv"))
    args = ap.parse_args()

    df = build_civic_labels(args.civic_dir)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    n = len(df)
    n_onc = int(df["civic_oncogenic"].sum())
    n_mech = int(df["civic_mechanism"].notna().sum())
    print("[civic_labels] missense pkey rows: %d" % n)
    print("  oncogenic (net support): %d" % n_onc)
    print("  with mechanism (GOF/LOF/DN): %d" % n_mech)
    if n_mech:
        print(df.loc[df["civic_mechanism"].notna(), "civic_mechanism"].value_counts().to_string())
    print("[civic_labels] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
