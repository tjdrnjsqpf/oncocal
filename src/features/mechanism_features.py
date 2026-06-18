"""기전 시그니처 feature — 서열 기반 생화학 변화 (구조 불요, Phase 2 차별점 1단계).

가설(research_plan §2, Part B 잔여 클래스):
  - germline 도구는 '보존+구조 파괴'를 봄 → **온건한(낮은 ΔΔG 경향) 활성화 치환**을 놓침.
    서열 수준 proxy: 작은 소수성/부피/극성 변화 = 비파괴적 = GOF 신호일 수 있음 (예: JAK2 V617F, V→F).
  - **neomorphic cysteine 생성** = 비정상 disulfide 이량체화(예: FGFR3 S249C/Y373C) → cysteine-gain 플래그.
  - proline 도입(구조 꺾임), glycine 관여(유연성), 전하 역전 등.

구조 기반(ΔΔG·activation domain·3D clustering)은 다음 단계(AlphaFold). 여기선 ref/alt AA만 사용.
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd

try:
    from ..data.variant_keys import parse_protein_change
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "data"))
    from variant_keys import parse_protein_change

# Kyte-Doolittle hydropathy
HYDRO = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
         "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
         "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# 잔기 부피 (Å^3)
VOL = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4,
       "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9,
       "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0}
# 전하 (생리적 pH)
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}
# Grantham polarity
POLAR = {"A": 8.1, "R": 10.5, "N": 11.6, "D": 13.0, "C": 5.5, "Q": 10.5, "E": 12.3,
         "G": 9.0, "H": 10.4, "I": 5.2, "L": 4.9, "K": 11.3, "M": 5.7, "F": 5.0,
         "P": 8.0, "S": 9.2, "T": 8.6, "W": 5.4, "Y": 6.2, "V": 5.9}


def _charge(a):
    return CHARGE.get(a, 0.0)


def features_for(ref: str, alt: str) -> dict:
    d = {
        "mech_d_hydropathy": HYDRO[alt] - HYDRO[ref],
        "mech_abs_d_hydropathy": abs(HYDRO[alt] - HYDRO[ref]),
        "mech_d_volume": VOL[alt] - VOL[ref],
        "mech_abs_d_volume": abs(VOL[alt] - VOL[ref]),
        "mech_d_polarity": POLAR[alt] - POLAR[ref],
        "mech_d_charge": _charge(alt) - _charge(ref),
        "mech_abs_d_charge": abs(_charge(alt) - _charge(ref)),
        "mech_cysteine_gain": int(alt == "C"),
        "mech_cysteine_loss": int(ref == "C"),
        "mech_proline_gain": int(alt == "P"),
        "mech_proline_loss": int(ref == "P"),
        "mech_glycine_involved": int(ref == "G" or alt == "G"),
        "mech_to_aromatic": int(alt in "FWY"),
    }
    # Grantham 거리 (조성 무시한 polarity+volume 근사)
    d["mech_grantham_like"] = float(np.sqrt(
        0.1018 * (POLAR[alt] - POLAR[ref]) ** 2 + 0.000399 * (VOL[alt] - VOL[ref]) ** 2))
    return d


MECH_COLS = list(features_for("V", "F").keys())


def add_mechanism_features(df: pd.DataFrame) -> pd.DataFrame:
    """df에 mech_* 컬럼 추가. protein_change(예 'V600E') 또는 ref_aa/alt_aa 사용."""
    recs = []
    for pc in df["protein_change"].astype(str):
        p = parse_protein_change(pc)
        if p is None:
            recs.append({c: np.nan for c in MECH_COLS})
            continue
        ref, _, alt = p
        try:
            recs.append(features_for(ref, alt))
        except KeyError:
            recs.append({c: np.nan for c in MECH_COLS})
    mech = pd.DataFrame(recs, index=df.index)
    return pd.concat([df, mech], axis=1)
