"""변이 키 유틸 — 모든 소스를 동일한 join 키로 묶기 위한 공용 함수.

핵심 결정 (research_plan.md §6, phase1_execution.md C):
  - CIViC/COSMIC/cancerhotspots/dbNSFP의 **유전체 좌표 빌드가 제각각**(CIViC=GRCh37).
  - 따라서 1차 join 키는 유전체 좌표가 아니라 **gene + protein change**:
        pkey = "<GENE>:p.<REF_AA><POS><ALT_AA>"   (1-letter, 예: "BRAF:p.V600E")
  - 좌표 기반 키(chr,pos,ref,alt)는 dbNSFP↔COSMIC 등 동일 빌드 소스 간 보조 키로만 사용.

이 모듈은 외부 의존성 없이 표준 라이브러리만 사용한다(파서 재현성).
"""
from __future__ import annotations
import re

# 3-letter -> 1-letter amino acid (Ter/*, 그리고 Sec/Pyl 포함)
AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*", "Xaa": "X",
}
AA1 = set("ACDEFGHIKLMNPQRSTVWY")  # 표준 20종 (missense 판정용)

# p.Val600Glu / p.V600E / Val600Glu / V600E  (단일 잔기 치환)
_RE_P3 = re.compile(r"p?\.?\s*([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\b")
_RE_P1 = re.compile(r"p?\.?\s*([A-Z])(\d+)([A-Z*])\b")


def aa_to1(aa: str) -> str | None:
    """아미노산 표기를 1-letter로. 이미 1-letter면 그대로. 모르면 None."""
    if not aa:
        return None
    aa = aa.strip()
    if len(aa) == 1:
        return aa.upper()
    return AA3TO1.get(aa.capitalize())


def parse_protein_change(s: str):
    """문자열에서 단일 잔기 치환 (ref_aa, pos, alt_aa) 를 1-letter로 추출.

    'V600E', 'p.V600E', 'p.Val600Glu', 'NP_004324.2:p.Val600Glu' 등 처리.
    매칭 실패 시 None. (frameshift/indel/fusion/splice 등은 None)
    """
    if not s or not isinstance(s, str):
        return None
    m = _RE_P3.search(s)
    if m:
        ref, pos, alt = AA3TO1.get(m.group(1)), m.group(2), AA3TO1.get(m.group(3))
    else:
        m = _RE_P1.search(s)
        if not m:
            return None
        ref, pos, alt = m.group(1), m.group(2), m.group(3)
    if ref is None or alt is None:
        return None
    return ref, int(pos), alt


def is_missense_change(ref_aa: str | None, alt_aa: str | None) -> bool:
    """단일 missense 여부: ref/alt 모두 표준 20 AA, 서로 다름, nonsense 아님."""
    if not ref_aa or not alt_aa:
        return False
    return ref_aa in AA1 and alt_aa in AA1 and ref_aa != alt_aa


def protein_key(gene: str, ref_aa: str, pos: int, alt_aa: str) -> str:
    """표준 protein-level join 키: 'GENE:p.<REF><POS><ALT>' (대문자 gene)."""
    return "%s:p.%s%d%s" % (str(gene).strip().upper(), ref_aa, int(pos), alt_aa)


def make_pkey(gene: str, change: str):
    """gene + 임의 protein-change 문자열 -> 표준 pkey. 미파싱/비missense면 None."""
    pc = parse_protein_change(change)
    if pc is None or not gene:
        return None
    ref, pos, alt = pc
    if not is_missense_change(ref, alt):
        return None
    return protein_key(gene, ref, pos, alt)
