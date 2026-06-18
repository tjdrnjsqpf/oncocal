"""Variant key utilities — shared functions to tie all sources together with the same join key.

Key decision (research_plan.md §6, phase1_execution.md C):
  - CIViC/COSMIC/cancerhotspots/dbNSFP **use different genome coordinate builds** (CIViC=GRCh37).
  - Therefore the primary join key is not genome coordinates but **gene + protein change**:
        pkey = "<GENE>:p.<REF_AA><POS><ALT_AA>"   (1-letter, e.g. "BRAF:p.V600E")
  - Coordinate-based keys (chr,pos,ref,alt) are used only as an auxiliary key between same-build sources such as dbNSFP↔COSMIC.

This module uses only the standard library, with no external dependencies (parser reproducibility).
"""
from __future__ import annotations
import re

# 3-letter -> 1-letter amino acid (including Ter/*, and Sec/Pyl)
AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*", "Xaa": "X",
}
AA1 = set("ACDEFGHIKLMNPQRSTVWY")  # standard 20 (for missense determination)

# p.Val600Glu / p.V600E / Val600Glu / V600E  (single-residue substitution)
_RE_P3 = re.compile(r"p?\.?\s*([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\b")
_RE_P1 = re.compile(r"p?\.?\s*([A-Z])(\d+)([A-Z*])\b")


def aa_to1(aa: str) -> str | None:
    """Convert an amino acid notation to 1-letter. If already 1-letter, keep it. None if unknown."""
    if not aa:
        return None
    aa = aa.strip()
    if len(aa) == 1:
        return aa.upper()
    return AA3TO1.get(aa.capitalize())


def parse_protein_change(s: str):
    """Extract a single-residue substitution (ref_aa, pos, alt_aa) from a string as 1-letter.

    Handles 'V600E', 'p.V600E', 'p.Val600Glu', 'NP_004324.2:p.Val600Glu', etc.
    None on match failure. (frameshift/indel/fusion/splice etc. return None)
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
    """Whether it is a single missense: ref/alt both standard 20 AA, different from each other, not nonsense."""
    if not ref_aa or not alt_aa:
        return False
    return ref_aa in AA1 and alt_aa in AA1 and ref_aa != alt_aa


def protein_key(gene: str, ref_aa: str, pos: int, alt_aa: str) -> str:
    """Standard protein-level join key: 'GENE:p.<REF><POS><ALT>' (uppercase gene)."""
    return "%s:p.%s%d%s" % (str(gene).strip().upper(), ref_aa, int(pos), alt_aa)


def make_pkey(gene: str, change: str):
    """gene + arbitrary protein-change string -> standard pkey. None if unparsed/non-missense."""
    pc = parse_protein_change(change)
    if pc is None or not gene:
        return None
    ref, pos, alt = pc
    if not is_missense_change(ref, alt):
        return None
    return protein_key(gene, ref, pos, alt)
