"""dbNSFP 5.3.1a (GRCh38, BGZF+tabix) -> 패널 missense 다도구 점수표.

설계 (PROGRESS.md §7, [[join-key-gene-protein]]):
  - dbNSFP BGZF는 GRCh38 좌표 정렬 + tabix 인덱스 → **패널 768 유전자 좌표영역만 tabix로 추출**(50GB 스트리밍 불요).
  - 점수 정규화: dbNSFP 제공 **`*_rankscore`(0–1, 높을수록 damaging, 방향 통일)** 컬럼을 그대로 사용 → 30여 도구 한 번에 비교 가능.
  - 보존도(GERP++_RS, phyloP/phastCons) + **gnomAD4.1 joint AF**(common variant=benign negative 보강).
  - 다중 transcript 값('; ' 구분)은 숫자=최댓값(가장 damaging), 결측='.'→ NaN.
  - 출력은 **GRCh38 좌표 키** → build_dataset 에서 AM universe 에 좌표 join.

⚠️ 컬럼은 이름으로 동적 선택(버전 견고). pysam.TabixFile 사용(CLI tabix 불요).
출력: data/interim/dbnsfp_panel.tsv
사용: python -m src.data.dbnsfp_join
"""
from __future__ import annotations
import argparse
import gzip
import os

import pandas as pd
import pysam

AA1 = set("ACDEFGHIKLMNPQRSTVWY")

KEY_COLS = ["#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt", "genename"]
# rankscore 외에 추가로 보존할 컬럼(원점수/보존도/빈도/메타)
EXTRA_COLS = [
    "AlphaMissense_score", "AlphaMissense_pred", "REVEL_score", "CADD_phred",
    "GERP++_RS", "phyloP100way_vertebrate", "phyloP470way_mammalian",
    "phastCons100way_vertebrate", "bStatistic",
    "gnomAD4.1_joint_AF", "gnomAD4.1_joint_POPMAX_AF",
    "clinvar_clnsig", "MANE",
]


def _repo_root() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _read_header(bgzf: str) -> list[str]:
    with gzip.open(bgzf, "rt", errors="replace") as f:
        return f.readline().rstrip("\n").split("\t")


def _num_max(cell: str):
    """다중 transcript('; ' 또는 ';') 숫자 셀 -> 최댓값(가장 damaging). '.'/빈값 -> None."""
    if cell is None or cell in (".", ""):
        return None
    best = None
    for tok in cell.replace(" ", "").split(";"):
        if tok in (".", ""):
            continue
        try:
            v = float(tok)
        except ValueError:
            continue
        best = v if best is None else max(best, v)
    return best


def _first_str(cell: str):
    if not cell or cell == ".":
        return ""
    for tok in cell.split(";"):
        if tok and tok != ".":
            return tok
    return ""


def build_dbnsfp_panel(bgzf: str, panel: pd.DataFrame):
    header = _read_header(bgzf)
    idx = {c: i for i, c in enumerate(header)}
    rank_cols = [c for c in header if c.endswith("_rankscore")]
    sel = [c for c in KEY_COLS + rank_cols + EXTRA_COLS if c in idx]
    numeric_cols = set(rank_cols) | {
        "AlphaMissense_score", "REVEL_score", "CADD_phred", "GERP++_RS",
        "phyloP100way_vertebrate", "phyloP470way_mammalian",
        "phastCons100way_vertebrate", "bStatistic",
        "gnomAD4.1_joint_AF", "gnomAD4.1_joint_POPMAX_AF"}
    print("[dbnsfp] rankscore 도구 %d개 + 추가 %d개 선택"
          % (len(rank_cols), len([c for c in EXTRA_COLS if c in idx])))

    # dbNSFP 일부 필드에 비ASCII 문자 → latin-1은 모든 바이트를 디코드(크래시 방지).
    # 우리가 쓰는 건 숫자 점수 컬럼이라 텍스트 글리프 손실은 무해.
    tbx = pysam.TabixFile(bgzf, encoding="latin-1")
    contigs = set(tbx.contigs)
    rows = []
    seen = set()
    n_genes_done = 0
    for _, g in panel.iterrows():
        ch = str(g["chromosome"]).replace("chr", "").strip()
        if ch not in contigs:
            continue
        try:
            start, stop = int(g["genome_start"]), int(g["genome_stop"])
        except (ValueError, TypeError):
            continue
        n_genes_done += 1
        if n_genes_done % 100 == 0:
            print("  ...%d genes, %d variants" % (n_genes_done, len(rows)), flush=True)
        for line in tbx.fetch(ch, max(0, start - 1), stop):
            p = line.split("\t")
            aaref, aaalt = p[idx["aaref"]], p[idx["aaalt"]]
            ref, alt = p[idx["ref"]], p[idx["alt"]]
            if len(ref) != 1 or len(alt) != 1:
                continue
            if aaref not in AA1 or aaalt not in AA1 or aaref == aaalt:
                continue  # 단일 missense 만
            gkey = (ch, p[idx["pos(1-based)"]], ref, alt)
            if gkey in seen:
                continue
            seen.add(gkey)
            rec = {
                "grch38_chr": ch, "grch38_pos": int(p[idx["pos(1-based)"]]),
                "grch38_ref": ref, "grch38_alt": alt,
                "aaref": aaref, "aaalt": aaalt,
                "dbnsfp_gene": _first_str(p[idx["genename"]]),
            }
            for c in sel:
                if c in KEY_COLS:
                    continue
                raw = p[idx[c]]
                key = "dbnsfp_" + c.replace("++", "pp").replace("(", "").replace(")", "")
                rec[key] = _num_max(raw) if c in numeric_cols else _first_str(raw)
            rows.append(rec)
    tbx.close()
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description="dbNSFP -> panel missense multi-tool scores (tabix).")
    root = _repo_root()
    ap.add_argument("--dbnsfp", default=os.path.join(
        root, "data", "external", "dbnsfp", "dbNSFP5.3.1a_grch38.gz"))
    ap.add_argument("--panel", default=os.path.join(root, "data", "interim", "cancer_gene_panel.tsv"))
    ap.add_argument("--out", default=os.path.join(root, "data", "interim", "dbnsfp_panel.tsv"))
    args = ap.parse_args()

    panel = pd.read_csv(args.panel, sep="\t")
    print("[dbnsfp] panel genes: %d" % len(panel))
    df = build_dbnsfp_panel(args.dbnsfp, panel)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print("[dbnsfp] panel missense variants: %d (genes=%d)"
          % (len(df), df["dbnsfp_gene"].nunique()))
    print("  컬럼 수: %d" % df.shape[1])
    # AlphaMissense 일치성 점검(브리지 vs dbNSFP)
    if "dbnsfp_AlphaMissense_score" in df.columns:
        print("  AlphaMissense_score 결측률: %.1f%%"
              % (df["dbnsfp_AlphaMissense_score"].isna().mean() * 100))
    print("[dbnsfp] wrote -> %s" % args.out)


if __name__ == "__main__":
    main()
