#!/usr/bin/env bash
# 등록 불필요한 공개 데이터 다운로드 (서버에서 실행)
# 사용법: repo 루트에서  bash scripts/server/02_download_open_data.sh
# 주: 일부 URL은 버전에 따라 바뀜 → 아래 변수에서 최신 링크 확인/수정.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
EXT="$ROOT/data/external"
log(){ echo "[$(date +%H:%M:%S)] $*"; }

# ---------------------------------------------------------------------------
# 1) ClinVar (NCBI, public domain) — 안정적 직링크
# ---------------------------------------------------------------------------
log "ClinVar variant_summary..."
wget -c -P "$EXT/clinvar" \
  "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

# ---------------------------------------------------------------------------
# 2) CIViC (CC0) — nightly TSV (안정 경로). 월별 release를 원하면 releases 페이지 URL로 교체.
# ---------------------------------------------------------------------------
log "CIViC nightly TSVs..."
CIVIC_BASE="https://civicdb.org/downloads/nightly"
for f in VariantSummaries MolecularProfileSummaries ClinicalEvidenceSummaries \
         AssertionSummaries FeatureSummaries VariantGroupSummaries; do
  wget -c -P "$EXT/civic" "${CIVIC_BASE}/nightly-${f}.tsv" || \
    log "  ! ${f} 실패 — https://civicdb.org/releases 에서 직접 확인 필요"
done

# ---------------------------------------------------------------------------
# 3) cancerhotspots.org (등록 불요) — 직링크가 막히면 페이지에서 링크 복사해 교체
# ---------------------------------------------------------------------------
log "cancerhotspots v2..."
wget -c -P "$EXT/cancerhotspots" \
  "https://www.cancerhotspots.org/files/hotspots_v2.xls" || \
  log "  ! 실패 — https://www.cancerhotspots.org/#/download 에서 링크 확인"

# ---------------------------------------------------------------------------
# 4) dbNSFP (academic, 무료) — ⚠️ 직링크가 페이지/버전마다 달라 변수로 지정
#    https://www.dbnsfp.org/download 에서 academic('a') zip 링크를 복사해 아래에 넣기.
# ---------------------------------------------------------------------------
DBNSFP_URL=""   # <-- 예: "https://.../dbNSFP4.9a.zip"  (academic 'a' 버전!)
if [[ -n "$DBNSFP_URL" ]]; then
  log "dbNSFP download (대용량 ~30GB+, 시간 오래 걸림)..."
  wget -c -P "$EXT/dbnsfp" "$DBNSFP_URL"
else
  log "dbNSFP: DBNSFP_URL 비어있음 → 스크립트 상단 변수에 academic zip 링크 채운 뒤 재실행."
fi

log "open-data 완료. 받은 파일:"
find "$EXT" -type f | sort
