#!/usr/bin/env bash
# COSMIC scripted download (학술 계정 필요). 서버에서 실행.
# 비밀번호는 코드에 넣지 말고 환경변수로:
#   export COSMIC_EMAIL='you@institution.edu'
#   export COSMIC_PASS='********'        # (반드시 변경한 새 비밀번호)
#   bash scripts/server/03_download_cosmic.sh
#
# ⚠️ COSMIC API 경로/버전은 release마다 바뀜. 정확한 path는 COSMIC download 페이지에서
#    각 파일의 "scripted download" 예시를 확인해 아래 FILES 배열을 맞출 것.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$ROOT/data/external/cosmic"; mkdir -p "$OUT"
: "${COSMIC_EMAIL:?COSMIC_EMAIL 환경변수 필요}"
: "${COSMIC_PASS:?COSMIC_PASS 환경변수 필요}"

VER="v100"        # <-- 최신 버전으로 수정
BUILD="GRCh38"
AUTH=$(printf '%s:%s' "$COSMIC_EMAIL" "$COSMIC_PASS" | base64)

# 받을 제품 파일들 (페이지의 scripted-download 경로와 일치시킬 것)
FILES=(
  "Cosmic_CancerGeneCensus_Tsv_${VER}_${BUILD}.tar"
  "Cosmic_MutantCensus_Tsv_${VER}_${BUILD}.tar"
)

API="https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted"
for f in "${FILES[@]}"; do
  echo "[req] $f"
  # 1) 서명된 다운로드 URL 요청 (JSON: {"url": "..."} 반환)
  url=$(curl -s -H "Authorization: Basic ${AUTH}" \
        "${API}?path=${BUILD}/cosmic/${VER}/${f}&bucket=downloads" \
        | sed -E 's/.*"url"\s*:\s*"([^"]+)".*/\1/')
  if [[ -z "$url" || "$url" == *'<'* ]]; then
    echo "  ! URL 못 받음. download 페이지에서 정확한 path/버전 확인 필요."
    continue
  fi
  # 2) 실제 다운로드
  curl -L "$url" -o "$OUT/$f"
  echo "  [ok] $OUT/$f"
done

echo "[done] COSMIC. 받은 파일:"; ls -lh "$OUT"
