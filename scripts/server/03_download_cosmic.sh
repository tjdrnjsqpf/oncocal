#!/usr/bin/env bash
# COSMIC scripted download (academic account required). Run on the server.
# Do not put the password in code; use environment variables:
#   export COSMIC_EMAIL='you@institution.edu'
#   export COSMIC_PASS='********'        # (always a freshly changed password)
#   bash scripts/server/03_download_cosmic.sh
#
# ⚠️ COSMIC API paths/versions change every release. For the exact path, check each file's
#    "scripted download" example on the COSMIC download page and match the FILES array below.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$ROOT/data/external/cosmic"; mkdir -p "$OUT"
: "${COSMIC_EMAIL:?COSMIC_EMAIL environment variable required}"
: "${COSMIC_PASS:?COSMIC_PASS environment variable required}"

VER="v100"        # <-- update to the latest version
BUILD="GRCh38"
AUTH=$(printf '%s:%s' "$COSMIC_EMAIL" "$COSMIC_PASS" | base64)

# Product files to download (match the scripted-download paths on the page)
FILES=(
  "Cosmic_CancerGeneCensus_Tsv_${VER}_${BUILD}.tar"
  "Cosmic_MutantCensus_Tsv_${VER}_${BUILD}.tar"
)

API="https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted"
for f in "${FILES[@]}"; do
  echo "[req] $f"
  # 1) Request a signed download URL (returns JSON: {"url": "..."})
  url=$(curl -s -H "Authorization: Basic ${AUTH}" \
        "${API}?path=${BUILD}/cosmic/${VER}/${f}&bucket=downloads" \
        | sed -E 's/.*"url"\s*:\s*"([^"]+)".*/\1/')
  if [[ -z "$url" || "$url" == *'<'* ]]; then
    echo "  ! Could not get URL. Check the exact path/version on the download page."
    continue
  fi
  # 2) Actual download
  curl -L "$url" -o "$OUT/$f"
  echo "  [ok] $OUT/$f"
done

echo "[done] COSMIC. Downloaded files:"; ls -lh "$OUT"
