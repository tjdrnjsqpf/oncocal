#!/usr/bin/env bash
# Download public data that requires no registration (run on the server)
# Usage: from the repo root, run  bash scripts/server/02_download_open_data.sh
# Note: some URLs change by version → check/update the latest links in the variables below.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
EXT="$ROOT/data/external"
log(){ echo "[$(date +%H:%M:%S)] $*"; }

# ---------------------------------------------------------------------------
# 1) ClinVar (NCBI, public domain) — stable direct link
# ---------------------------------------------------------------------------
log "ClinVar variant_summary..."
wget -c -P "$EXT/clinvar" \
  "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

# ---------------------------------------------------------------------------
# 2) CIViC (CC0) — nightly TSV (stable path). For a monthly release, swap in the releases page URL.
# ---------------------------------------------------------------------------
log "CIViC nightly TSVs..."
CIVIC_BASE="https://civicdb.org/downloads/nightly"
for f in VariantSummaries MolecularProfileSummaries ClinicalEvidenceSummaries \
         AssertionSummaries FeatureSummaries VariantGroupSummaries; do
  wget -c -P "$EXT/civic" "${CIVIC_BASE}/nightly-${f}.tsv" || \
    log "  ! ${f} failed — check directly at https://civicdb.org/releases"
done

# ---------------------------------------------------------------------------
# 3) cancerhotspots.org (no registration) — if the direct link is blocked, copy the link from the page and replace it
# ---------------------------------------------------------------------------
log "cancerhotspots v2..."
wget -c -P "$EXT/cancerhotspots" \
  "https://www.cancerhotspots.org/files/hotspots_v2.xls" || \
  log "  ! failed — check the link at https://www.cancerhotspots.org/#/download"

# ---------------------------------------------------------------------------
# 4) dbNSFP (academic, free) — ⚠️ the direct link differs by page/version, so set it via a variable
#    Copy the academic ('a') zip link from https://www.dbnsfp.org/download and put it below.
# ---------------------------------------------------------------------------
DBNSFP_URL=""   # <-- e.g.: "https://.../dbNSFP4.9a.zip"  (academic 'a' version!)
if [[ -n "$DBNSFP_URL" ]]; then
  log "dbNSFP download (large ~30GB+, takes a long time)..."
  wget -c -P "$EXT/dbnsfp" "$DBNSFP_URL"
else
  log "dbNSFP: DBNSFP_URL is empty → fill the academic zip link in the variable at the top of the script and rerun."
fi

log "open-data complete. Downloaded files:"
find "$EXT" -type f | sort
