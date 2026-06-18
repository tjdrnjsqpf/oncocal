#!/usr/bin/env bash
# Initial server setup: directory structure + conda environment
# Usage: from the repo root, run  bash scripts/server/01_setup_env.sh
set -euo pipefail

# repo root = two levels above this script
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
echo "[repo root] $ROOT"

# 1) Create data/results directories
mkdir -p data/external/{dbnsfp,civic,cosmic,oncokb,cancerhotspots,clinvar}
mkdir -p data/{raw,interim,processed} results/{figures,tables} logs
echo "[ok] directory tree ready"

# 2) conda environment (requires environment.yml — assumed shipped with the repo)
# ⚠️ This server's default conda (<system anaconda>) is conda 4.8.3 / Python 2.7, too old
#    to solve python 3.11 → must use an up-to-date miniforge mamba.
MAMBA="${MAMBA:-/var2/lsg/miniforge3/bin/mamba}"
if [[ -x "$MAMBA" ]]; then
  if "$MAMBA" env list | grep -q '/mutfunc$'; then
    echo "[skip] conda env 'mutfunc' already exists"
  else
    "$MAMBA" env create -f environment.yml -y
    "$("$MAMBA" env list | awk '/\/mutfunc$/{print $NF}')/bin/pip" install -q "xlrd>=2.0.1" || true
    echo "[ok] env 'mutfunc' created. Usage: /var2/lsg/miniforge3/envs/mutfunc/bin/python"
  fi
else
  echo "[warn] miniforge mamba not found($MAMBA). Set the MAMBA env var to its path or install miniforge."
  echo "       (Note: anaconda2's conda/mamba is py2.7 and will fail.)"
fi

echo "[done] setup complete."
