#!/usr/bin/env bash
# 서버 초기 셋업: 디렉터리 구조 + conda 환경
# 사용법: repo 루트에서  bash scripts/server/01_setup_env.sh
set -euo pipefail

# repo 루트 = 이 스크립트의 두 단계 상위
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
echo "[repo root] $ROOT"

# 1) 데이터/결과 디렉터리 생성
mkdir -p data/external/{dbnsfp,civic,cosmic,oncokb,cancerhotspots,clinvar}
mkdir -p data/{raw,interim,processed} results/{figures,tables} logs
echo "[ok] directory tree ready"

# 2) conda 환경 (environment.yml 필요 — repo와 함께 전송됐다고 가정)
# ⚠️ 이 서버의 기본 conda(<system anaconda>)는 conda 4.8.3 / Python 2.7 로 너무 오래되어
#    python 3.11 을 solve 하지 못한다 → 반드시 최신 miniforge mamba 사용.
MAMBA="${MAMBA:-/var2/lsg/miniforge3/bin/mamba}"
if [[ -x "$MAMBA" ]]; then
  if "$MAMBA" env list | grep -q '/mutfunc$'; then
    echo "[skip] conda env 'mutfunc' already exists"
  else
    "$MAMBA" env create -f environment.yml -y
    "$("$MAMBA" env list | awk '/\/mutfunc$/{print $NF}')/bin/pip" install -q "xlrd>=2.0.1" || true
    echo "[ok] env 'mutfunc' created. 사용: /var2/lsg/miniforge3/envs/mutfunc/bin/python"
  fi
else
  echo "[warn] miniforge mamba 없음($MAMBA). MAMBA 환경변수로 경로 지정 또는 miniforge 설치."
  echo "       (주의: anaconda2의 conda/mamba 는 py2.7 이라 실패함.)"
fi

echo "[done] setup complete."
