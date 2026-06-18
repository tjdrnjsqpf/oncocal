#!/usr/bin/env bash
# dbNSFP 50GB 이어받기 루프 (연결 끊겨도 size 도달까지 재개)
set -u
URL="https://dist.genos.us/academic/e55b09/dbNSFP5.3.1a_grch38.gz"
OUT="data/external/dbnsfp/dbNSFP5.3.1a_grch38.gz"
EXP=50271494764
for i in $(seq 1 60); do
  sz=$(stat -c%s "$OUT" 2>/dev/null || echo 0)
  if [ "$sz" -ge "$EXP" ]; then echo "[$i] reached $sz >= $EXP"; break; fi
  echo "[$i] $(date +%H:%M:%S) resume from $sz ..."
  curl --http1.1 -C - --retry 5 --retry-delay 10 -s -o "$OUT" "$URL"
  sleep 8
done
sz=$(stat -c%s "$OUT" 2>/dev/null || echo 0)
[ "$sz" -ge "$EXP" ] && echo "DBNSFP_DL_DONE" || echo "DBNSFP_DL_INCOMPLETE($sz)"
