#!/usr/bin/env bash
# Download the GLIMPSE2 v2.0.0 static Linux x86-64 binaries (same release as the docker image used in testing)
# into resources/glimpse_bin. Run once on a node with internet access (e.g. the head node).
# Usage: fetch_glimpse.sh
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); D=$ROOT/resources/glimpse_bin; mkdir -p "$D"
URL=https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.0
for t in chunk split_reference phase ligate; do
  f=$D/GLIMPSE2_${t}_static
  [ -x "$f" ] || { curl -fsSL -o "$f" "$URL/GLIMPSE2_${t}_static"; chmod +x "$f"; }
done
"$D/GLIMPSE2_phase_static" --help 2>&1 | head -5
