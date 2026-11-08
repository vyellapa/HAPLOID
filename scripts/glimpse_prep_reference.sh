#!/usr/bin/env bash
# One-time GLIMPSE2 reference prep for one chromosome: chunk the lifted hg38 panel, then convert each
# chunk to GLIMPSE2's binary reference format.
# Usage: glimpse_prep_reference.sh <chrom number> [threads]
set -euo pipefail
CHR=$1; THREADS=${2:-8}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
. "$ROOT/scripts/glimpse_env.sh"
PANEL=$ROOT/resources/panel_hg38/1kGP_phase3.chr$CHR.hg38.biSNP_maf0.01.bcf
MAP=$ROOT/resources/maps_b38/chr$CHR.b38.gmap.gz
OUT=$ROOT/work/glimpse/ref; mkdir -p "$OUT"

glimpse GLIMPSE2_chunk --input "$PANEL" --region "chr$CHR" --map "$MAP" --sequential \
  --threads "$THREADS" --output "$OUT/chunks.chr$CHR.txt" --log "$OUT/chunks.chr$CHR.log" > /dev/null

while IFS=$'\t' read -r IDX _ IRG ORG _; do
  glimpse GLIMPSE2_split_reference --reference "$PANEL" --map "$MAP" --input-region "$IRG" --output-region "$ORG" \
    --threads "$THREADS" --output "$OUT/1kGP.chr$CHR" --log "$OUT/split.chr$CHR.$IDX.log" > /dev/null
done < "$OUT/chunks.chr$CHR.txt"
echo "chr$CHR: $(wc -l < "$OUT/chunks.chr$CHR.txt") chunks, $(ls "$OUT"/1kGP.chr$CHR*.bin | wc -l) binary reference files"
