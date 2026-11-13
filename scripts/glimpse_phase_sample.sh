#!/usr/bin/env bash
# Phase + impute one sample on one chromosome with GLIMPSE2 (BAM mode), then ligate the chunks.
# BAM mode computes genotype likelihoods from reads; flagged duplicates are skipped by GLIMPSE2's default.
# Single-sample GLIMPSE2_phase runs on ~1 core regardless of --threads (threads parallelise over target samples),
# so chunks are run concurrently instead: JOBS at a time, 1 thread each (runner: glimpse_env.sh).
# Usage: glimpse_phase_sample.sh <sample> <chrom number> [jobs]
set -euo pipefail
S=$1; CHR=$2; JOBS=${3:-4}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
. "$ROOT/scripts/glimpse_env.sh"
BAM=$ROOT/work/bam/$S.auto.md.bam
REF=${HAPLOID_GLIMPSE_REF:-$ROOT/work/glimpse/ref}
OUT=$ROOT/work/glimpse/$S; mkdir -p "$OUT/chunks"
export ROOT BAM REF OUT S CHR

phase_chunk() {   # args: chunk index, buffered input region
  local O=$OUT/chunks/chr$CHR.$1.bcf BIN=$REF/1kGP.chr${CHR}_${2//[:\-]/_}.bin   # bins named by input region
  [ -s "$O" ] && [ -s "$O.csi" ] && return 0                                       # resume: chunk already done
  glimpse GLIMPSE2_phase --bam-file "$BAM" \
    --ind-name "$S" --reference "$BIN" --threads 1 --output "$O" --log "$O.log" > /dev/null
  [ -s "$O" ] || { echo "GLIMPSE2_phase produced no output for $O (see $O.log)" >&2; return 1; }
}
export -f phase_chunk

cut -f1,3 "$REF/chunks.chr$CHR.txt" | xargs -P "$JOBS" -n 2 bash -c 'phase_chunk "$0" "$1"'

LIST=$OUT/chunks/chr$CHR.list
cut -f1 "$REF/chunks.chr$CHR.txt" | while read -r IDX; do echo "$OUT/chunks/chr$CHR.$IDX.bcf"; done > "$LIST"

# this GLIMPSE2 build's ligate has no --log and exits 0 on errors: capture output, check the file
glimpse GLIMPSE2_ligate --input "$LIST" \
  --output "$OUT/$S.chr$CHR.phased.bcf" --threads 2 > "$OUT/ligate.chr$CHR.log" 2>&1
[ -s "$OUT/$S.chr$CHR.phased.bcf" ] || { echo "GLIMPSE2_ligate failed (see $OUT/ligate.chr$CHR.log)" >&2; exit 1; }
[ -s "$OUT/$S.chr$CHR.phased.bcf.csi" ] || bcftools index -f "$OUT/$S.chr$CHR.phased.bcf"
echo "$S chr$CHR: $(bcftools index -n "$OUT/$S.chr$CHR.phased.bcf" 2>/dev/null) sites"
