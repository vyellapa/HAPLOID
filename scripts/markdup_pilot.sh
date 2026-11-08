#!/usr/bin/env bash
# Mark duplicates on the pilot chromosomes only (input BAMs are bwa -> sort, no markdup).
# Usage: markdup_pilot.sh <sample> [chroms...]   env: BAM_IN (default <project>/<sample>.bam), THREADS (4), TAG
set -euo pipefail
S=$1; shift; CHROMS=${*:-chr16 chr17 chr22}; TAG=${TAG:-pilot}
ROOT=$(cd "$(dirname "$0")/.." && pwd); OUT=$ROOT/work/bam; T=$OUT/tmp_${S}_$TAG; mkdir -p "$T"
TH=${THREADS:-4}
samtools view -@"$TH" -u "${BAM_IN:-$ROOT/$S.bam}" $CHROMS \
  | samtools collate -@"$TH" -O -u - "$T/col" \
  | samtools fixmate -@"$TH" -m -u - - \
  | samtools sort -@"$TH" -m 2G -T "$T/srt" -u - \
  | samtools markdup -@"$TH" -s -f "$OUT/$S.$TAG.markdup_stats.txt" - "$OUT/$S.$TAG.md.bam"
samtools index "$OUT/$S.$TAG.md.bam"; rm -rf "$T"
