#!/usr/bin/env bash
# In-silico dilution: mix unique (non-duplicate, primary, mapped) read pairs from a tumor BAM and a same-patient
# low-TF BAM (same germline, so allelic signals stay valid). Pairs are subsampled by read-name hash (mates kept
# together). Expected TF of the mixture = f * TF_tumor + (1 - f) * TF_diluent, f = tumor share of fragments.
# Inputs are duplicate-marked BAMs (work/bam/<s>.auto.md.bam); the output has no duplicates.
# Usage: make_mixture.sh <tumor.md.bam> <diluent.md.bam> <f> <total_pairs> <out.bam> [seed]
set -euo pipefail
T=$1; D=$2; F=$3; N=$4; OUT=$5; SEED=${6:-11}
TH=${THREADS:-4}
pairs() {   # unique read-1 count, cached next to the BAM
  [ -s "$1.uniq_pairs" ] || samtools view -@"$TH" -c -F 0xD04 -f 0x40 "$1" > "$1.uniq_pairs"
  cat "$1.uniq_pairs"
}
UT=$(pairs "$T"); UD=$(pairs "$D")
read -r FT FD <<< "$(awk -v f="$F" -v n="$N" -v ut="$UT" -v ud="$UD" 'BEGIN { printf "%.7f %.7f", f * n / ut, (1 - f) * n / ud }')"
awk -v a="$FT" -v b="$FD" 'BEGIN { exit !(a <= 1 && b <= 1) }' \
  || { echo "not enough unique pairs: need fractions $FT (tumor) and $FD (diluent) <= 1" >&2; exit 1; }
sub() {   # sub <bam> <fraction> <out>: all unique pairs when fraction ~1, else hash-subsample
  if awk -v x="$2" 'BEGIN { exit !(x > 0.9999999) }'; then
    samtools view -@"$TH" -u -F 0xD04 -o "$3" "$1"
  elif awk -v x="$2" 'BEGIN { exit !(x > 0) }'; then
    samtools view -@"$TH" -u -F 0xD04 -s "$SEED${2#0}" -o "$3" "$1"
  fi
}
TMP=$(mktemp -d "${OUT%.bam}.tmp.XXXX")
sub "$T" "$FT" "$TMP/t.bam"; sub "$D" "$FD" "$TMP/d.bam"
samtools merge -@"$TH" -f -o "$OUT" "$TMP"/*.bam
samtools index "$OUT"; rm -rf "$TMP"
printf "tumor_bam\tdiluent_bam\tf\ttotal_pairs\tfrac_tumor\tfrac_diluent\ttumor_uniq_pairs\tdiluent_uniq_pairs\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$T" "$D" "$F" "$N" "$FT" "$FD" "$UT" "$UD" > "${OUT%.bam}.design.tsv"
echo "$(basename "$OUT"): f=$F, $(samtools view -c -f 0x40 "$OUT") pairs (tumor frac $FT of $UT, diluent frac $FD of $UD)"
