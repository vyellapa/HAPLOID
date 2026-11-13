#!/usr/bin/env bash
# Assemble the build context for the full image (container/Dockerfile.full) from prepared resources.
# Allow-list: only public reference files go in. Never added: control_bin_mask / pipeline_bias (built from healthy
# controls), cohort_bias, fit_overrides / user constraints, panels of normals, the 1kGP panel BCFs (only needed to
# prepare the GLIMPSE2 reference, which is included ready-made).
# Usage: bundle_resources.sh <resources dir> <prepared GLIMPSE2 reference dir> <hg38.fa> <out dir>
# Files are hard-linked (no extra disk space; <out dir> must be on the same filesystem) and must be world-readable:
# Dockerfile.full COPYs them as they are (a chmod -R layer would store every file twice).
set -euo pipefail
RES=$1; GREF=$2; FA=$3; OUT=$4
die() { echo "ERROR: $*" >&2; exit 1; }
for c in $(seq 1 22); do
  [ -e "$GREF/.done.chr$c" ] || die "$GREF is not fully prepared (no .done.chr$c)"
  for f in sites_hg38/chr$c.sites.vcf.gz sites_hg38/chr$c.sites.tsv.gz maps_b38/chr$c.b38.gmap.gz; do
    [ -s "$RES/$f" ] || die "$RES/$f missing"
  done
done
for f in hg38-blacklist.v2.bed cytoband_b38.bed exclusion/sites_hg38.bed.gz exclusion/bins_500kb.tsv; do
  [ -s "$RES/$f" ] || die "$RES/$f missing"
done
grep -q dragen "$RES/exclusion/SOURCES.txt" && die "$RES/exclusion was built with the DRAGEN file; rebuild it (build_exclusion.py)"
[ -s "$FA" ] && [ -s "$FA.fai" ] || die "$FA (+ .fai) missing"

rm -rf "$OUT"; mkdir -p "$OUT/small/exclusion" "$OUT/hg38"
link() { mkdir -p "$(dirname "$2")"; ln "$1" "$2"; }
link "$FA" "$OUT/hg38/hg38.fa"; link "$FA.fai" "$OUT/hg38/hg38.fa.fai"
for f in hg38-blacklist.v2.bed cytoband_b38.bed; do cp "$RES/$f" "$OUT/small/$f"; chmod 644 "$OUT/small/$f"; done
for f in "$RES"/exclusion/*; do link "$f" "$OUT/small/exclusion/$(basename "$f")"; done
for c in $(seq 1 22); do
  for f in "$RES"/sites_hg38/chr$c.sites.tsv.gz*; do link "$f" "$OUT/small/sites_hg38/$(basename "$f")"; done
  # the sites VCF is rewritten without its ##bcftools_* history lines, which record local paths of the build machine
  v=$OUT/small/sites_hg38/chr$c.sites.vcf.gz
  gzip -dc "$RES/sites_hg38/chr$c.sites.vcf.gz" | grep -v -E '^##(bcftools|tabix)_' | bgzip -c > "$v"
  chmod 644 "$v"; bcftools index -f "$v"
  link "$RES/maps_b38/chr$c.b38.gmap.gz" "$OUT/small/maps_b38/chr$c.b38.gmap.gz"
  case $c in [1-4]) part=a ;; [5-9]) part=b ;; 1[0-5]) part=c ;; *) part=d ;; esac   # four ~3 GB image layers
  for f in "$GREF"/1kGP.chr${c}_*.bin "$GREF/chunks.chr$c.txt" "$GREF/.done.chr$c"; do
    link "$f" "$OUT/glimpse_ref_$part/$(basename "$f")"
  done
done
du -sh "$OUT"/* | sed 's/^/  /'
echo "build context ready: $OUT (no control masks, cohort bias, constraints or panels of normals)"
