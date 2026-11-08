#!/usr/bin/env bash
# Genotype likelihoods (GLIMPSE2 input) + allele depths at panel sites for one sample/chromosome.
# Duplicates (0x400) are skipped by mpileup's default filter; overlapping mates are
# counted once (mpileup overlap detection stays on).
# Usage: gl_at_panel_sites.sh <bam> <chrom number> <out prefix>
set -euo pipefail
BAM=$1; CHR=$2; OUT=$3
ROOT=$(cd "$(dirname "$0")/.." && pwd)
REF=${REF:-$ROOT/resources/hg38.fa}
SITES_VCF=$ROOT/resources/sites_hg38/chr$CHR.sites.vcf.gz
SITES_TSV=$ROOT/resources/sites_hg38/chr$CHR.sites.tsv.gz

bcftools mpileup -f "$REF" -I -E -q 20 -Q 20 -a 'FORMAT/DP,FORMAT/AD' \
    -T "$SITES_VCF" -r "chr$CHR" -Ou "$BAM" \
  | bcftools call -Aim -C alleles -T "$SITES_TSV" -Ob -o "$OUT.chr$CHR.gl.bcf"
bcftools index -f "$OUT.chr$CHR.gl.bcf"
