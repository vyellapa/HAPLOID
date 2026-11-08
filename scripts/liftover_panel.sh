#!/usr/bin/env bash
# Lift one chromosome of the 1kGP Phase 3 phased panel (b37) to hg38.
#
#   1. keep biallelic SNPs with MAF >= MIN_MAF (rare sites carry no LOH signal at 2x)
#   2. lift positions with UCSC liftOver (hg19ToHg38 chain; b37 == hg19 on chr1-22)
#   3. drop sites that move to another chromosome / alt contig, or collide after lifting
#   4. drop sites whose REF no longer matches hg38.fa (ref changes, strand-flipped segments)
#
# Usage: liftover_panel.sh <chrom number, e.g. 22>
set -euo pipefail

CHR=$1
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PANEL_B37=${PANEL_B37:?set PANEL_B37 (folder of 1000 Genomes phase 3 GRCh37 VCFs)}
REF=${REF:-$ROOT/resources/hg38.fa}
CHAIN=$ROOT/resources/chain/hg19ToHg38.over.chain.gz
MIN_MAF=${MIN_MAF:-0.01}
OUT=$ROOT/resources/panel_hg38
TMP=$OUT/tmp_chr$CHR
mkdir -p "$TMP"

IN=$PANEL_B37/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
FINAL=$OUT/1kGP_phase3.chr$CHR.hg38.biSNP_maf${MIN_MAF}.bcf

# 1. filter; keep population AFs, tag each record with its b37 position
bcftools view -m2 -M2 -v snps -q "$MIN_MAF:minor" -Ou "$IN" \
  | bcftools annotate -x '^INFO/AC,INFO/AN,INFO/AF,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/EUR_AF,INFO/SAS_AF' \
      --set-id 'b37_%POS' -Ob -o "$TMP/b37.bcf"
N_IN=$(bcftools view -H "$TMP/b37.bcf" | wc -l | tr -d ' ')

# 2. lift positions
bcftools query -f 'chr%CHROM\t%POS0\t%POS\t%ID\t0\t+\n' "$TMP/b37.bcf" > "$TMP/b37.bed"
liftOver "$TMP/b37.bed" "$CHAIN" "$TMP/hg38.bed" "$TMP/unmapped.bed" 2> "$TMP/liftover.log"

# 3. same chromosome, single base, + strand, unique target position.
#    Minus-strand sites must go: palindromic SNPs (A/T, C/G) there still pass the REF
#    check in step 4 but have REF/ALT (and so every genotype) swapped.
awk -v c="chr$CHR" '$1 == c && $3 - $2 == 1 && $6 == "+"' "$TMP/hg38.bed" \
  | awk '{n[$3]++; row[NR]=$0; p[NR]=$3} END {for (i=1; i<=NR; i++) if (n[p[i]] == 1) print row[i]}' \
  > "$TMP/hg38.ok.bed"

# rewrite CHROM/POS; header gets the hg38 contig line
{
  bcftools view -h "$TMP/b37.bcf" | grep -v '^##contig' | sed '$d'
  awk -v c="chr$CHR" '$1 == c {printf "##contig=<ID=%s,length=%s,assembly=hg38>\n", $1, $2}' "$REF.fai"
  echo '##INFO_NOTE=<Description="Lifted from 1kGP phase3 b37 with UCSC liftOver hg19ToHg38; ID holds the b37 position">'
  bcftools view -h "$TMP/b37.bcf" | tail -1
  bcftools view -H "$TMP/b37.bcf" \
    | awk -v c="chr$CHR" 'BEGIN {FS = OFS = "\t"}
        NR == FNR {map[$4] = $3; next}
        ($3 in map) {$1 = c; $2 = map[$3]; print}' "$TMP/hg38.ok.bed" -
} | bcftools sort -m 4G -T "$TMP" -Ob -o "$TMP/hg38.unchecked.bcf"
N_LIFT=$(bcftools view -H "$TMP/hg38.unchecked.bcf" | wc -l | tr -d ' ')

# 4. REF must match hg38. `norm --check-ref x` segfaults in bcftools 1.22, so list
#    mismatches with `w` (on a sites-only stream) and exclude them explicitly.
bcftools view -G -Ou "$TMP/hg38.unchecked.bcf" \
  | bcftools norm -f "$REF" --check-ref w -Ou -o /dev/null 2>&1 \
  | awk 'BEGIN {OFS = "\t"} $1 == "REF_MISMATCH" {print $2, $3}' > "$TMP/ref_mismatch.tsv"
bcftools index -f "$TMP/hg38.unchecked.bcf"
if [ -s "$TMP/ref_mismatch.tsv" ]; then
  bcftools view -T "^$TMP/ref_mismatch.tsv" -Ob -o "$FINAL" "$TMP/hg38.unchecked.bcf"
else
  cp "$TMP/hg38.unchecked.bcf" "$FINAL"
fi
bcftools index -f "$FINAL"
N_OUT=$(bcftools index -n "$FINAL")

printf "chr%s\tb37_biSNP_maf%s=%s\tlifted_unique=%s\tref_ok=%s\tkept=%.2f%%\n" \
  "$CHR" "$MIN_MAF" "$N_IN" "$N_LIFT" "$N_OUT" "$(echo "100*$N_OUT/$N_IN" | bc -l)" \
  | tee "$OUT/chr$CHR.liftover_stats.tsv"
rm -rf "$TMP"
