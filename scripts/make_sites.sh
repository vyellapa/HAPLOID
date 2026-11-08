#!/usr/bin/env bash
# Site tables for one chromosome from the lifted panel (liftover_panel.sh): the genotype likelihood step calls exactly
# these sites, and the joint fit reads their population allele frequencies.
#   resources/sites_hg38/chr<N>.sites.vcf.gz   panel sites without genotypes, INFO = AF and the 1kGP super-population AFs
#   resources/sites_hg38/chr<N>.sites.tsv.gz   CHROM POS REF,ALT (tabix-indexed; bcftools call -C alleles -T)
# Usage: make_sites.sh <chrom number>
set -euo pipefail
CHR=$1
ROOT=$(cd "$(dirname "$0")/.." && pwd)
P=$ROOT/resources/panel_hg38/1kGP_phase3.chr$CHR.hg38.biSNP_maf0.01.bcf
O=$ROOT/resources/sites_hg38; mkdir -p "$O"
bcftools view -G "$P" -Ou | bcftools annotate -x '^INFO/AF,INFO/EAS_AF,INFO/EUR_AF,INFO/AFR_AF,INFO/AMR_AF,INFO/SAS_AF' \
  -Oz -o "$O/chr$CHR.sites.vcf.gz"
bcftools index -f "$O/chr$CHR.sites.vcf.gz"
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' "$O/chr$CHR.sites.vcf.gz" | bgzip -c > "$O/chr$CHR.sites.tsv.gz"
tabix -f -s1 -b2 -e2 "$O/chr$CHR.sites.tsv.gz"
