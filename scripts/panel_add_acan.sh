#!/usr/bin/env bash
# Add AC/AN (required by GLIMPSE2) to an already-lifted panel file, in place.
# Usage: panel_add_acan.sh <chrom number>
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
P=$ROOT/resources/panel_hg38/1kGP_phase3.chr$1.hg38.biSNP_maf0.01.bcf
bcftools +fill-tags "$P" -Ob -o "$P.tmp" -- -t AC,AN && mv "$P.tmp" "$P" && bcftools index -f "$P"
