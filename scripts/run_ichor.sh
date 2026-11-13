#!/usr/bin/env bash
# ichorCNA for one BAM at 500 kb, with the settings the pipeline was developed on: HMMcopy readCounter (MAPQ >= 20,
# chr1-22, X, Y) -> runIchorCNA.R with ichorCNA's own hg38 GC / mappability / centromere files; no panel of normals
# unless ICHOR_PON names one (build_controls.sh: <controls>/ichor_pon_median.rds, from the lab's own healthy controls).
# Usage: run_ichor.sh <sample> <bam> <outdir>
#   env: ICHOR_PON (panel of normals .rds), ICHOR_SCRIPT (runIchorCNA.R; default: found on PATH), READCOUNTER (readCounter)
# Output: <outdir>/<sample>_Y.correctedDepth.txt, _Y.params.txt, _Y.seg.txt (+ ichorCNA's plots in <outdir>/<sample>_Y/)
set -euo pipefail
S=$1; BAM=$2; OUT=$3
. "$(dirname "$0")/bam_lib.sh"
ICHOR_SCRIPT=${ICHOR_SCRIPT:-$(command -v runIchorCNA.R)}
mkdir -p "$OUT"
EXT=$(ichor_extdata)
PON=()
if [ -n "${ICHOR_PON:-}" ]; then
  [ -s "$ICHOR_PON" ] || { echo "panel of normals not found: $ICHOR_PON" >&2; exit 1; }
  PON=(--normalPanel "$ICHOR_PON")
fi

count_reads "$BAM" "$OUT/$S.wig"

Rscript "$ICHOR_SCRIPT" --id "${S}_Y" --WIG "$OUT/$S.wig" ${PON[@]+"${PON[@]}"} --ploidy "c(2,3)" \
  --normal "c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)" --maxCN 5 \
  --gcWig "$EXT/gc_hg38_500kb.wig" --mapWig "$EXT/map_hg38_500kb.wig" \
  --centromere "$EXT/GRCh38.GCA_000001405.2_centromere_acen.txt" \
  --includeHOMD False --chrs 'c(1:22, "X", "Y")' --chrTrain 'c(1:22)' \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates 'c(1,3)' --txnE 0.99999 --txnStrength 100000 --outDir "$OUT"
rm -f "$OUT"/*.RData
for f in correctedDepth.txt params.txt seg.txt; do
  [ -s "$OUT/${S}_Y.$f" ] || { echo "ichorCNA did not write $OUT/${S}_Y.$f" >&2; exit 1; }
done
echo "$S depth HMM${ICHOR_PON:+ (panel of normals)}: TF $(awk -F'\t' '/^Tumor Fraction:/ {print $2}' "$OUT/${S}_Y.params.txt"), ploidy $(awk -F'\t' '/^Ploidy:/ {print $2}' "$OUT/${S}_Y.params.txt")"
