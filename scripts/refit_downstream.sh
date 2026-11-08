#!/usr/bin/env bash
# Re-run everything after phasing (bins -> ASCAT/CBS segmentation -> joint fit in several segment modes) for samples
# whose phased haplotype counts already exist, e.g. after changing masks, exclusions or the call filter.
# ichorCNA files are looked up in $BAM_DIR/<s>_500k (cluster), work/dilution/<s>_500k, <s>_500k, then controls/<s>_500k;
# without a correctedDepth
# file the depth falls back to our fragment counts.
# Usage: refit_downstream.sh <sample> [...]      env: MODES ("ascat arms cbs_union"), JOBS (4), PY,
#        STEPS ("bins seg joint"; e.g. STEPS=joint to refit only), BAM_DIR, HAPLOID_OVERRIDES (joint_fit.py constraints),
#        OUT_SUFFIX (output to work/joint_modes/<mode><suffix>, e.g. "_constrained", keeping the free fits to compare)
# Run with bash (zsh does not split "$@" lists given as one string).  Whole cohort:  refit_downstream.sh $(cat list)
# Output: work/joint_modes/<mode>${OUT_SUFFIX}/
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); SC=$ROOT/scripts
export PY=${PY:-python3} MODES=${MODES:-arms cbs_union ascat} STEPS=${STEPS:-bins seg joint} ROOT SC
one() {
  local s=$1 d=""
  for c in ${BAM_DIR:+"$BAM_DIR/${s}_500k"} "$ROOT/work/dilution/${s}_500k" "$ROOT/${s}_500k" "$ROOT/controls/${s}_500k"; do [ -d "$c" ] && { d=$c; break; }; done
  export ICHOR_DEPTH=$d/${s}_Y.correctedDepth.txt ICHOR_PARAMS=$d/${s}_Y.params.txt ICHOR_SEG=$d/${s}_Y.seg.txt
  local sex; sex=$(awk -F'\t' '/^Gender:/ {print ($2 == "male") ? "XY" : "XX"}' "$ICHOR_PARAMS" 2>/dev/null || true)
  cd "$SC"
  : > "$ROOT/logs/samples/$s.refit.log"
  case " $STEPS " in *" bins "*) "$PY" build_bins.py "$s" >> "$ROOT/logs/samples/$s.refit.log" 2>&1 ;; esac
  case " $STEPS " in *" seg "*) Rscript run_ascat.R "$s" "${sex:-XX}" logR 70 >> "$ROOT/logs/samples/$s.refit.log" 2>&1 ;; esac
  case " $STEPS " in *" joint "*) ;; *) echo "[$(date '+%T')] $s refit done"; return 0 ;; esac
  for m in $MODES; do
    mkdir -p "$ROOT/work/joint_modes/$m${OUT_SUFFIX:-}"
    HAPLOID_SEG_MODE=$m HAPLOID_JOINT_OUT=$ROOT/work/joint_modes/$m${OUT_SUFFIX:-} "$PY" joint_fit.py "$s" >> "$ROOT/logs/samples/$s.refit.log" 2>&1
  done
  echo "[$(date '+%T')] $s refit done"
}
i=0
for s in "$@"; do                          # JOBS samples at a time; a failure is logged, the others continue
  ( one "$s" || echo "[$(date '+%T')] $s refit FAILED (see logs/samples/$s.refit.log)" ) &
  i=$((i + 1)); [ $((i % ${JOBS:-4})) -eq 0 ] && wait
done
wait
