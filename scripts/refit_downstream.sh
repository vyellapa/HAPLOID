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
mkdir -p "$ROOT/logs/samples"
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
pids=(); failed=0
wait_batch() {
  local pid
  for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
  pids=()
}
for s in "$@"; do                          # JOBS samples at a time; a failure is logged, the others continue
  (
    trap 'status=$?; if [ "$status" -ne 0 ]; then echo "[$(date '+%T')] $s refit FAILED (see logs/samples/$s.refit.log)" >&2; fi' EXIT
    one "$s"
  ) &
  pids+=($!)
  if [ "${#pids[@]}" -ge "${JOBS:-4}" ]; then wait_batch; fi
done
wait_batch
exit "$failed"
