#!/usr/bin/env bash
# Full pipeline for one sample:
#   markdup (chr1-22) -> genotype likelihoods + allele depths at panel sites -> GLIMPSE2 impute/phase per chromosome
#   -> haplotype read counts -> depth (ichorCNA correctedDepth) + phased BAF bins -> ASCAT segmentation -> joint fit
# Resumable: each finished step leaves work/done/<sample>.<step> (holding its runtime in s) and is skipped on rerun.
# On failure, work/done/<sample>.FAILED names the step; delete it (and fix the cause) to retry.
# Usage: run_sample.sh <sample> <bam>
#   env: CORES (4)  REF  PY (python3)  KEEP_BAM (0: delete the markdup BAM at the end)  HAPLOID_RES
#        MODES ("arms cbs_union ascat": joint-fit segment modes -> work/joint_modes/<mode>/)
#        ICHOR_PARAMS / ICHOR_SEG / ICHOR_DEPTH (default: <bam dir>/<sample>_Y.params.txt / _Y.seg.txt / _Y.correctedDepth.txt)   GLIMPSE_* (glimpse_env.sh)
set -Eeuo pipefail      # -E: the ERR trap also fires inside step functions
S=$1; BAM_IN=$2
ROOT=$(cd "$(dirname "$0")/.." && pwd); SC=$ROOT/scripts
CORES=${CORES:-4}; PY=${PY:-python3}; KEEP_BAM=${KEEP_BAM:-0}
export REF=${REF:-$ROOT/resources/hg38.fa} HAPLOID_RES=${HAPLOID_RES:-$ROOT/resources}
BD=$(dirname "$BAM_IN")
ICHOR_PARAMS=${ICHOR_PARAMS:-$BD/${S}_Y.params.txt}; ICHOR_SEG=${ICHOR_SEG:-$BD/${S}_Y.seg.txt}
ICHOR_DEPTH=${ICHOR_DEPTH:-$BD/${S}_Y.correctedDepth.txt}
export ICHOR_PARAMS=${ICHOR_PARAMS//\{s\}/$S} ICHOR_SEG=${ICHOR_SEG//\{s\}/$S} ICHOR_DEPTH=${ICHOR_DEPTH//\{s\}/$S}  # "{s}" templates allowed
DONE=$ROOT/work/done; LOGD=$ROOT/logs/samples; MD=$ROOT/work/bam/$S.auto.md.bam
mkdir -p "$DONE" "$LOGD" "$ROOT/work/gl" "$ROOT/work/bam"
rm -f "$DONE/$S.FAILED"
CUR=start
trap 'echo "$CUR" > "$DONE/$S.FAILED"; echo "[$(date "+%F %T")] $S FAILED at step $CUR" >&2' ERR

step() {   # step <name> <function>: run once, record runtime
  CUR=$1
  [ -e "$DONE/$S.$1" ] && return 0
  echo "[$(date '+%F %T')] $S $1"
  local t=$SECONDS
  "$2"
  echo $((SECONDS - t)) > "$DONE/$S.$1"
}

do_markdup() {    # input BAMs are bwa -> sort with no duplicate marking
  BAM_IN=$BAM_IN THREADS=$CORES TAG=auto "$SC/markdup_pilot.sh" "$S" $(printf 'chr%s ' $(seq 1 22))
}
do_gl() {
  seq 1 22 | xargs -P "$CORES" -I{} sh -c '"$0" "$1" {} "$2" 2>> "$3"' \
    "$SC/gl_at_panel_sites.sh" "$MD" "$ROOT/work/gl/$S" "$LOGD/$S.gl.log"
}
do_phase() {      # chromosomes in turn, CORES single-threaded chunks at a time
  for c in $(seq 1 22); do
    [ -s "$ROOT/work/glimpse/$S/$S.chr$c.phased.bcf.csi" ] || "$SC/glimpse_phase_sample.sh" "$S" "$c" "$CORES"
  done
}
do_hap() {
  cd "$SC" && seq 1 22 | xargs -P "$CORES" -I{} "$PY" haplotype_counts.py "$S" {}
}
do_depth() {      # depth = ichorCNA correctedDepth; our own fragment counts only when that file is missing
  cd "$SC"
  if [ ! -s "$ICHOR_DEPTH" ]; then
    "$PY" fragment_counts.py "$S" > "$LOGD/$S.frag.log" && "$PY" depth_correction.py "$S"
  fi
  "$PY" build_bins.py "$S"
}
do_ascat() {
  local sex
  sex=$(awk -F'\t' '/^Gender:/ {print ($2 == "male") ? "XY" : "XX"}' "$ICHOR_PARAMS" 2>/dev/null || true)
  cd "$SC" && Rscript run_ascat.R "$S" "${sex:-XX}" logR 70 > "$LOGD/$S.ascat.log" 2>&1
}
do_joint() {      # one fit per segment mode: arms = detection + TF, cbs_union = events, ascat = reference
  cd "$SC"
  local m pids=() fail=0
  for m in ${MODES:-arms cbs_union ascat}; do     # modes are independent and single-threaded: run them side by side
    mkdir -p "$ROOT/work/joint_modes/$m"
    HAPLOID_SEG_MODE=$m HAPLOID_JOINT_OUT=$ROOT/work/joint_modes/$m "$PY" joint_fit.py "$S" > "$LOGD/$S.joint.$m.log" 2>&1 &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || fail=1; done
  for m in ${MODES:-arms cbs_union ascat}; do grep "===" "$LOGD/$S.joint.$m.log" | sed "s/^/  [$m] /" || true; done
  [ "$fail" = 0 ] || { echo "  joint fit failed (see $LOGD/$S.joint.<mode>.log)" >&2; return 1; }
}
do_cleanup() {
  rm -rf "$ROOT/work/glimpse/$S/chunks"
  [ "$KEEP_BAM" = 1 ] || rm -f "$MD" "$MD.bai"
}

echo "[$(date '+%F %T')] $S start ($BAM_IN, $CORES cores, GLIMPSE runner $(. "$SC/glimpse_env.sh"; echo "$GLIMPSE_RUNNER"))"
step markdup do_markdup
step gl      do_gl
step phase   do_phase
step hap     do_hap
step depth   do_depth
step ascat   do_ascat
step joint   do_joint
step cleanup do_cleanup
echo "[$(date '+%F %T')] $S DONE"
