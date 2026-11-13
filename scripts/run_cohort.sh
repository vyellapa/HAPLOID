#!/usr/bin/env bash
# Run the HAPLOID pipeline (scripts/run_sample.sh) on a cohort, JOBS samples at a time with CORES cores each.
#   1. preflight: tools, python/R packages, reference, project resources, GLIMPSE2 runner
#   2. one-time shared prep: panel AC/AN tags, GLIMPSE2 binary reference (22 chromosomes), GC 'possible' table
#   3. samples in parallel; each is resumable, so rerunning this script picks up where it stopped
#   4. cohort table (cohort_table.py: work/cohort_summary.tsv, work/cohort_events.tsv), a line per failed sample, and
#      the second pass: cohort_bias.py (residual CSF-specific depth bias) -> refit_downstream.sh -> cohort_table.py
# Usage: run_cohort.sh [bam ...]          (no arguments: every $BAM_DIR/*_500k/*.bam)
# env: BAM_DIR, REF, JOBS (15), CORES (4), PREP_JOBS (16), PY (python3), KEEP_BAM (0), GLIMPSE_RUNNER / GLIMPSE_BIN,
#      MODES (joint-fit segment modes, "arms cbs_union ascat"), ICHOR_* (run_sample.sh), ALLOW_FRAGMENT_DEPTH (0),
#      HAPLOID_GLIMPSE_REF (prepared GLIMPSE2 reference; default work/glimpse/ref), HAPLOID_CONTROLS (control mask dir)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); SC=$ROOT/scripts
BAM_DIR=${BAM_DIR:?set BAM_DIR (folder of <sample>_500k/<sample>.bam)}
export REF=${REF:-$ROOT/resources/hg38.fa}
export CORES=${CORES:-4} PY=${PY:-python3} KEEP_BAM=${KEEP_BAM:-0} HAPLOID_RES=${HAPLOID_RES:-$ROOT/resources}
JOBS=${JOBS:-15}; PREP_JOBS=${PREP_JOBS:-16}
GREF=${HAPLOID_GLIMPSE_REF:-$ROOT/work/glimpse/ref}   # prepared GLIMPSE2 reference (built in the container)
CTRL=${HAPLOID_CONTROLS:-$ROOT/resources}
. "$SC/glimpse_env.sh"
LOG=$ROOT/logs; DONE=$ROOT/work/done; mkdir -p "$LOG/samples" "$DONE"
ts() { date '+%F %T'; }
die() { echo "ERROR: $*" >&2; exit 1; }

# ---- 1. preflight -------------------------------------------------------------------------------------------
echo "[$(ts)] preflight"
for t in samtools bcftools Rscript "$PY"; do command -v "$t" > /dev/null || die "$t not in PATH"; done
"$PY" -c 'import numpy, pandas, pysam, matplotlib' 2> /dev/null || die "$PY needs numpy pandas pysam matplotlib"
Rscript -e 'suppressMessages({library(ASCAT); library(DNAcopy)})' > /dev/null 2>&1 || die "R needs ASCAT and DNAcopy"
[ -s "$REF" ] && [ -s "$REF.fai" ] || die "reference $REF (+ .fai) not found"
for f in hg38-blacklist.v2.bed cytoband_b38.bed; do [ -s "$HAPLOID_RES/$f" ] || die "$HAPLOID_RES/$f missing"; done
# region exclusions: without these files the pipeline silently runs unmasked, so they are required (build_exclusion.py)
for f in resources/exclusion/sites_hg38.bed.gz resources/exclusion/bins_500kb.tsv; do
  [ -s "$ROOT/$f" ] || die "$f missing (see resources/README.md)"
done
if [ -s "$CTRL/control_bin_mask_500kb.tsv" ]; then echo "  bin mask: healthy controls ($CTRL/control_bin_mask_500kb.tsv)"
elif [ "${HAPLOID_PUBLIC_MASK:-1}" != 0 ] && [ -s "$ROOT/resources/exclusion/public_bin_mask_500kb.tsv" ]; then
  echo "  bin mask: public (hg38 alt-haplotype loci + top-1% GC bins; resources/exclusion/public_bin_mask_500kb.tsv)"
else echo "  WARNING: no bin mask - bias-prone bins (MHC and other alt-haplotype regions) enter the depth fit"; fi
[ -s "$ROOT/resources/cohort_bias_500kb.tsv" ] && echo "  cohort bias present: second pass (resources/cohort_bias_500kb.tsv applied)"
[ -s "$ROOT/work/depth/ichor_pon_hg38_500kb.tsv" ] \
  || echo "  note: work/depth/ichor_pon_hg38_500kb.tsv missing - only needed when a sample has no ichorCNA correctedDepth"
ref_todo=$(for c in $(seq 1 22); do [ -e "$GREF/.done.chr$c" ] || echo "$c"; done)
for c in $(seq 1 22); do
  for f in sites_hg38/chr$c.sites.vcf.gz sites_hg38/chr$c.sites.tsv.gz maps_b38/chr$c.b38.gmap.gz \
           $(grep -qx "$c" <<< "$ref_todo" && echo panel_hg38/1kGP_phase3.chr$c.hg38.biSNP_maf0.01.bcf); do
    [ -s "$ROOT/resources/$f" ] || die "resources/$f missing"     # the panel only while its GLIMPSE2 reference is unprepared
  done
done
ghelp=$(glimpse GLIMPSE2_phase --help 2>&1 || true)
grep -q -- --bam-file <<< "$ghelp" \
  || die "GLIMPSE2 not runnable with GLIMPSE_RUNNER=$GLIMPSE_RUNNER (run scripts/fetch_glimpse.sh on a node with internet)"
echo "  GLIMPSE2 runner: $GLIMPSE_RUNNER;  free space under work/: $(df -h "$ROOT/work" | awk 'NR==2 {print $4}')"

# ---- samples and their ichorCNA depth -----------------------------------------------------------------------
if [ $# -gt 0 ]; then BAMS=("$@"); else BAMS=("$BAM_DIR"/*_500k/*.bam); fi
[ -e "${BAMS[0]}" ] || die "no BAMs found"
n_noichor=0
for b in "${BAMS[@]}"; do
  s=$(basename "$b" .bam); f=${ICHOR_DEPTH:-$(dirname "$b")/${s}_Y.correctedDepth.txt}; f=${f//\{s\}/$s}
  [ -s "$f" ] || { n_noichor=$((n_noichor + 1)); echo "  no ichorCNA correctedDepth for $s ($f)"; }
done
[ "$n_noichor" -eq 0 ] || [ "${ALLOW_FRAGMENT_DEPTH:-0}" = 1 ] \
  || die "$n_noichor samples lack ichorCNA correctedDepth (their depth would come from a different method); fix the paths or set ALLOW_FRAGMENT_DEPTH=1"

# ---- 2. one-time shared prep --------------------------------------------------------------------------------
for c in $ref_todo; do
  p=$ROOT/resources/panel_hg38/1kGP_phase3.chr$c.hg38.biSNP_maf0.01.bcf
  bcftools view -h "$p" | grep '^##INFO=<ID=AC,' > /dev/null || echo "$c"
done | xargs -r -P 8 -I{} "$SC/panel_add_acan.sh" {}

todo=$ref_todo
if [ -n "$todo" ]; then
  echo "[$(ts)] GLIMPSE2 reference prep for chromosomes: $(echo $todo)"
  echo "$todo" | xargs -P "$PREP_JOBS" -I{} sh -c \
    '"$0" {} 4 > "$1/glimpse_prep.chr{}.log" 2>&1 && touch "$2/.done.chr{}"' \
    "$SC/glimpse_prep_reference.sh" "$LOG" "$GREF" \
    || die "reference prep failed (see $LOG/glimpse_prep.chr*.log)"
fi

if [ "$n_noichor" -gt 0 ] && [ ! -s "$ROOT/work/depth/possible.npz" ]; then   # only for the fragment-count fallback
  echo "[$(ts)] GC 'possible fragments' table"
  (cd "$SC" && "$PY" fragment_counts.py --possible > "$LOG/frag_possible.log" 2>&1) || die "possible table failed"
fi

# ---- 3. samples ---------------------------------------------------------------------------------------------
n_done=0; LIST=$LOG/cohort_todo.txt; : > "$LIST"
for b in "${BAMS[@]}"; do
  s=$(basename "$b" .bam)
  if [ -e "$DONE/$s.cleanup" ]; then n_done=$((n_done + 1)); else echo "$b" >> "$LIST"; fi
done
echo "[$(ts)] ${#BAMS[@]} samples: $n_done already done, $(wc -l < "$LIST") to run; $JOBS at a time x $CORES cores"
echo "  progress:  ls $DONE/*.cleanup | wc -l      failures: ls $DONE/*.FAILED      per-sample logs: $LOG/samples/"

xargs -P "$JOBS" -I{} bash -c 'b={}; s=$(basename "$b" .bam); "$0" "$s" "$b" >> "$1/$s.log" 2>&1 || true' \
  "$SC/run_sample.sh" "$LOG/samples" < "$LIST"

# ---- 4. summary ---------------------------------------------------------------------------------------------
n_ok=0; n_fail=0
for b in "${BAMS[@]}"; do
  s=$(basename "$b" .bam)
  if [ -e "$DONE/$s.cleanup" ]; then n_ok=$((n_ok + 1))
  elif [ -e "$DONE/$s.FAILED" ]; then n_fail=$((n_fail + 1)); echo "  FAILED $s at step $(cat "$DONE/$s.FAILED") (log: $LOG/samples/$s.log)"
  fi
done
echo "[$(ts)] finished: $n_ok done, $n_fail failed of ${#BAMS[@]}"
(cd "$SC" && "$PY" clonal_tf.py > /dev/null && "$PY" cohort_table.py) || echo "  clonal_tf.py / cohort_table.py failed"
if [ ! -s "$ROOT/resources/cohort_bias_500kb.tsv" ] && [ -z "${HAPLOID_NO_NEXT:-}" ]; then cat << NEXT
Next (second pass, ~minutes per sample, no re-phasing):
  1. $PY $SC/cohort_bias.py        # residual bias from the tumor-negative samples; REVIEW the printed blocks + work/cohort_bias.png
  2. ls $DONE/*.cleanup | xargs -n1 basename | sed 's/.cleanup\$//' > $LOG/cohort_samples.txt
     BAM_DIR=$BAM_DIR JOBS=$JOBS bash $SC/refit_downstream.sh \$(cat $LOG/cohort_samples.txt)
  3. $PY $SC/clonal_tf.py && $PY $SC/cohort_table.py       # work/cohort_summary.tsv + work/cohort_events.tsv
  4. BAM_DIR=$BAM_DIR JOBS=$JOBS bash $SC/resolve_discordant.sh   # cbs_union ploidy off ichorCNA's by > 0.5: refit with
     ploidy held at ichorCNA's, TF free; add sample names to leave out suspected germline-CNV samples
NEXT
fi
