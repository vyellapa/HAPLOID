#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SC=$ROOT/scripts
export REF=${REF:-$ROOT/resources/hg38.fa} HAPLOID_RES=${HAPLOID_RES:-$ROOT/resources}
export GLIMPSE_RUNNER=${GLIMPSE_RUNNER:-native} PY=${PY:-python3}
export BAM_DIR=${BAM_DIR:-/data/bams} CORES=${CORES:-4} JOBS=${JOBS:-1} PREP_JOBS=${PREP_JOBS:-1}
for d in work logs; do                    # full image: work/ and logs/ link into the /out mount
  if [ -L "$ROOT/$d" ] && [ ! -e "$ROOT/$d" ]; then mkdir -p "$(readlink "$ROOT/$d")" 2>/dev/null \
    || { echo "ERROR: mount a writable folder at /out (e.g. -v \"$PWD/results:/out\")" >&2; exit 1; }; fi
done
cd "$ROOT"

usage() {
  cat <<'HELP'
Usage: haploid COMMAND [arguments...]

  run [--controls DIR] --bam BAM [--sample NAME] [--ichor DIR] [--bam ...]
  run [--controls DIR] --sheet SAMPLES.tsv
                     Whole analysis: read depth + depth HMM (unless --ichor gives
                     existing 500 kb results), phasing, joint fit, ploidy check,
                     plots and a PDF. Sheet: tab-separated, header sample, bam[,
                     ichor]. Results in $HAPLOID_FINAL (full image: /out/final).
                     --controls: a folder made by 'controls' (its bin mask replaces
                     the public one; its panel of normals is used for depth).
  controls --bam CTRL.bam ... [--wig CTRL.wig ...] [--list FILE] [--out DIR]
                     Your own healthy-control references from healthy-control BAMs
                     aligned exactly like the samples (>= 10 recommended): depth
                     panel of normals + HAPLOID bin mask. Default output:
                     $HAPLOID_CONTROLS_OUT (full image: /out/controls).
  cohort [BAM ...]    Prepare shared references, run samples, and build reports.
                     Without arguments, use $BAM_DIR/*_500k/*.bam.
  sample BAM         Run one BAM through the same preparation and full pipeline.
                     Its filename, without .bam, is the sample name.
  report             Rebuild clonal estimates and cohort summary/event tables.
  bias [options...]  Build a cohort depth-bias table; review it before refitting.
  refit SAMPLE ...   Rebuild bins, segmentation, fits, and reports after phasing.
  constrain ARGS...  Refit with known values (refit_with.sh).
  resolve [SAMPLE...] Refit discordant ploidies; arguments exclude samples.
  plot SAMPLE ...    Write genome-wide comparison plots.
  check              Verify installed tools and Python/R packages (no data needed).
  exec COMMAND ...   Run any command, including reference-preparation scripts.

Mount prepared reference resources at /opt/haploid/resources, BAMs and optional
depth results at /data/bams, and persistent output directories at
/opt/haploid/work and /opt/haploid/logs. Set REF if the FASTA has another name.
Defaults: JOBS=1, CORES=4, PREP_JOBS=1. See README.md for run examples.
HELP
}

report() {
  "$PY" "$SC/clonal_tf.py"
  "$PY" "$SC/cohort_table.py"
}

cohort() {
  local bam sample failed=0
  local bams=("$@")
  if [ "${#bams[@]}" -eq 0 ]; then
    shopt -s nullglob
    bams=("$BAM_DIR"/*_500k/*.bam)
    shopt -u nullglob
  fi
  if [ "${#bams[@]}" -eq 0 ]; then
    echo "ERROR: no BAMs found under $BAM_DIR/*_500k/*.bam" >&2
    return 1
  fi
  for bam in "${bams[@]}"; do
    [ -s "$bam" ] || { echo "ERROR: BAM not found or empty: $bam" >&2; return 1; }
  done
  bash "$SC/run_cohort.sh" "${bams[@]}"
  # The existing runner logs failed samples and continues; propagate failures to
  # Docker/SLURM instead of letting a partial cohort look like a successful run.
  for bam in "${bams[@]}"; do
    sample=$(basename "$bam" .bam)
    if [ ! -e "$ROOT/work/done/$sample.cleanup" ] || [ -e "$ROOT/work/done/$sample.FAILED" ]; then
      echo "ERROR: $sample did not finish; see logs/samples/$sample.log" >&2
      failed=1
    fi
  done
  [ "$failed" -eq 0 ] || return 1
  report
}

command=${1:---help}
[ "$#" -eq 0 ] || shift
case "$command" in
  --help|-h|help) usage ;;
  cohort) cohort "$@" ;;
  sample)
    [ "$#" -eq 1 ] || { echo "Usage: haploid sample BAM" >&2; exit 2; }
    cohort "$1"
    ;;
  report) report ;;
  bias) exec "$PY" "$SC/cohort_bias.py" "$@" ;;
  refit)
    [ "$#" -gt 0 ] || { echo "Usage: haploid refit SAMPLE ..." >&2; exit 2; }
    bash "$SC/refit_downstream.sh" "$@"
    report
    ;;
  constrain) exec bash "$SC/refit_with.sh" "$@" ;;
  resolve) exec bash "$SC/resolve_discordant.sh" "$@" ;;
  run) exec bash "$SC/run_bam.sh" "$@" ;;
  controls) exec bash "$SC/build_controls.sh" "$@" ;;
  plot) exec "$PY" "$SC/plot_cna_compare.py" "$@" ;;
  check)
    for tool in samtools bcftools bgzip tabix liftOver readCounter runIchorCNA.R createPanelOfNormals.R Rscript "$PY"; do
      command -v "$tool" >/dev/null || { echo "ERROR: $tool missing" >&2; exit 1; }
    done
    "$PY" -c 'import numpy, pandas, matplotlib, pysam; print("Python packages: OK")'
    Rscript -e 'suppressMessages({library(ASCAT); library(DNAcopy); library(ichorCNA); library(optparse)}); stopifnot(packageVersion("ASCAT") >= "3.1.0"); stopifnot(file.exists(system.file("extdata", "gc_hg38_500kb.wig", package = "ichorCNA"))); cat("R packages: OK\n")'
    . "$SC/glimpse_env.sh"
    for tool in GLIMPSE2_chunk GLIMPSE2_split_reference GLIMPSE2_phase GLIMPSE2_ligate; do
      help=$(glimpse "$tool" --help 2>&1 || true)
      [[ "$help" == *--input* || "$help" == *--reference* || "$help" == *--bam-file* ]] \
        || { echo "ERROR: $tool is not runnable: $help" >&2; exit 1; }
    done
    echo 'Phasing tools: OK'
    echo 'HAPLOID dependencies: OK'
    ;;
  exec)
    [ "$#" -gt 0 ] || { echo "Usage: haploid exec COMMAND ..." >&2; exit 2; }
    exec "$@"
    ;;
  *) echo "ERROR: unknown command: $command (use --help)" >&2; exit 2 ;;
esac
