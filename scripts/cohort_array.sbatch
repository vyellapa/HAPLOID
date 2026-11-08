#!/bin/bash
# One SLURM array task per sample (scripts/run_sample.sh). Run scripts/run_cohort.sh once first (on 1-2 samples) so the
# preflight and the one-time shared prep (GLIMPSE2 reference) are done; every task here only does its own sample.
#   cd <project>; mkdir -p logs/array
#   ls <bam_dir>/*_500k/*.bam > logs/bam_list.txt
#   sbatch --export=ALL,PROJECT=$PWD,REF=/path/to/hg38.fa --array=1-$(wc -l < logs/bam_list.txt)%40 scripts/cohort_array.sbatch
# Finished samples exit immediately, so the same command resubmits only what is left. Afterwards:
#   python3 scripts/cohort_table.py, then the second pass that run_cohort.sh prints.
# Check memory / time of the first finished tasks with  seff <jobid>_<task>  and adjust --mem / -t.
#SBATCH -J haploid
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o logs/array/%A_%a.log
# activate your environment here (python3 with numpy/pandas/matplotlib, R with ASCAT + DNAcopy, bcftools, samtools)
set -euo pipefail
P=${PROJECT:?submit with --export=ALL,PROJECT=<project dir>,REF=<hg38.fa>}
b=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$P/logs/bam_list.txt")
s=$(basename "$b" .bam)
[ -e "$P/work/done/$s.cleanup" ] && { echo "$s already done"; exit 0; }
rm -f "$P/work/done/$s.FAILED"
export CORES=$SLURM_CPUS_PER_TASK PY=${PY:-python3} REF=${REF:?set REF}
export HAPLOID_RES=$P/resources KEEP_BAM=0
mkdir -p "$P/logs/samples"
bash "$P/scripts/run_sample.sh" "$s" "$b" >> "$P/logs/samples/$s.log" 2>&1
