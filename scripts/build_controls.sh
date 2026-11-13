#!/usr/bin/env bash
# A lab's own healthy-control references, from healthy-control BAMs (or their 500 kb readCounter wigs):
#   ichor_pon_median.rds (+ .txt)  ichorCNA panel of normals (createPanelOfNormals.R): ichorCNA normalises each sample's
#                                  depth with it when ichorCNA is run here (run_bam.sh --controls)
#   control_bin_mask_500kb.tsv     bins that misbehave in the controls (control_bin_mask.py), left out of the fit;
#                                  replaces the public bin mask
#   pipeline_bias_500kb.tsv        per-bin median control logR (subtracted only with HAPLOID_CONTROL_BIAS=1)
#   controls.tsv, wigs/            the controls used and their read counts
# The controls must be processed exactly like the samples: same aligner, same hg38 (with or without _alt contigs), same
# read filtering. BAMs are counted here with the samples' readCounter settings (bam_lib.sh). Healthy people only; >= 10
# recommended (at least 3 required).
# Usage: build_controls.sh [--bam CTRL.bam ...] [--wig CTRL.wig ...] [--list FILE] [--out DIR]
#   --list FILE: one BAM or wig path per line.  Control names are the file names without .bam / .wig.
#   env: JOBS (BAMs counted in parallel, 4), REF, HAPLOID_CONTROLS_OUT (default output folder, work/controls)
# Then: run_bam.sh --controls DIR ...   (run_cohort.sh: HAPLOID_CONTROLS=DIR, ichorCNA PoN not applied there)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); SC=$ROOT/scripts
. "$SC/bam_lib.sh"
export REF=${REF:-$ROOT/resources/hg38.fa}
PY=${PY:-python3}; JOBS=${JOBS:-4}
OUT=${HAPLOID_CONTROLS_OUT:-$ROOT/work/controls}
die() { echo "ERROR: $*" >&2; exit 1; }
ts() { date '+%F %T'; }

files=()
while [ $# -gt 0 ]; do
  case $1 in
    --bam|--wig) [ $# -ge 2 ] || die "$1 needs a path"; files+=("$2"); shift 2 ;;
    --list) [ -s "${2:-}" ] || die "list ${2:-} not found"
            while IFS= read -r f; do f=${f%$'\r'}; [ -n "$f" ] && files+=("$f"); done < "$2"; shift 2 ;;
    --out)  [ $# -ge 2 ] || die "--out needs a folder"; OUT=$2; shift 2 ;;
    -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
    *) die "unknown argument: $1 (see --help)" ;;
  esac
done
[ "${#files[@]}" -ge 3 ] || die "${#files[@]} controls given; at least 3 are needed (>= 10 recommended)"
[ "${#files[@]}" -ge 10 ] || echo "WARNING: ${#files[@]} controls; with fewer than 10 the panel and the variability mask are noisy" >&2
mkdir -p "$OUT/wigs"; OUT=$(cd "$OUT" && pwd)

# ---- 1. one 500 kb wig per control (BAMs counted like the samples; wigs checked for 500 kb, chr-named bins) ----
names=(); todo=()
for f in "${files[@]}"; do
  case $f in *.bam) n=$(basename "$f" .bam) ;; *.wig) n=$(basename "$f" .wig) ;; *) die "$f: expected a .bam or .wig" ;; esac
  [[ "$n" =~ ^[A-Za-z0-9._-]+$ ]] || die "control name '$n' may only use letters, digits, '.', '_' and '-'"
  for m in ${names[@]+"${names[@]}"}; do [ "$m" != "$n" ] || die "two controls are named $n"; done
  names+=("$n")
  w=$OUT/wigs/$n.wig
  if [[ $f == *.bam ]]; then
    check_bam "$f" || die "control $n: see above"
    [ -s "$w" ] && [ "$w" -nt "$f" ] || todo+=("$f")
  else
    [ -s "$f" ] || die "wig not found: $f"
    head -1 "$f" | grep -q '^fixedStep chrom=chr[0-9XY]* start=1 step=500000 ' \
      || die "$f: not a 500 kb readCounter wig with chr-named chromosomes ($(head -1 "$f"))"
    cp "$f" "$w"
  fi
done
if [ "${#todo[@]}" -gt 0 ]; then
  echo "[$(ts)] counting reads in ${#todo[@]} control BAMs ($JOBS at a time)"
  printf '%s\n' "${todo[@]}" | xargs -P "$JOBS" -I{} bash -c \
    '. "$1/bam_lib.sh"; b=$2; count_reads "$b" "$3/$(basename "$b" .bam).wig" || { echo "readCounter failed: $b" >&2; exit 255; }' \
    _ "$SC" {} "$OUT/wigs" || die "read counting failed"
fi
wigs=(); for n in "${names[@]}"; do wigs+=("$OUT/wigs/$n.wig"); done
for w in "${wigs[@]}"; do [ -s "$w" ] || die "$w missing"; done
# sex as ichorCNA infers it (fraction of reads on chrY >= 0.001). Only chrX/chrY depend on it: the bin mask uses chr1-22,
# the panel centres each control's chrX, and HAPLOID fits chr1-22 only. Panel chrY = median over the controls' mix.
yfrac() { awk '/^fixedStep/ {split($2, a, "="); c = a[2]; next} {t += $1; if (c == "chrY") y += $1} END {printf "%.5f", t ? y / t : 0}' "$1"; }
{ printf 'control\tsource\tchrY_read_fraction\tsex\n'
  for k in "${!names[@]}"; do
    f=$(yfrac "${wigs[$k]}"); printf '%s\t%s\t%s\t%s\n' "${names[$k]}" "${files[$k]}" "$f" "$(awk -v f="$f" 'BEGIN {print (f >= 0.001) ? "male" : "female"}')"
  done; } > "$OUT/controls.tsv"
n_male=$(awk -F'\t' 'NR > 1 && $4 == "male"' "$OUT/controls.tsv" | wc -l | tr -d ' ')
echo "[$(ts)] ${#wigs[@]} controls: $n_male male, $((${#wigs[@]} - n_male)) female"

# ---- 2. ichorCNA panel of normals (chr1-22, X, Y: the chromosomes run_ichor.sh analyses) ----
EXT=$(ichor_extdata) || die "ichorCNA not installed"
echo "[$(ts)] ichorCNA panel of normals from ${#wigs[@]} controls"
printf '%s\n' "${wigs[@]}" > "$OUT/wigs.txt"
rm -f "$OUT"/ichor_pon_median.*
Rscript "${ICHOR_PON_SCRIPT:-$(command -v createPanelOfNormals.R)}" --filelist "$OUT/wigs.txt" \
  --gcWig "$EXT/gc_hg38_500kb.wig" --mapWig "$EXT/map_hg38_500kb.wig" \
  --centromere "$EXT/GRCh38.GCA_000001405.2_centromere_acen.txt" --chrs 'c(1:22, "X", "Y")' \
  --outfile "$OUT/ichor_pon" > "$OUT/ichor_pon.log" 2>&1 \
  || die "createPanelOfNormals.R failed (see $OUT/ichor_pon.log)"
[ -s "$OUT/ichor_pon_median.rds" ] || die "no panel written (see $OUT/ichor_pon.log)"

# ---- 3. HAPLOID control bin mask ----
echo "[$(ts)] control bin mask"
ICHOR_EXTDATA=$EXT "$PY" "$SC/control_bin_mask.py" --pipeline "${wigs[@]}" --out "$OUT" | tee "$OUT/control_bin_mask.log"

cat << DONE
[$(ts)] controls ready in $OUT:
  ichor_pon_median.rds        ichorCNA panel of normals (${#wigs[@]} controls: $n_male male, $((${#wigs[@]} - n_male)) female)
  control_bin_mask_500kb.tsv  $(($(wc -l < "$OUT/control_bin_mask_500kb.tsv") - 1)) masked 500 kb bins
Use them with: run --controls <this folder> ...  (it replaces the public bin mask; supplied --ichor results are used as given)
DONE
