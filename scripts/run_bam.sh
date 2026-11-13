#!/usr/bin/env bash
# Whole analysis from BAMs to final results: ichorCNA (unless given) -> run_cohort.sh (markdup ... joint fit, plots)
# -> resolve_discordant.sh -> final_report.py.
# Usage: run_bam.sh [--controls DIR] --bam X.bam [--sample NAME] [--ichor DIR] [--bam Y.bam ...]
#        run_bam.sh [--controls DIR] --sheet samples.tsv   (tab-separated with a header: sample, bam, and optionally ichor)
#   --controls DIR: the lab's own healthy-control references (build_controls.sh): its bin mask replaces the public one,
#                and ichorCNA run here normalises with its panel of normals (supplied --ichor results are used as given)
#   --ichor DIR: existing ichorCNA results for that BAM; the folder may hold one sample's files or several samples'
#                (<name>[_Y].correctedDepth.txt / .params.txt / .seg.txt, matched by the sample name). Used only at
#                500 kb bins (the bins every later step is built on); at any other bin size ichorCNA is run here.
#   env: CORES (4), JOBS (1), REF, HAPLOID_FINAL (work/final), HAPLOID_CONTROLS (same as --controls)
# Output: $HAPLOID_FINAL/summary.tsv, events.tsv, plots/, report.pdf; everything else under work/ and logs/
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); SC=$ROOT/scripts
. "$SC/bam_lib.sh"
export PY=${PY:-python3} CORES=${CORES:-4} JOBS=${JOBS:-1} PREP_JOBS=${PREP_JOBS:-1}
export REF=${REF:-$ROOT/resources/hg38.fa}
IN=$ROOT/work/inputs; mkdir -p "$IN" "$ROOT/logs"
die() { echo "ERROR: $*" >&2; exit 1; }
ts() { date '+%F %T'; }

names=(); bams=(); ichors=()
add() { names+=("$1"); bams+=("$2"); ichors+=("${3:-}"); }
while [ $# -gt 0 ]; do
  case $1 in
    --bam)    [ $# -ge 2 ] || die "--bam needs a path"; add "" "$2" ""; shift 2 ;;
    --sample) [ "${#bams[@]}" -gt 0 ] || die "--sample goes after its --bam"; names[${#names[@]}-1]=$2; shift 2 ;;
    --ichor)  [ "${#bams[@]}" -gt 0 ] || die "--ichor goes after its --bam"; ichors[${#ichors[@]}-1]=$2; shift 2 ;;
    --sheet)
      [ -s "$2" ] || die "sample sheet $2 not found"
      while IFS=$'\t' read -r s b i _; do
        [ -n "$s" ] && [ "$s" != sample ] && add "$s" "$b" "$i"
      done < <(tr -d '\r' < "$2"); shift 2 ;;
    --controls) [ -d "${2:-}" ] || die "--controls: folder ${2:-} not found"; HAPLOID_CONTROLS=$(cd "$2" && pwd); shift 2 ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
    *) die "unknown argument: $1 (see --help)" ;;
  esac
done
[ "${#bams[@]}" -gt 0 ] || die "no BAMs given (--bam or --sheet)"

# ---- healthy controls: their bin mask (read by every later step) and ichorCNA panel of normals ----
pon_note=""
if [ -n "${HAPLOID_CONTROLS:-}" ]; then
  [ -s "$HAPLOID_CONTROLS/control_bin_mask_500kb.tsv" ] \
    || die "$HAPLOID_CONTROLS has no control_bin_mask_500kb.tsv (build it with: controls --bam ...)"
  export HAPLOID_CONTROLS
  if [ -s "$HAPLOID_CONTROLS/ichor_pon_median.rds" ]; then
    export ICHOR_PON=$HAPLOID_CONTROLS/ichor_pon_median.rds; pon_note=" with panel of normals ($ICHOR_PON)"
  fi
  echo "[$(ts)] healthy controls: $HAPLOID_CONTROLS (bin mask${pon_note:+ + ichorCNA panel of normals})"
fi

ichor_files() {   # <dir> <sample> -> prints the correctedDepth path prefix ("<dir>/<name>.") or fails
  local d=$1 s=$2 hits
  [ -d "$d" ] || { echo "ichorCNA folder not found: $d" >&2; return 1; }
  mapfile -t hits < <(find "$d" -maxdepth 2 -name '*correctedDepth.txt' | sort)
  if [ "${#hits[@]}" -gt 1 ]; then
    mapfile -t hits < <(printf '%s\n' "${hits[@]}" | awk -v s="$s" '{n = $0; sub(/.*\//, "", n)} index(n, s) == 1')
  fi
  [ "${#hits[@]}" -eq 1 ] || { echo "ichorCNA: expected one *correctedDepth.txt for $s in $d, found ${#hits[@]}" >&2; return 1; }
  local p=${hits[0]%correctedDepth.txt}
  [ -s "${p}params.txt" ] && [ -s "${p}seg.txt" ] || { echo "ichorCNA: ${p}params.txt or ${p}seg.txt missing" >&2; return 1; }
  echo "$p"
}

bin_width() { awk -F'\t' 'NR == 2 {print $3 - $2 + 1; exit}' "$1"; }

# ---- stage each sample as work/inputs/<s>_500k/<s>.bam + <s>_Y.* (the layout run_cohort.sh reads) ----
staged=()
for k in "${!bams[@]}"; do
  b=$(cd "$(dirname "${bams[$k]}")" && pwd)/$(basename "${bams[$k]}")
  s=${names[$k]:-$(basename "$b" .bam)}
  [[ "$s" =~ ^[A-Za-z0-9._-]+$ ]] || die "sample name '$s' may only use letters, digits, '.', '_' and '-'"
  check_bam "$b" || exit 1
  d=$IN/${s}_500k; mkdir -p "$d"
  ln -sfn "$b" "$d/$s.bam"
  for x in "$b.bai" "${b%.bam}.bai" "$b.csi"; do [ -e "$x" ] && { ln -sfn "$x" "$d/$s.bam.${x##*.}"; break; }; done
  src=""
  if [ -n "${ichors[$k]}" ]; then
    p=$(ichor_files "${ichors[$k]}" "$s") || die "sample $s: see above"
    w=$(bin_width "${p}correctedDepth.txt")
    if [ "$w" = 500000 ]; then
      for f in correctedDepth params seg; do                  # copied (small), so reruns need only the BAM mount
        rm -f "$d/${s}_Y.$f.txt"; cp "${p}$f.txt" "$d/${s}_Y.$f.txt"; done
      src="supplied ($(dirname "$p"))"
    else
      echo "[$(ts)] $s: supplied ichorCNA bins are $w bp, not 500 kb - running ichorCNA at 500 kb instead"
      src="computed$pon_note (supplied results at $w bp bins not usable)"
    fi
  fi
  if [ -z "$src" ] || [[ "$src" == computed* ]]; then
    src=${src:-computed$pon_note}         # rerun when an earlier run here used a different panel of normals (or none)
    if [ ! -s "$d/${s}_Y.correctedDepth.txt" ] || [ "$(cat "$d/ichor_source.txt" 2>/dev/null)" != "$src" ]; then
      rm -f "$d/${s}_Y."*.txt
      echo "[$(ts)] $s: ichorCNA"
      bash "$SC/run_ichor.sh" "$s" "$d/$s.bam" "$d" > "$ROOT/logs/$s.ichor.log" 2>&1 \
        || die "ichorCNA failed for $s (see logs/$s.ichor.log)"
      tail -1 "$ROOT/logs/$s.ichor.log"
    fi
  fi
  echo "$src" > "$d/ichor_source.txt"
  staged+=("$d/$s.bam")
done

# ---- pipeline, then the ichorCNA-ploidy check, then final results ----
BAM_DIR=$IN HAPLOID_NO_NEXT=1 bash "$SC/run_cohort.sh" "${staged[@]}"
failed=0; ok=()
for b in "${staged[@]}"; do
  s=$(basename "$b" .bam)
  if [ -e "$ROOT/work/done/$s.cleanup" ]; then ok+=("$s"); else failed=1; echo "ERROR: $s did not finish (logs/samples/$s.log)" >&2; fi
done
[ "${#ok[@]}" -gt 0 ] || exit 1
BAM_DIR=$IN JOBS=$JOBS bash "$SC/resolve_discordant.sh" > "$ROOT/logs/resolve_discordant.log" 2>&1 \
  || echo "  resolve_discordant.sh failed (logs/resolve_discordant.log); reporting the free fits" >&2
(cd "$SC" && "$PY" final_report.py "${ok[@]}")
exit "$failed"
