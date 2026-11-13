# Shared by run_bam.sh, run_ichor.sh and build_controls.sh (sourced), so samples and healthy controls are checked and
# counted identically: a control processed differently from the samples adds bias instead of removing it.

ICHOR_CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

check_bam() {   # <bam>: indexed, chr-named, and chr1-22 lengths match $REF.fai; prints the problem and fails otherwise
  local b=$1 bad
  [ -s "$b" ] || { echo "BAM not found: $b" >&2; return 1; }
  [ -e "$b.bai" ] || [ -e "${b%.bam}.bai" ] || [ -e "$b.csi" ] || { echo "BAM not indexed: $b (samtools index it)" >&2; return 1; }
  bad=$(samtools view -H "$b" | awk -F'\t' '$1 == "@SQ" {for (i = 2; i <= NF; i++) {if ($i ~ /^SN:/) n = substr($i, 4); if ($i ~ /^LN:/) l = substr($i, 4)}; print n "\t" l}' \
        | awk -F'\t' 'NR == FNR {len[$1] = $2; next} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {n++; if (len[$1] != $2) bad = bad " " $1} END {if (n < 22) print "missing chr1-22 (" n " found; chr-prefixed names needed)"; else if (bad) print "lengths differ from the reference on" bad}' \
          <(cut -f1,2 "$REF.fai") -)
  [ -z "$bad" ] || { echo "$b: $bad - align to hg38 (chr1, chr2, ...)" >&2; return 1; }
}

count_reads() {   # <bam> <out.wig>: HMMcopy readCounter, 500 kb bins, MAPQ >= 20, chr1-22, X, Y (ichorCNA's input)
  "${READCOUNTER:-readCounter}" --window 500000 --quality 20 --chromosome "$ICHOR_CHRS" "$1" > "$2.tmp" && mv "$2.tmp" "$2"
}

ichor_extdata() {   # ichorCNA's hg38 GC / mappability / centromere files (installed R package)
  local e
  e=$(Rscript -e 'cat(system.file("extdata", package = "ichorCNA"))')
  [ -s "$e/gc_hg38_500kb.wig" ] || { echo "ichorCNA hg38 500 kb files not found in '$e'" >&2; return 1; }
  echo "$e"
}
