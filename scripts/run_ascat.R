# Pipeline B: real ASCAT (3.2.0) on 500 kb bins (logR = CNACS-style corrected depth, BAF = phased H1 fraction),
# plus DNAcopy CBS segmentation of logR and mirrored BAF for pipeline C.
# Every bin is treated as a germline heterozygous "probe" (BAF comes only from phased het sites), so ASCAT
# mirrors it and ASPCF segments logR + BAF jointly. gamma = 1 for sequencing.
# Usage: Rscript run_ascat.R <sample> <gender XX|XY> [logR column] [aspcf penalty]
suppressMessages({library(ASCAT); library(DNAcopy)})
a <- commandArgs(trailingOnly = TRUE)
s <- a[1]; gender <- a[2]
lrcol <- ifelse(length(a) >= 3, a[3], "logR_cnacs")
penalty <- ifelse(length(a) >= 4, as.numeric(a[4]), 70)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
out <- file.path(root, "work/cnv", paste0(s, "_ascat")); dir.create(out, showWarnings = FALSE, recursive = TRUE)
f_seg <- file.path(out, paste0(s, ".ascat_segments.tsv"))
unlink(f_seg)   # written only when this run's ASCAT fit succeeds; an earlier run's would hide the ASPCF fallback

d <- read.delim(file.path(root, "work/cnv", paste0(s, ".bins.tsv")))
d <- d[is.finite(d[[lrcol]]), ]
d$chrs <- sub("chr", "", d$chrom)
id <- paste0("b", seq_len(nrow(d)))
lr <- data.frame(chrs = d$chrs, pos = d$pos, x = d[[lrcol]], row.names = id); names(lr)[3] <- s
bf <- data.frame(chrs = d$chrs, pos = d$pos, x = d$baf, row.names = id); names(bf)[3] <- s
f_lr <- file.path(out, "Tumor_LogR.txt"); f_bf <- file.path(out, "Tumor_BAF.txt")
write.table(lr, f_lr, sep = "\t", quote = FALSE, col.names = NA)
write.table(bf, f_bf, sep = "\t", quote = FALSE, col.names = NA)

obj <- ascat.loadData(Tumor_LogR_file = f_lr, Tumor_BAF_file = f_bf, chrs = c(1:22), gender = gender,
                      genomeVersion = "hg38")
gg <- list(germlinegenotypes = matrix(!is.finite(d$baf), ncol = 1, dimnames = list(id, s)),  # TRUE = homozygous
           failedarrays = NULL)
obj <- ascat.aspcf(obj, ascat.gg = gg, penalty = penalty, out.dir = out, seed = 42L)
ascat.plotSegmentedData(obj, img.dir = out)
# ASPCF breakpoints alone, written before the ASCAT fit (which finds no solution when the tumor fraction is ~0);
# joint_fit.py falls back to these when the ASCAT segments are missing
lrs <- obj$Tumor_LogR_segmented[, 1]; ch <- as.character(obj$SNPpos[, 1]); ps <- obj$SNPpos[, 2]
run <- cumsum(c(TRUE, ch[-1] != ch[-length(ch)] | lrs[-1] != lrs[-length(lrs)]))
aspcf <- do.call(rbind, lapply(split(seq_along(run), run), function(i)
  data.frame(chr = ch[i[1]], startpos = min(ps[i]), endpos = max(ps[i]), logR = lrs[i[1]], bins = length(i))))
write.table(aspcf, file.path(out, paste0(s, ".aspcf_segments.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
res <- tryCatch(ascat.runAscat(obj, gamma = 1, img.dir = out, write_segments = FALSE, min_purity = 0.05),
                error = function(e) { message("ASCAT fit failed: ", conditionMessage(e)); NULL })
seg <- if (is.null(res)) NULL else res$segments
if (!is.null(seg) && nrow(seg) > 0)
  write.table(seg, f_seg, sep = "\t", quote = FALSE, row.names = FALSE)
na1 <- function(x) if (is.null(x) || length(x) == 0) NA else x
summ <- data.frame(sample = s, purity = na1(res$purity), ploidy = na1(res$ploidy), goodness = na1(res$goodnessOfFit),
                   penalty = penalty, logR = lrcol)
write.table(summ, file.path(out, paste0(s, ".ascat_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

# --- CBS for pipeline C (CNACS segments depth and BAF separately with DNAcopy) ---
cbs <- function(v, name, track) {   # NULL (track skipped) when too few bins: DNAcopy's noise estimate needs >= 2
  ok <- is.finite(v)
  if (sum(ok) < 50) { message("CBS: only ", sum(ok), " finite ", track, " bins - track skipped"); return(NULL) }
  cna <- CNA(v[ok], d$chrs[ok], d$pos[ok], data.type = "logratio", sampleid = name)
  sg <- segment(smooth.CNA(cna), alpha = 0.01, undo.splits = "prune", undo.prune = 0.01, min.width = 2, verbose = 0)
  cbind(track = track, sg$output)
}
mb <- abs(d$baf - 0.5) + 0.5
write.table(rbind(cbs(d[[lrcol]], s, "logR"), cbs(mb, s, "mBAF")),
            file.path(root, "work/cnv", paste0(s, ".cbs_segments.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
