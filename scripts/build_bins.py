"""Join 500 kb depth logR with per-bin phased BAF for the CNV pipelines.

Depth (column `logR`, source in `logR_source`):
  ichorCNA  : <sample>_Y.correctedDepth.txt (log2_TNratio_corrected: GC + mappability corrected; the cohort ran without a PoN) - default;
              path from ICHOR_DEPTH ("{s}" = sample). Used as is: the healthy-plasma control bias
              (resources/pipeline_bias_500kb.tsv) is NOT subtracted by default (env HAPLOID_CONTROL_BIAS=1 to subtract it) -
              its chromosome-scale part does not transfer to CSF (chr19 -0.11 in the plasma controls vs -0.03 in tumor-free
              CSF, so subtracting it made a false chr19 gain); bin-level control artifacts are masked instead (below).
  fragments : fallback when that file is missing - our fragment counts, GC-corrected, no PoN
              (work/depth/<sample>.logR.tsv from fragment_counts.py + depth_correction.py)
  cohort    : when resources/cohort_bias_500kb.tsv exists (cohort_bias.py, built from this cohort's tumor-negative samples),
              its residual per-bin bias is subtracted too and its flagged bins dropped; `logR_precohort` keeps the value
              before that step (cohort_bias.py reads it, so rebuilding the table does not compound).
Bins within +/- 2 Mb of a centromere, bins in resources/exclusion/bins_500kb.tsv (GIAB segdups / VDJ + UCSC gap
exclusions covering >= 50% of the bin, no padding; build_exclusion.py), and bins in resources/control_bin_mask_500kb.tsv (built by control_bin_mask.py from
bwa-aligned healthy controls: e.g. MHC and other alt-contig regions where reads get MAPQ 0), are dropped.
BAF per bin = H1 read fraction over GLIMPSE2-phased hets (phase orientation is arbitrary between switch errors,
so pipelines must mirror it, as ASCAT does for heterozygous probes). Bins with < MIN_READS haplotype reads
get BAF = NA.
Usage: build_bins.py [samples...]   Output: work/cnv/<sample>.bins.tsv
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
from allelic_lib import CEN_PAD, CONTROLS, RES, bin_mask, exclusion_bins

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = sys.argv[1:] or sys.exit("usage: build_bins.py <samples...>")
BIN, MIN_READS = 500_000, 100
ICHOR_DEPTH = os.environ.get("ICHOR_DEPTH", str(ROOT / "{s}_500k/{s}_Y.correctedDepth.txt"))
OUT = ROOT / "work/cnv"; OUT.mkdir(parents=True, exist_ok=True)
cb = pd.read_csv(RES / "cytoband_b38.bed", sep="\t", header=None)
cen = cb[cb[4] == "acen"].groupby(0).agg({1: "min", 2: "max"})


def depth(s):
    f = Path(ICHOR_DEPTH.replace("{s}", s))
    if f.exists():
        d = pd.read_csv(f, sep="\t")
        d = d[d.chr.astype(str).isin([str(i) for i in range(1, 23)])].copy()
        d["chrom"] = "chr" + d.chr.astype(str)
        d["bin"] = (d.start - 1) // BIN
        d = d.rename(columns={"log2_TNratio_corrected": "logR"})[["chrom", "bin", "start", "logR"]]
        src = "ichorCNA"
        bf = CONTROLS / "pipeline_bias_500kb.tsv"
        if os.environ.get("HAPLOID_CONTROL_BIAS", "0") == "1" and bf.exists():
            d = d.merge(pd.read_csv(bf, sep="\t")[["chrom", "start", "bias"]], on=["chrom", "start"], how="left")
            d["logR"] = d.logR - d.bias.fillna(0)
            src = "ichorCNA-biascorrected"
    else:
        d = pd.read_csv(ROOT / f"work/depth/{s}.logR.tsv", sep="\t")
        d = d.rename(columns={"logR_nopon": "logR"})[["chrom", "bin", "start", "logR"]]
        src = "fragments"
    d["logR_precohort"] = d.logR
    cf = Path(os.environ.get("HAPLOID_COHORT_BIAS", ROOT / "resources/cohort_bias_500kb.tsv"))
    if cf.exists() and src != "fragments":
        cbias = pd.read_csv(cf, sep="\t")
        d = d.merge(cbias[["chrom", "start", "cohort_bias", "cohort_mask"]], on=["chrom", "start"], how="left")
        d = d[~d.cohort_mask.fillna(False).astype(bool)]
        d["logR"] = d.logR - d.cohort_bias.fillna(0)
        d = d.drop(columns=["cohort_bias", "cohort_mask"]); src += "+cohort"
    near_cen = [cen.loc[c, 1] - CEN_PAD <= st <= cen.loc[c, 2] + CEN_PAD for c, st in zip(d.chrom, d.start)]
    key = d.set_index(["chrom", "start"]).index
    masked = key.isin(bin_mask().set_index(["chrom", "start"]).index) | \
        key.isin(exclusion_bins().set_index(["chrom", "start"]).index)
    return d[~np.array(near_cen) & ~masked & np.isfinite(d.logR)].assign(logR_source=src), src


for s in SAMPLES:
    lr, src = depth(s)
    hc = pd.concat([pd.read_csv(f, sep="\t") for f in sorted((ROOT / "work/phased").glob(f"{s}.chr*.hapcounts.tsv.gz"))])
    hc["bin"] = hc.pos // BIN
    g = hc.groupby(["chrom", "bin"]).agg(h1=("h1", "sum"), h2=("h2", "sum"), hets=("pos", "size")).reset_index()
    g["n"] = g.h1 + g.h2
    d = lr.merge(g, on=["chrom", "bin"], how="left")
    d[["h1", "h2", "n", "hets"]] = d[["h1", "h2", "n", "hets"]].fillna(0).astype(int)
    d["baf"] = np.where(d.n >= MIN_READS, d.h1 / d.n.where(d.n > 0, 1), np.nan)
    d["pos"] = d.start + BIN // 2
    d.to_csv(OUT / f"{s}.bins.tsv", sep="\t", index=False)
    print(f"{s}: {len(d)} bins (depth: {src}), {d.baf.notna().sum()} with BAF, median reads/bin={int(d.n.median())}")
