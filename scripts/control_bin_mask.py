"""Bin mask + per-bin bias correction from healthy controls processed exactly like the samples.

Useful when the samples' ichorCNA correctedDepth was made WITHOUT a panel of normals. With bwa alignment to an hg38 that
has _alt contigs, reads in alt-haplotype regions (MHC, ...) get MAPQ 0 and fail readCounter's MAPQ >= 20, and other
reproducible per-bin biases remain. With healthy controls run through the same pipeline:
  bias[bin] = median GC-corrected logR across those controls
(with --pon <wigs>, the median of a second control set is subtracted instead: use only if the samples' depth was already
normalised with a PoN built from those wigs).
Each control: counts -> map >= 0.9 bins -> log2(count / median) -> GC correction (running median) -> centre.
Both sets are summed to the coarser bin size first; values are then assigned to every 500 kb bin inside.
Masked (reason column): |bias| > BIAS_MAX (biased_in_controls), sample-pipeline controls variable between people
(robust SD > SD_FOLD x genome median: variable_in_controls, e.g. copy-number polymorphisms), or no coverage.
Usage: control_bin_mask.py --pipeline <sample-pipeline control wigs...> [--pon <PoN-pipeline control wigs...>] [--out DIR]
Output: DIR (default resources/)/control_bin_mask_500kb.tsv (masked bins), pipeline_bias_500kb.tsv (bias for all bins)
ichorCNA's hg38 GC / mappability wigs: $ICHOR_EXTDATA, else ichorCNA/inst/extdata in the project, else the R package's.
build_controls.sh runs this on control BAMs counted exactly like the samples.
"""
import argparse, os, subprocess
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXT = Path(os.environ.get("ICHOR_EXTDATA", ROOT / "ichorCNA/inst/extdata"))
if not (EXT / "gc_hg38_500kb.wig").exists():
    EXT = Path(subprocess.run(["Rscript", "-e", 'cat(system.file("extdata", package = "ichorCNA"))'],
                              capture_output=True, text=True).stdout.strip())
SD_FOLD, BIAS_MAX, BIN = 3.0, 0.3, 500_000
AUTO = [str(i) for i in range(1, 23)]


def read_wig(f):
    rows, chrom, start, step, i = [], None, 1, BIN, 0
    for line in open(f):
        if line.startswith("fixedStep"):
            kv = dict(x.split("=") for x in line.split()[1:])
            chrom, start, step, i = kv["chrom"].replace("chr", ""), int(kv["start"]), int(kv["step"]), 0
        else:
            rows.append((chrom, start + i * step, float(line), step)); i += 1
    return pd.DataFrame(rows, columns=["chr", "start", "v", "step"])


def coarsen(w, step):
    w = w.assign(start=(w.start - 1) // step * step + 1)
    return w.groupby(["chr", "start"], sort=False, as_index=False).v.sum()


def curve(x, y, w=0.02):
    o = np.argsort(x); xs, ys = x[o], y[o]
    g = np.linspace(np.quantile(xs, .01), np.quantile(xs, .99), 101)
    m = np.array([np.median(ys[np.abs(xs - v) <= w]) if np.sum(np.abs(xs - v) <= w) >= 20 else np.nan for v in g])
    ok = np.isfinite(m)
    return np.interp(x, g[ok], m[ok])


def control_logr(wigs, step):
    kb = step // 1000
    gc = read_wig(EXT / f"gc_hg38_{kb}kb.wig").rename(columns={"v": "gc"})
    mp = read_wig(EXT / f"map_hg38_{kb}kb.wig").rename(columns={"v": "map"})
    base = gc[["chr", "start", "gc"]].merge(mp[["chr", "start", "map"]], on=["chr", "start"])
    base = base[base.chr.isin(AUTO) & (base["map"] >= 0.9) & (base.gc > 0)].reset_index(drop=True)
    L = []
    for f in wigs:
        d = base.merge(coarsen(read_wig(f), step), on=["chr", "start"], how="left")
        x = np.log2(np.where(d.v > 0, d.v, np.nan) / np.nanmedian(d.v)); ok = np.isfinite(x)
        x[ok] -= curve(d.gc.to_numpy()[ok], x[ok])
        L.append(x - np.nanmedian(x))
    M = np.vstack(L)
    med = np.nanmedian(M, axis=0)
    return base.assign(med=med, sd=1.4826 * np.nanmedian(np.abs(M - med), axis=0))


ap = argparse.ArgumentParser(description="healthy-control bin mask (see the module docstring)")
ap.add_argument("--pipeline", nargs="+", required=True, help="control wigs counted like the samples")
ap.add_argument("--pon", nargs="+", default=[], help="control wigs of the PoN the samples were normalised with")
ap.add_argument("--out", type=Path, default=ROOT / "resources", help="output folder (default resources/)")
args = ap.parse_args()
pipe, pon, OUT = args.pipeline, args.pon, args.out
OUT.mkdir(parents=True, exist_ok=True)
step = max(int(read_wig(f).step.iat[0]) for f in pipe[:1] + pon[:1])
m = control_logr(pipe, step)
if pon:
    m = m.merge(control_logr(pon, step)[["chr", "start", "med"]], on=["chr", "start"], suffixes=("", "_pon"))
    m["bias"] = m.med - m.med_pon
else:
    m["bias"] = m.med
sd0 = np.nanmedian(m.sd)
m["reason"] = np.where(~np.isfinite(m.bias), "no_coverage",
                       np.where(np.abs(m.bias) > BIAS_MAX, "biased_in_controls",
                                np.where(m.sd > SD_FOLD * sd0, "variable_in_controls", "")))
# expand to 500 kb bins
k = step // BIN
x = pd.concat([m.assign(start=m.start + i * BIN) for i in range(k)]).sort_values(["chr", "start"])
x["chrom"], x["end"] = "chr" + x["chr"], x["start"] + BIN - 1
x[["chrom", "start", "end", "bias", "sd"]].round(4).to_csv(OUT / "pipeline_bias_500kb.tsv", sep="\t", index=False)
x[x.reason != ""][["chrom", "start", "end", "reason", "bias", "sd"]].round(4).to_csv(
    OUT / "control_bin_mask_500kb.tsv", sep="\t", index=False)
ab = np.abs(m.bias[np.isfinite(m.bias)])
print(f"{len(pipe)} sample-pipeline controls{f' minus {len(pon)} PoN-pipeline controls' if pon else ''}, {len(m)} bins of {step // 1000} kb; "
      f"|bias| median {ab.median():.3f}, 90% {ab.quantile(.9):.3f}, >0.1 in {int((ab > 0.1).sum())} bins; "
      f"masked {int((x.reason != '').sum())} x 500 kb: {x[x.reason != ''].reason.value_counts().to_dict()}")
