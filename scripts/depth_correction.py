"""Bin-level logR (500 kb) from the fragment table, two ways, normalised with ichorCNA's healthy-donor PoN.

  ichor-style : total fragments per bin -> GC correction with a running-median curve over bin GC
                (ichorCNA uses loess on bin GC) -> log2 ratio to the median -> minus PoN median logR
  cnacs-style : per fragment-length class k, genome-wide fragment rate by fragment GC
                rate_k(g) = sum_b obs[b,k,g] / sum_b possible[b,k,g]; expected[b,k] = sum_g possible[b,k,g] rate_k(g);
                ratio[b] = sum_k obs[b,k] / sum_k expected[b,k] (CNACS correct_gc.sh + correct_length.sh)
                -> log2 ratio to the median -> minus PoN median logR. Per-class logR is also written.
Bins: those in the ichorCNA PoN (well-mappable autosomes), minus centromere +/- 2 Mb.
Noise metric: MAPD = median |logR[i+1] - logR[i]| over adjacent bins on the same chromosome.

Usage: depth_correction.py [samples...]   Output: work/depth/<sample>.logR.tsv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from allelic_lib import RES

ROOT = Path(__file__).resolve().parents[1]
DEP = ROOT / "work/depth"
SAMPLES = sys.argv[1:] or sys.exit("usage: depth_correction.py <samples...>")
BIN = 500_000
CEN_PAD = 2_000_000

pon = pd.read_csv(DEP / "ichor_pon_hg38_500kb.tsv", sep="\t")
pon = pon[pon.seqnames.astype(str).isin([str(i) for i in range(1, 23)])]
pon["chrom"] = "chr" + pon.seqnames.astype(str)
pon["bin"] = (pon.start - 1) // BIN
cb = pd.read_csv(RES / "cytoband_b38.bed", sep="\t", header=None)
cen = cb[cb[4] == "acen"].groupby(0).agg({1: "min", 2: "max"})
pon = pon[[not (cen.loc[c, 1] - CEN_PAD <= s <= cen.loc[c, 2] + CEN_PAD) for c, s in zip(pon.chrom, pon.start)]]
possible = np.load(DEP / "possible.npz")


def running_median_curve(x, y, w=0.02):
    """Smooth y against x with a sliding-window median (stand-in for loess)."""
    order = np.argsort(x); xs, ys = x[order], y[order]
    grid = np.linspace(np.quantile(xs, 0.01), np.quantile(xs, 0.99), 101)
    med = np.array([np.median(ys[np.abs(xs - g) <= w]) if np.any(np.abs(xs - g) <= w) else np.nan for g in grid])
    ok = np.isfinite(med)
    return np.interp(x, grid[ok], med[ok])


def mapd(df, col):
    d = []
    for _, g in df.groupby("chrom"):
        v = g.sort_values("bin")[col].to_numpy()
        d.append(np.abs(np.diff(v)))
    return float(np.nanmedian(np.concatenate(d)))


for s in SAMPLES:
    obs = np.load(DEP / f"{s}.frag_counts.npz")
    rows = []
    for _, p in pon.iterrows():
        o = obs[p.chrom][p.bin]; q = possible[p.chrom][p.bin]
        rows.append((p.chrom, p.bin, p.start, p.gc, p.Median, o, q))
    df = pd.DataFrame({"chrom": [r[0] for r in rows], "bin": [r[1] for r in rows], "start": [r[2] for r in rows],
                       "gc": [r[3] for r in rows], "pon": [r[4] for r in rows]})
    O = np.stack([r[5] for r in rows]); Q = np.stack([r[6] for r in rows])        # [bin, class, gc]

    # ichor-style
    tot = O.sum(axis=(1, 2)).astype(float)
    ok = tot > 0
    lr = np.full(len(df), np.nan); lr[ok] = np.log2(tot[ok] / np.median(tot[ok]))
    curve = running_median_curve(df.gc.to_numpy()[ok], lr[ok])
    lr[ok] = lr[ok] - curve
    df["logR_nopon"] = lr - np.nanmedian(lr)              # no PoN: fallback depth track in build_bins.py
    df["logR_ichorstyle"] = lr - np.nanmedian(lr) - df.pon

    # cnacs-style
    rate = O.sum(axis=0) / np.maximum(Q.sum(axis=0), 1)            # [class, gc]
    rate[Q.sum(axis=0) < 1000] = 0                                    # unreliable GC strata
    E = (Q * rate[None]).sum(axis=2)                                  # [bin, class]
    ratio = O.sum(axis=2).sum(axis=1) / np.maximum(E.sum(axis=1), 1e-9)
    lc = np.log2(np.where(ratio > 0, ratio, np.nan))
    df["logR_cnacs"] = lc - np.nanmedian(lc) - df.pon
    for k, name in enumerate(["lt130", "130_160", "160_190", "ge190"]):
        rk = O[:, k].sum(axis=1) / np.maximum(E[:, k], 1e-9)
        lk = np.log2(np.where(rk > 0, rk, np.nan))
        df[f"logR_cnacs_{name}"] = lk - np.nanmedian(lk) - df.pon
    df["frags"] = tot.astype(int)
    df.drop(columns="pon").to_csv(DEP / f"{s}.logR.tsv", sep="\t", index=False)
    fr = O.sum(axis=(0, 2)) / O.sum()
    print(f"{s}: {len(df)} bins  median frags/bin={int(np.median(tot))}  "
          f"MAPD ichor-style={mapd(df, 'logR_ichorstyle'):.4f}  cnacs-style={mapd(df, 'logR_cnacs'):.4f}  "
          f"length mix <130/130-160/160-190/>=190 = " + "/".join(f"{x:.2f}" for x in fr))
