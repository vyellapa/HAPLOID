"""Per-haplotype read counts at GLIMPSE2-phased heterozygous sites.

Joins the ligated GLIMPSE2 output (phased GT + genotype posteriors) with the allele depths from our own
pileup (work/gl/<sample>.chr<N>.gl.bcf: dup-filtered, MAPQ/BQ >= 20), keeps confident hets
(P(het) >= MIN_GP_HET), and assigns reads to haplotypes: GT 0|1 -> H1 carries REF, 1|0 -> H1 carries ALT.
Adds the genetic position (cM) for the phase-switch model.

Usage: haplotype_counts.py <sample> <chrom number> [min_gp_het]
Output: work/phased/<sample>.chr<N>.hapcounts.tsv.gz
"""
import subprocess, sys
import numpy as np, pandas as pd
from allelic_lib import ROOT, bq, in_mask

S, CHR = sys.argv[1], int(sys.argv[2])
MIN_GP_HET = float(sys.argv[3]) if len(sys.argv) > 3 else 0.9

ph = bq(ROOT / f"work/glimpse/{S}/{S}.chr{CHR}.phased.bcf", "%POS\t%REF\t%ALT\t[%GT]\t[%GP]\n")
ph.columns = ["pos", "ref", "alt", "gt", "gp"]
ph = ph[ph["gt"].isin(["0|1", "1|0"])].copy()
ph["gp_het"] = ph.gp.str.split(",", expand=True)[1].astype(float)
ph = ph[ph.gp_het >= MIN_GP_HET]

ad = bq(ROOT / f"work/gl/{S}.chr{CHR}.gl.bcf", "%POS\t%REF\t%ALT\t[%AD]\n").dropna()
ad.columns = ["pos", "ref", "alt", "ad"]
ad[["r", "a"]] = ad["ad"].str.split(",", expand=True).astype(int)

d = ph.merge(ad[["pos", "ref", "alt", "r", "a"]], on=["pos", "ref", "alt"])
d = d[~in_mask(CHR, d.pos.to_numpy()) & ((d.r + d.a) > 0)].copy()
d["h1"] = np.where(d["gt"] == "0|1", d.r, d.a)
d["h2"] = np.where(d["gt"] == "0|1", d.a, d.r)

gm = pd.read_csv(ROOT / f"resources/maps_b38/chr{CHR}.b38.gmap.gz", sep="\t")
d["cm"] = np.interp(d.pos, gm["pos"], gm["cM"])
d.insert(0, "chrom", f"chr{CHR}")

out = ROOT / "work/phased"; out.mkdir(parents=True, exist_ok=True)
f = out / f"{S}.chr{CHR}.hapcounts.tsv.gz"
d[["chrom", "pos", "ref", "alt", "gt", "gp_het", "r", "a", "h1", "h2", "cm"]].to_csv(f, sep="\t", index=False)
n = d.r + d.a
print(f"{S} chr{CHR}: phased hets={len(ph)}  with reads={len(d)}  reads={int(n.sum())}  "
      f"mean reads/het={n.mean():.2f}  H1 fraction={d.h1.sum() / n.sum():.3f}  -> {f.name}")
