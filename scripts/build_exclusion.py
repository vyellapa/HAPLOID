"""Build the region exclusions from published hg38 resources (replaces the earlier hand-made CNV exclusion BED).

Sources (copied into resources/exclusion/ with their provenance):
  GIAB genome stratifications (NIST)  GRCh38_segdups.bed.gz           segmental duplications (from UCSC genomicSuperDups)
                                      GRCh38_AllTandemRepeats_ge101bp_slop5.bed.gz  tandem repeats >= 101 bp (~ read
                                      length: reads cannot span them, so alleles inside are unreliable)
  (deliberately NOT used: GRCh38_lowmappabilityall, GRCh38_alldifficultregions - too broad, 220 / 572 Mb)
  GIAB HG002 SV v0.6                  vdj_loci.GRCh38.bed             IG / TCR loci (somatic V(D)J rearrangement)
  UCSC hg38 database                  gap.txt.gz, centromeres.txt.gz  telomeres, heterochromatin (gap table) and centromeres
                                      (per chromosome, first to last modelled centromere piece) -> ucsc_gap_centromeres_hg38.bed
                                      (replaces DRAGEN 4.4 human.hg38.excl.bed: same regions on chr1-22, public)
Outputs:
  resources/exclusion/sites_hg38.bed.gz  merged intervals for allele sites: segdups + tandem repeats >= 101 bp + VDJ + UCSC gap,
                                         each +/- SITE_PAD (50 bp: a read starting just outside can still mismap)
  resources/exclusion/bins_500kb.tsv     500 kb bins whose coverage by segdups + VDJ + UCSC gap is >= BIN_FRAC (no padding:
                                         neighbouring bins are depth-corrected and mappability-filtered on their own) (tandem repeats are short and scattered, and
                                         mappability is handled per bin by ichorCNA's map >= 0.9 filter, so neither is used)
  resources/exclusion/public_bin_mask_500kb.tsv  depth/allele bin mask from public data, used by allelic_lib.bin_mask() only when no
                                         healthy-control mask is given: bins touched by an hg38 alternate-haplotype locus (UCSC
                                         altLocations.txt.gz; reads can map to the _alt contig with MAPQ 0, so depth is biased) and
                                         the 1% most GC-rich autosomal bins (ichorCNA gc_hg38_500kb.wig; GC bias varies between
                                         libraries and per-sample correction leaves residue there)
  resources/exclusion/review_*.tsv       per-source size, panel sites lost, and every masked bin block (printed)
Usage: build_exclusion.py <dir>     (<dir>/NIST/...: GIAB BEDs as below; <dir>/ucsc/: gap.txt.gz, centromeres.txt.gz, altLocations.txt.gz from
       hgdownload.soe.ucsc.edu/goldenPath/hg38/database/; <dir>/ichorCNA/gc_hg38_500kb.wig from the ichorCNA repository)
"""
import gzip, os, shutil, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "resources/exclusion"; OUT.mkdir(parents=True, exist_ok=True)
BIN, BIN_FRAC, PAD, SITE_PAD = 500_000, 0.5, 0, 50    # no bin padding; 50 bp around site intervals
AUTO = [f"chr{i}" for i in range(1, 23)]
G = Path(sys.argv[1]).expanduser()
SRC = {"segdups": G / "NIST/giab_strat/GRCh38_segdups.bed.gz",
       "tandem_repeats_ge101bp": G / "NIST/giab_strat/GRCh38_AllTandemRepeats_ge101bp_slop5.bed.gz",
       "vdj": G / "NIST/v06_shortread/vdj_loci.GRCh38.bed",
       "ucsc_gap": G / "ucsc/gap.txt.gz"}
UCSC_CEN = G / "ucsc/centromeres.txt.gz"
UCSC_ALT = G / "ucsc/altLocations.txt.gz"
ICHOR_GC = G / "ichorCNA/gc_hg38_500kb.wig"
GC_TOP = 0.01


def read(f):
    d = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"], comment="#",
                    compression="infer")
    d = d[d.chrom.isin(AUTO)].dropna()
    return d.astype({"start": np.int64, "end": np.int64})


def merge(d):
    out = []
    for c, g in d.sort_values(["chrom", "start"]).groupby("chrom", sort=False):
        s, e = g.start.to_numpy(), g.end.to_numpy()
        cs, ce = s[0], e[0]
        for a, b in zip(s[1:], e[1:]):
            if a <= ce:
                ce = max(ce, b)
            else:
                out.append((c, cs, ce)); cs, ce = a, b
        out.append((c, cs, ce))
    return pd.DataFrame(out, columns=["chrom", "start", "end"]).astype({"start": np.int64, "end": np.int64})


def read_ucsc():
    """Telomeres + heterochromatin from the gap table; one interval per chromosome spanning its centromere model."""
    g = pd.read_csv(SRC["ucsc_gap"], sep="\t", header=None, usecols=[1, 2, 3, 7], names=["chrom", "start", "end", "type"])
    g = g[g.chrom.isin(AUTO) & g.type.isin(["telomere", "heterochromatin"])]
    c = pd.read_csv(UCSC_CEN, sep="\t", header=None, usecols=[1, 2, 3], names=["chrom", "start", "end"])
    c = c[c.chrom.isin(AUTO)].groupby("chrom").agg(start=("start", "min"), end=("end", "max")).reset_index()
    d = pd.concat([g, c.assign(type="centromere")]).sort_values(["chrom", "start"])
    d.to_csv(OUT / "ucsc_gap_centromeres_hg38.bed", sep="\t", header=False, index=False)
    return d[["chrom", "start", "end"]].astype({"start": np.int64, "end": np.int64})


R = {k: read(f) for k, f in SRC.items() if k != "ucsc_gap"}
R["ucsc_gap"] = read_ucsc()
for k, f in SRC.items():
    if k != "ucsc_gap":
        shutil.copy(f, OUT / f.name)
allsrc = pd.concat(R.values())
sites = merge(allsrc.assign(start=(allsrc.start - SITE_PAD).clip(lower=0), end=allsrc.end + SITE_PAD))
sites.to_csv(OUT / "sites_hg38.bed.gz", sep="\t", header=False, index=False, compression="gzip")

fai = pd.read_csv(Path(os.environ.get("REF", ROOT / "resources/hg38.fa")).with_suffix(".fa.fai"), sep="\t", header=None, index_col=0)[1]
binsrc = merge(pd.concat([R["segdups"], R["vdj"], R["ucsc_gap"]]))
rows = []
for c in AUTO:
    n = int(fai[c]) // BIN + 1
    cov = np.zeros(n)
    for r in binsrc[binsrc.chrom == c].itertuples():
        for b in range(r.start // BIN, (r.end - 1) // BIN + 1):
            cov[b] += min(r.end, (b + 1) * BIN) - max(r.start, b * BIN)
    hit = cov / BIN >= BIN_FRAC
    pad = hit.copy(); k = PAD // BIN
    for i in np.where(hit)[0]:
        pad[max(0, i - k):i + k + 1] = True
    rows += [(c, b * BIN + 1, (b + 1) * BIN, round(cov[b] / BIN, 3), "covered" if hit[b] else "pad") for b in np.where(pad)[0]]
B = pd.DataFrame(rows, columns=["chrom", "start", "end", "frac_excluded", "reason"])
B.to_csv(OUT / "bins_500kb.tsv", sep="\t", index=False)
alt = pd.read_csv(UCSC_ALT, sep="\t", header=None, usecols=[1, 2, 3], names=["chrom", "start", "end"])
alt = alt[alt.chrom.isin(AUTO)]
ab = sorted({(c, b) for c, s, e in alt.itertuples(index=False) for b in range(s // BIN, (e - 1) // BIN + 1)},
            key=lambda x: (int(x[0][3:]), x[1]))
gc, chrom = [], None
for line in ICHOR_GC.read_text().splitlines():                # fixedStep wig, one value per 500 kb bin
    if line.startswith("fixedStep"):
        kv = dict(x.split("=") for x in line.split()[1:]); chrom = "chr" + kv["chrom"].replace("chr", "")
        pos, step = int(kv["start"]), int(kv["step"])
    else:
        gc.append((chrom, pos, float(line))); pos += step
gc = pd.DataFrame(gc, columns=["chrom", "start", "gc"]); gc = gc[gc.chrom.isin(AUTO) & (gc.gc > 0)]
top = gc[gc.gc >= gc.gc.quantile(1 - GC_TOP)]
pm = pd.DataFrame([(c, b * BIN + 1, (b + 1) * BIN, "ucsc_alt_locus") for c, b in ab]
                  + [(c, s, s + BIN - 1, "gc_top1pct") for c, s in zip(top.chrom, top.start) if (c, (s - 1) // BIN) not in set(ab)],
                  columns=["chrom", "start", "end", "reason"])
pm.assign(o=pm.chrom.str[3:].astype(int)).sort_values(["o", "start"]).drop(columns="o") \
    .to_csv(OUT / "public_bin_mask_500kb.tsv", sep="\t", index=False)
# ---- review: per-source size and 1kGP panel sites removed; every masked bin block with its composition ----
panel = pd.concat([pd.read_csv(ROOT / f"resources/sites_hg38/{c}.sites.tsv.gz", sep="\t", header=None, usecols=[0, 1],
                               names=["chrom", "pos"]) for c in AUTO])


def hits(iv):
    n = 0
    for c, g in panel.groupby("chrom"):
        m = merge(iv[iv.chrom == c]) if (iv.chrom == c).any() else None
        if m is None:
            continue
        s_, e_ = m.start.to_numpy(), m.end.to_numpy()
        p = g.pos.to_numpy() - 1
        i = np.searchsorted(s_, p, side="right") - 1
        n += int(((i >= 0) & (p < e_[np.clip(i, 0, None)])).sum())
    return n


stats = [dict(source=k, intervals=len(v), Mb=round((v.end - v.start).sum() / 1e6, 1), panel_sites_removed=hits(v))
         for k, v in R.items()]
stats.append(dict(source="UNION (sites mask)", intervals=len(sites), Mb=round((sites.end - sites.start).sum() / 1e6, 1),
                  panel_sites_removed=hits(sites)))
S_ = pd.DataFrame(stats); S_["pct_panel"] = (100 * S_.panel_sites_removed / len(panel)).round(1)
S_.to_csv(OUT / "review_sources.tsv", sep="\t", index=False)
print(f"1kGP panel sites (autosomes): {len(panel):,}"); print(S_.to_string(index=False))
B["blk"] = ((B.chrom != B.chrom.shift()) | (B.start != B.end.shift() + 1)).cumsum()
comp = {k: R[k] for k in ["segdups", "vdj", "ucsc_gap"]}


def composition(c, a, b):
    out = []
    for k, v in comp.items():
        x = v[(v.chrom == c) & (v.end > a) & (v.start < b)]
        mb = sum(min(b, r.end) - max(a, r.start) for r in x.itertuples()) / 1e6
        if mb > 0:
            out.append(f"{k} {mb:.2f}")
    return "; ".join(out)


blocks = B.groupby("blk").agg(chrom=("chrom", "first"), start=("start", "min"), end=("end", "max"),
                              bins=("start", "size")).reset_index(drop=True)
blocks["Mb"] = blocks.bins * BIN / 1e6
blocks["contains (Mb, before padding)"] = [composition(r.chrom, r.start, r.end) for r in blocks.itertuples()]
blocks.to_csv(OUT / "review_bin_blocks.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_colwidth", 80)
print(f"\nmasked depth-bin blocks ({len(blocks)}):")
print(blocks.assign(start=(blocks.start / 1e6).round(2), end=(blocks.end / 1e6).round(2)).to_string(index=False))
ORIGIN = {"segdups": "GIAB GRCh38 genome stratifications (NIST)", "tandem_repeats_ge101bp": "GIAB GRCh38 genome stratifications (NIST)",
          "vdj": "GIAB HG002 SV v0.6 (NIST)", "ucsc_gap": "UCSC hg38 database", "ucsc_centromeres": "UCSC hg38 database", "ucsc_alt": "UCSC hg38 database", "ichor_gc": "ichorCNA inst/extdata"}
(OUT / "SOURCES.txt").write_text("".join(f"{k}\t{f.name}\t{ORIGIN[k]}\n"     # file names, not local paths
                                         for k, f in {**SRC, "ucsc_centromeres": UCSC_CEN, "ucsc_alt": UCSC_ALT, "ichor_gc": ICHOR_GC}.items()))
print(f"sites: {len(sites)} merged intervals, {(sites.end - sites.start).sum() / 1e6:.0f} Mb (autosomes)")
print(f"public bin mask (used without a healthy-control mask): {len(pm)} bins = {len(ab)} alt-locus + "
      f"{int((pm.reason == 'gc_top1pct').sum())} top-1% GC (GC >= {top.gc.min():.3f})")
print(f"bins: {int((B.reason == 'covered').sum())} bins >= {BIN_FRAC:.0%} excluded + {int((B.reason == 'pad').sum())} padding bins "
      f"({len(B) * BIN / 1e6:.0f} Mb)")
