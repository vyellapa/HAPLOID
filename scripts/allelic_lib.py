"""Shared pieces for the unphased allelic ("both-allele rate") vs ichorCNA comparison.

Prediction per site with n reads and panel AF p: P(both alleles seen) = 2pq * (1 - b^n - (1-b)^n),
where b = ((1-TF)*1 + TF*minor) / ((1-TF)*2 + TF*CN). Minor copies are not known from ichorCNA,
so minor = 0 for CN <= 1, else floor(CN/2). Subclonal segments use TF * subclone fraction.
"""
import os, subprocess
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
# env overrides for running elsewhere (HPC): HAPLOID_RES = dir with hg38-blacklist.v2.bed + cytoband_b38.bed;
# ICHOR_PARAMS / ICHOR_SEG = ichorCNA file paths, "{s}" is replaced by the sample name
# HAPLOID_CONTROLS = dir with control_bin_mask_500kb.tsv / pipeline_bias_500kb.tsv (default resources/)
RES = Path(os.environ.get("HAPLOID_RES", ROOT / "resources"))
ICHOR_PARAMS = os.environ.get("ICHOR_PARAMS", str(ROOT / "ichor/{s}_Y.params.txt"))
ICHOR_SEG = os.environ.get("ICHOR_SEG", str(ROOT / "ichor/{s}_Y.seg.txt"))
CEN_PAD = 2_000_000   # exclude +/- 2 Mb around each centromere (acen bands)
EXCL = ROOT / "resources/exclusion"   # build_exclusion.py: GIAB segdups / tandem repeats >= 101 bp / VDJ + UCSC gap exclusions, +/- 50 bp
MIN_SITES = 300       # bins with fewer informative sites are dropped

# reference palette (dataviz skill), light mode
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
CALL_COL = {"HETD": "#2a78d6", "HOMD": "#2a78d6", "NEUT": "#8f8e88", "GAIN": "#eb6834", "AMP": "#e34948",
            "HLAMP": "#e34948"}


def ichor_kv(s):
    kv = {}
    for line in Path(ICHOR_PARAMS.replace("{s}", s)).read_text().splitlines():
        if ":\t" in line:
            k, v = line.split(":\t", 1); kv[k.strip()] = v.strip()
    return kv


def params(s):
    kv = ichor_kv(s)
    sf = kv["Subclone Fraction"]                          # ichorCNA writes "NA" when there is no subclone
    return float(kv["Tumor Fraction"]), (np.nan if sf == "NA" else float(sf))


def ichor_ploidy(s):
    return float(ichor_kv(s)["Ploidy"])


def segments(s):
    return pd.read_csv(ICHOR_SEG.replace("{s}", s), sep="\t")


def bq(f, fmt):
    p = subprocess.Popen(["bcftools", "query", "-f", fmt, str(f)], stdout=subprocess.PIPE)
    return pd.read_csv(p.stdout, sep="\t", header=None, na_values=".")


CONTROLS = Path(os.environ.get("HAPLOID_CONTROLS", ROOT / "resources"))   # control_bin_mask / pipeline_bias (control_bin_mask.py)


def bin_mask():
    """500 kb bins that misbehave in healthy controls (control_bin_mask.py). Without that file, the public mask
    (exclusion/public_bin_mask_500kb.tsv, build_exclusion.py: hg38 alternate-haplotype loci + the 1% most GC-rich bins),
    unless HAPLOID_PUBLIC_MASK=0. Empty if neither applies."""
    f = CONTROLS / "control_bin_mask_500kb.tsv"
    if not f.exists() and os.environ.get("HAPLOID_PUBLIC_MASK", "1") != "0":
        f = EXCL / "public_bin_mask_500kb.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame(columns=["chrom", "start", "end"])


def exclusion_bins():
    """500 kb bins (chrom, start) excluded for depth: >= 50% segdup/VDJ/UCSC-gap-excluded, no padding (build_exclusion.py)."""
    f = EXCL / "bins_500kb.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame(columns=["chrom", "start", "end"])


_SITE_EXCL, _MASK = None, {}


def mask_intervals(c):
    """Merged, sorted (n, 2) intervals masked for allele sites on chromosome c: ENCODE blacklist, heterochromatin
    cytobands, centromeres +/- CEN_PAD, control bin mask, and GIAB segdups / tandem repeats >= 101 bp / VDJ + UCSC gap exclusions, +/- 50 bp
    (mismapping -> false hets or false allelic loss)."""
    global _SITE_EXCL
    if c in _MASK:
        return _MASK[c]
    if _SITE_EXCL is None:
        f = EXCL / "sites_hg38.bed.gz"
        _SITE_EXCL = (pd.read_csv(f, sep="\t", header=None, names=[0, 1, 2]) if f.exists()
                      else pd.DataFrame(columns=[0, 1, 2]))
    bl = pd.read_csv(RES / "hg38-blacklist.v2.bed", sep="\t", header=None, usecols=[0, 1, 2])
    cb = pd.read_csv(RES / "cytoband_b38.bed", sep="\t", header=None)
    het = cb[cb[4].isin(["gvar", "stalk"])][[0, 1, 2]]
    cen = cb[cb[4] == "acen"].groupby(0).agg({1: "min", 2: "max"}).reset_index()
    cen[1], cen[2] = (cen[1] - CEN_PAD).clip(lower=0), cen[2] + CEN_PAD
    bm = bin_mask()[["chrom", "start", "end"]].set_axis([0, 1, 2], axis=1)
    iv = pd.concat([bl, het, cen, bm, _SITE_EXCL])
    iv = iv[iv[0] == f"chr{c}"][[1, 2]].astype(np.int64).sort_values(1).to_numpy()
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    _MASK[c] = np.array(merged, dtype=np.int64).reshape(-1, 2)
    return _MASK[c]


def in_mask(c, pos):
    """Boolean array: positions on chromosome c (int or 'chrN') that fall in a masked interval."""
    iv = mask_intervals(int(str(c).replace("chr", "")))
    pos = np.asarray(pos)
    if not len(iv):
        return np.zeros(len(pos), bool)
    i = np.searchsorted(iv[:, 0], pos, side="right") - 1
    return (i >= 0) & (pos <= iv[np.clip(i, 0, None), 1])


def site_table(s, c):
    ad = bq(ROOT / f"work/gl/{s}.chr{c}.gl.bcf", "%POS\t[%AD]\n").dropna(); ad.columns = ["pos", "ad"]
    af = bq(ROOT / f"resources/sites_hg38/chr{c}.sites.vcf.gz", "%POS\t%INFO/AF\n"); af.columns = ["pos", "af"]
    d = ad.merge(af, on="pos")
    d = d[~in_mask(c, d.pos.to_numpy())]
    d[["r", "a"]] = d["ad"].str.split(",", expand=True).astype(int)
    d["n"] = d.r + d.a
    d = d[(d.n >= 2) & (d.n <= 8) & d.af.between(0.05, 0.95)].copy()
    d["obs"] = ((d.r > 0) & (d.a > 0)).astype(float)
    d["exp"] = 2 * d.af * (1 - d.af) * (1 - 2 * 0.5 ** d.n)
    return d


def annotate_cn(d, seg, tf, sf):
    d["cn"], d["minor"], d["tf"] = 2.0, 1.0, tf
    for _, g in seg.iterrows():
        m = d.pos.between(g.start, g.end)
        cn = float(g.Corrected_Copy_Number)
        d.loc[m, "cn"] = cn
        d.loc[m, "minor"] = 0.0 if cn <= 1 else np.floor(cn / 2)
        d.loc[m, "tf"] = tf * sf if str(g["subclone.status"]).upper() == "TRUE" else tf
    b = ((1 - d.tf) + d.tf * d.minor) / ((1 - d.tf) * 2 + d.tf * d.cn)
    d["pred"] = 2 * d.af * (1 - d.af) * (1 - b ** d.n - (1 - b) ** d.n)
    return d


def binned(d, bin_mb):
    d["bin"] = (d.pos // (bin_mb * 1e6)).astype(int)
    g = d.groupby("bin").agg(mid=("pos", "median"), sites=("obs", "size"), obs=("obs", "sum"),
                             exp=("exp", "sum"), pred=("pred", "sum"), cn=("cn", "median"))
    g = g[g.sites >= MIN_SITES].copy()
    g["ratio"], g["pred_ratio"] = g.obs / g.exp, g.pred / g.exp
    return g


def style(ax):
    ax.set_facecolor(SURF)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=2)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)


# figure style for publication plots: white background, black axes, grey data, saturated lines
# (red / blue / orange / teal pass the colour-vision-deficiency check for every pair, so all four can share one panel)
BW_DOTS, BW_SEP = "#8c8c8c", "#c8c8c8"
LINE = {"red": "#d7191c", "blue": "#1f5fbf", "orange": "#e69f00", "teal": "#00a6a6"}


def style_bw(ax):
    ax.set_facecolor("white"); ax.grid(False)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color("black"); sp.set_linewidth(0.8)
    ax.tick_params(colors="black", labelsize=8, length=3, width=0.8)


def segment_types(seg):
    """Per segment: the call type of a passing call (LOH / gain / loss / imbalance against the most common clonal
    total copy number), "" for segments that are not passing calls."""
    base = seg[seg.sub_to.isna()].groupby(["nA", "nB"]).mb.sum().idxmax(); bt = sum(base)
    def typ(r):
        a, b = r.nA, r.nB
        if isinstance(r.sub_to, str) and r.sub_to:
            st_ = r.sub_to.split("@")[0].split(":"); a, b = int(st_[0]), int(st_[1])
        if b == 0 and a + b >= 1:
            return "LOH"
        return "gain" if a + b > bt else ("loss" if a + b < bt else "imbalance")
    return [typ(r) if str(r.call_pass) == "True" else "" for r in seg.itertuples()]


def call_events(seg):
    """Passing calls merged into events: consecutive passing segments on a chromosome with the same type."""
    x = seg.copy()
    x["ok"] = x.call_pass.astype(str) == "True"
    if not x.ok.any():
        return pd.DataFrame(columns=["chrom", "start", "end", "type", "mb", "focal"])
    x["type"] = segment_types(x)
    x["_c"] = x.chrom.str[3:].astype(int); x = x.sort_values(["_c", "start"])
    ev, prev_ok = [], False
    for r in x.itertuples():
        if not r.ok:                     # a non-call segment ends the current event
            prev_ok = False; continue
        if prev_ok and ev[-1]["chrom"] == r.chrom and ev[-1]["type"] == r.type:
            ev[-1]["end"] = r.end; ev[-1]["mb"] += r.mb
        else:
            ev.append(dict(chrom=r.chrom, start=r.start, end=r.end, type=r.type, mb=r.mb,
                           focal=str(getattr(r, "focal", False)) == "True"))
        prev_ok = True
    return pd.DataFrame(ev)
