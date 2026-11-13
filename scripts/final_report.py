"""Final results for a run (run_bam.sh): the samples' rows of work/cohort_summary.tsv / cohort_events.tsv, their plots,
and a PDF (A4 landscape: a summary table, then one page per sample with work/plots/cna_<sample>.png).
ichor_TF / ichor_ploidy = the depth-only HMM estimate; ichorCNA_source = whether those results were supplied or computed
here (work/inputs/<s>_500k/ichor_source.txt).
segments.tsv: every segment of the fit each sample is reported from (cohort_summary segments_from), with its depth,
allele balance (as the major-allele fraction, as in the plots' BAF panel), fitted copy numbers and call.
Usage: final_report.py <samples...>     env: HAPLOID_FINAL (output folder, default work/final)
Output: summary.tsv, events.tsv, segments.tsv, plots/cna_<s>.png + cna_compare_<s>.png, report.pdf
"""
import os, shutil, sys, textwrap
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from allelic_lib import segment_types, effective

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = sys.argv[1:] or sys.exit("usage: final_report.py <samples...>")
OUT = Path(os.environ.get("HAPLOID_FINAL", ROOT / "work/final")); (OUT / "plots").mkdir(parents=True, exist_ok=True)
A4 = (11.69, 8.27)

S = pd.read_csv(ROOT / "work/cohort_summary.tsv", sep="\t", dtype={"sample": str})
S = S[S["sample"].isin(SAMPLES)].copy()
src = {s: (ROOT / f"work/inputs/{s}_500k/ichor_source.txt") for s in SAMPLES}
S.insert(S.columns.get_loc("ichor_ploidy") + 1, "ichorCNA_source",
         [src[s].read_text().strip() if src[s].exists() else "" for s in S["sample"]])
S.to_csv(OUT / "summary.tsv", sep="\t", index=False)
E = pd.read_csv(ROOT / "work/cohort_events.tsv", sep="\t", dtype={"sample": str})
E = E[E["sample"].isin(SAMPLES)]
E.to_csv(OUT / "events.tsv", sep="\t", index=False)
SEG_COLS = {"chrom": "chrom", "start": "start", "end": "end", "mb": "size_mb", "bins": "bins", "r_obs": "logR",
            "hets_mb": "het_sites_per_mb", "nA": "major_cn", "nB": "minor_cn",
            "sub_to": "subclonal_state", "focal": "focal", "filter_reason": "filter_reason"}


def baf_median(g, bins):
    """Median over the segment's 500 kb bins of the phased BAF mirrored to >= 0.5 (the grey dots, folded onto the upper
    band: phase switches put an imbalanced segment's bins on both sides of 0.5, so the raw median drifts to 0.5)."""
    out = []
    for c, s0, e0 in zip(g.chrom, g.start, g.end):
        x = bins[(bins.chrom == c) & (bins.start >= s0) & (bins.start < e0)].baf.dropna()
        out.append(0.5 + float(np.median(np.abs(x - 0.5))) if len(x) else np.nan)
    return out


segs = []
for _, r in S.iterrows():
    fit = r.get("segments_from") if isinstance(r.get("segments_from"), str) else "arms"
    f = ROOT / f"work/joint_modes/{fit}/{r['sample']}.joint_segments.tsv"
    if not f.exists():
        continue
    seg = pd.read_csv(f, sep="\t")
    det = str(r.tumor_detected) == "True"
    g = seg[list(SEG_COLS)].rename(columns=SEG_COLS)
    g.insert(0, "fit", fit); g.insert(0, "sample", r["sample"])
    g.insert(g.columns.get_loc("minor_cn") + 1, "total_cn", g.major_cn + g.minor_cn)
    g.insert(g.columns.get_loc("subclonal_state") + 1, "call", segment_types(seg) if det else [""] * len(g))
    bf = ROOT / f"work/cnv/{r['sample']}.bins.tsv"
    g.insert(g.columns.get_loc("logR") + 1, "baf_median",
             baf_median(g, pd.read_csv(bf, sep="\t")) if bf.exists() else np.nan)
    g.insert(g.columns.get_loc("baf_median") + 1, "baf_hmm", 0.5 + seg.theta_hat)
    rho = float(r.TF) if det else 0.0                      # the fit's tumor fraction; 0 = all normal
    tot, mnr = effective(seg)
    g.insert(g.columns.get_loc("baf_hmm") + 1, "baf_fit", 1 - (1 - rho + rho * mnr) / (2 * (1 - rho) + rho * tot))
    segs.append(g.round({"size_mb": 2, "logR": 4, "baf_median": 3, "baf_hmm": 3, "baf_fit": 3, "het_sites_per_mb": 1}))
(pd.concat(segs) if segs else pd.DataFrame(columns=["sample", "fit", *SEG_COLS.values()])).to_csv(
    OUT / "segments.tsv", sep="\t", index=False)
for s in SAMPLES:
    for p in (f"cna_{s}.png", f"cna_compare_{s}.png"):
        if (ROOT / "work/plots" / p).exists():
            shutil.copy(ROOT / "work/plots" / p, OUT / "plots" / p)


def num(v, d=2):
    try:
        return f"{float(v):.{d}f}"
    except (TypeError, ValueError):
        return ""


def row(r):
    det = str(r.tumor_detected) == "True"
    flags = [f for f in (r.fit_flag if isinstance(r.fit_flag, str) else "",
                         "allele QC low" if r.allele_qc == "low" else "") if f]
    ev = E[E["sample"] == r["sample"]]
    evs = ", ".join(f"{c.replace('chr', '')}:{t}" for c, t in zip(ev.chrom, ev.type))
    return [r["sample"], "detected" if det else "not detected",
            f"{num(r.TF)} ({num(r.TF_lo)}-{num(r.TF_hi)})" if det else "",
            num(r.ploidy) if det else "", num(r.ctDNA_fraction) if det else "",
            f"{num(r.ichor_TF)} / {num(r.ichor_ploidy)}", r.ichorCNA_source.split(" ")[0],
            "\n".join(flags), textwrap.fill(evs, 48)]


cols = ["Sample", "Status", "TF (95%)", "Ploidy", "ctDNA", "Depth HMM\nTF / ploidy", "Depth\nresults", "QC", "Events"]
WIDTHS = [0.25, 0.075, 0.09, 0.045, 0.045, 0.075, 0.07, 0.09, 0.26]
PER_PAGE_LINES = 34
rows = [row(r) for _, r in S.iterrows()]
with PdfPages(OUT / "report.pdf") as pdf:
    pages, cur, n = [], [], 0                             # rows of several lines (wrapped events) fill a page sooner
    for r in rows:
        k = max(c.count("\n") + 1 for c in r)
        if cur and n + k > PER_PAGE_LINES:
            pages.append(cur); cur, n = [], 0
        cur.append(r); n += k
    pages.append(cur)
    for part in pages:
        fig = plt.figure(figsize=A4); ax = fig.add_axes([0.03, 0.05, 0.94, 0.82]); ax.axis("off")
        fig.text(0.03, 0.94, "HAPLOID results", fontsize=16, fontweight="bold")
        fig.text(0.03, 0.905, "TF = tumor cell fraction with its 95% segment-bootstrap interval; ctDNA = tumor share of the "
                 "cfDNA. QC: uniform_imbalance = likely a second person's DNA or noisy allele counts, not tumor; "
                 "allele QC low = too few allele reads, so 'no LOH' is not a negative result.", fontsize=7.5, wrap=True)
        if part:
            t = ax.table(cellText=part, colLabels=cols, loc="upper left", cellLoc="left", colWidths=WIDTHS)
            t.auto_set_font_size(False); t.set_fontsize(7)
            lines = [2] + [max(c.count("\n") + 1 for c in r) for r in part]
            for (rr, _), c in t.get_celld().items():
                c.set_edgecolor("#cccccc"); c.set_height(0.024 * lines[rr] + 0.008); c.PAD = 0.03
                if rr == 0:
                    c.set_text_props(fontweight="bold")
        pdf.savefig(fig); plt.close(fig)
    for s in SAMPLES:
        p = ROOT / f"work/plots/cna_{s}.png"
        if not p.exists():
            continue
        r = S[S["sample"] == s].iloc[0]
        fig = plt.figure(figsize=A4)
        fig.text(0.03, 0.95, s, fontsize=12, fontweight="bold")
        det = str(r.tumor_detected) == "True"
        fig.text(0.03, 0.92, (f"detected · TF {num(r.TF)} ({num(r.TF_lo)}-{num(r.TF_hi)}) · ploidy {num(r.ploidy)} · "
                              f"ctDNA {num(r.ctDNA_fraction)}" if det else "not detected")
                 + f" · depth HMM {num(r.ichor_TF)} / {num(r.ichor_ploidy)}"
                 + (f" ({r.ichorCNA_source})" if r.ichorCNA_source else "")
                 + (f" · {r.tf_basis}" if isinstance(r.tf_basis, str) and "depth HMM" in r.tf_basis else ""), fontsize=8.5)
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.87]); ax.axis("off")
        ax.imshow(plt.imread(p), interpolation="lanczos"); ax.set_anchor("N")
        pdf.savefig(fig, dpi=200); plt.close(fig)
print(f"{len(S)} samples ({int((S.tumor_detected.astype(str) == 'True').sum())} tumor detected) -> {OUT}/summary.tsv, "
      f"events.tsv, segments.tsv, plots/, report.pdf")
