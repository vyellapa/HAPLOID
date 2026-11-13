"""Final results for a run (run_bam.sh): the samples' rows of work/cohort_summary.tsv / cohort_events.tsv, their plots,
and a PDF (A4 landscape: a summary table, then one page per sample with work/plots/cna_<sample>.png).
ichorCNA_source says whether ichorCNA results were supplied or computed here (work/inputs/<s>_500k/ichor_source.txt).
Usage: final_report.py <samples...>     env: HAPLOID_FINAL (output folder, default work/final)
Output: summary.tsv, events.tsv, plots/cna_<s>.png + cna_compare_<s>.png, report.pdf
"""
import os, shutil, sys, textwrap
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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


cols = ["Sample", "Status", "TF (95%)", "Ploidy", "ctDNA", "ichorCNA\nTF / ploidy", "ichorCNA\nresults", "QC", "Events"]
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
                 + f" · ichorCNA {num(r.ichor_TF)} / {num(r.ichor_ploidy)}"
                 + (f" ({r.ichorCNA_source})" if r.ichorCNA_source else "")
                 + (f" · {r.tf_basis}" if isinstance(r.tf_basis, str) and "ichorCNA" in r.tf_basis else ""), fontsize=8.5)
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.87]); ax.axis("off")
        ax.imshow(plt.imread(p), interpolation="lanczos"); ax.set_anchor("N")
        pdf.savefig(fig, dpi=200); plt.close(fig)
print(f"{len(S)} samples ({int((S.tumor_detected.astype(str) == 'True').sum())} tumor detected) -> {OUT}/summary.tsv, "
      f"events.tsv, plots/, report.pdf")
