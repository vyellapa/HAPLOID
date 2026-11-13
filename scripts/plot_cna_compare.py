"""Genome-wide CNA comparison per sample: the data, then each method's copy-number profile on the same axis.
  depth   : 500 kb logR bins (work/cnv/<s>.bins.tsv, the depth the joint fit uses); ichorCNA segment logR (orange)
            (red) and the joint fit's (cbs_union) predicted segment logR (blue)
  alleles : phased BAF per bin (faint; haplotype orientation flips at phase switches, so dots sit at both b and 1 - b),
            segment allelic imbalance 0.5 +/- theta_hat from the allele HMM (orange) and the cbs_union prediction (blue);
            y 0.25-0.75, widened to 0-1 when a predicted allele fraction falls outside
  CN      : ichorCNA (total copy number; it has no allele information), then joint fit in ascat / arms / cbs_union
            mode: total CN (black) and minor-allele CN (red); subclonal states drawn at their effective copy number.
            Background shading = passing calls by copy change (gain orange / loss blue); copy-neutral LOH teal. LOH from a
            deletion is shaded as a loss - the minor-allele line at 0 shows the LOH. ichorCNA shading = its own calls.
Two figures per sample:
  work/plots/cna_compare_<sample>.png : all of the above
  work/plots/cna_<sample>.png         : depth, alleles (taller) and the cbs_union fit, plus its refits when present
                                        (cbs_union_ichorploidy / cbs_union_user); sized to fill an A4 landscape page
Usage: plot_cna_compare.py <samples...>     run by run_sample.sh and refit_downstream.sh after the joint fit
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from allelic_lib import RES, ICHOR_PARAMS, ICHOR_SEG, LINE, BW_DOTS, BW_SEP, style_bw

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = sys.argv[1:] or sys.exit("usage: plot_cna_compare.py <samples...>")
MODES = ["ascat", "arms", "cbs_union", "cbs_union_ichorploidy", "cbs_union_user"]   # refit rows only when present
COMPACT_MODES = ["cbs_union", "cbs_union_ichorploidy", "cbs_union_user"]
LABEL = {"ascat": "ASPCF segments", "arms": "chromosome arms", "cbs_union": "cbs_union",    # shown in plot titles
         "cbs_union_ichorploidy": "cbs_union, ploidy held at the depth HMM's", "cbs_union_user": "cbs_union, your values"}
CALL = {"gain": LINE["orange"], "loss": LINE["blue"], "CN-LOH": LINE["teal"]}
ICHOR, JOINT, HMM, MINOR = LINE["red"], LINE["blue"], LINE["orange"], LINE["red"]
FAI = RES / "hg38.fa.fai"
if not FAI.exists() and os.environ.get("REF"):            # cluster runs: the reference index sits next to $REF
    FAI = Path(os.environ["REF"] + ".fai")
fai = pd.read_csv(FAI, sep="\t", header=None, index_col=0)[1]
CHR = [f"chr{i}" for i in range(1, 23)]
OFF, acc = {}, 0
for c in CHR:
    OFF[c] = acc; acc += int(fai[c])
GLEN = acc
X = lambda c, p: np.asarray([OFF[x] for x in c]) + np.asarray(p)


def ichor_files(s):
    """ICHOR_SEG / ICHOR_PARAMS (set by run_sample.sh), else <s>_500k/ or ichor/ in the project, else $BAM_DIR/<s>_500k/."""
    cand = [(Path(ICHOR_SEG.replace("{s}", s)), Path(ICHOR_PARAMS.replace("{s}", s)))]
    cand += [(d / f"{s}_Y.seg.txt", d / f"{s}_Y.params.txt")
             for d in (ROOT / f"{s}_500k", ROOT / "ichor", Path(os.environ.get("BAM_DIR", ROOT)) / f"{s}_500k")]
    return next(((g, p) for g, p in cand if g.exists() and p.exists()), (None, None))


def ichor(s):
    g, p = ichor_files(s)
    if g is None:                             # ichorCNA results not copied here: plot without them
        return None, float("nan"), float("nan")
    seg = pd.read_csv(g, sep="\t")
    seg = seg[seg.chrom.astype(str).str.replace("chr", "").isin([str(i) for i in range(1, 23)])].copy()
    seg["chrom"] = "chr" + seg.chrom.astype(str).str.replace("chr", "")
    kv = dict(l.split(":\t", 1) for l in p.read_text().splitlines() if ":\t" in l)
    return seg, float(kv["Tumor Fraction"]), float(kv["Ploidy"])


def effective(seg):
    """Effective total and minor copy number per segment (subclonal mixes weighted by their cell fraction)."""
    tot, mnr = [], []
    for r in seg.itertuples():
        a, b = r.nA, r.nB
        if isinstance(r.sub_to, str) and r.sub_to:
            ab, f = r.sub_to.split("@"); a2, b2 = map(int, ab.split(":")); f = float(f)
            a, b = a + f * (a2 - a), b + f * (b2 - b)
        tot.append(a + b); mnr.append(min(a, b))
    return np.array(tot), np.array(mnr)


def call_type(r, tot, mnr, base_tot):
    if str(r.call_pass) != "True":
        return None
    if tot > base_tot + 0.05:
        return "gain"
    if tot < base_tot - 0.05:
        return "loss"
    return "CN-LOH" if mnr < 0.05 else None


def frame(ax, ylab):
    style_bw(ax)
    for c in CHR[1:]:
        ax.axvline(OFF[c], color=BW_SEP, lw=0.5, zorder=0)
    ax.set_xlim(0, GLEN); ax.set_ylabel(ylab, fontsize=8.5, color="black")


def steps(ax, seg, y, **kw):
    for r, v in zip(seg.itertuples(), y):
        ax.plot([OFF[r.chrom] + r.start, OFF[r.chrom] + r.end], [v, v], solid_capstyle="butt", **kw)


def draw(s, bins, iseg, itf, ipl, J, compact):
    if compact:                               # no ichorCNA row; data panels about twice a copy-number row
        J = {m: v for m, v in J.items() if m in COMPACT_MODES}
        n = 2 + len(J); ratios = [2.4, 2.0] + [1.2] * len(J)
        figsize = (17, 1.4 * sum(ratios) + 1.1)
    else:
        n = 3 + len(J); ratios = [1.3, 1.1] + [1] * (n - 2)
        figsize = (17, 1.55 * n + 1.6)
    fig, ax = plt.subplots(n, 1, figsize=figsize, facecolor="white", sharex=True,
                           gridspec_kw=dict(height_ratios=ratios, hspace=0.42))
    # ---- depth ----
    a = ax[0]; frame(a, "logR")
    a.scatter(X(bins.chrom, bins.pos), bins.logR, s=4, color=BW_DOTS, alpha=0.4, lw=0, rasterized=True)
    if iseg is not None:
        steps(a, iseg, iseg["seg.median.logR"], color=ICHOR, lw=2)
    if "cbs_union" in J:
        seg, sm = J["cbs_union"]; rho, psi = sm.joint_TF, sm.joint_ploidy
        tot, mnr = effective(seg)
        pred = np.log2((2 * (1 - rho) + rho * tot) / (2 * (1 - rho) + rho * psi))
        steps(a, seg, pred, color=JOINT, lw=1.6)
    a.set_ylim(-1.2, 1.2); a.axhline(0, color="black", lw=0.6, ls=(0, (2, 2)), zorder=0)
    a.set_title("depth: 500 kb bins (grey) · depth HMM segments (red) · joint fit cbs_union prediction (blue)",
                fontsize=9, color="black", loc="left")
    # ---- alleles ----
    a = ax[1]; frame(a, "BAF")
    ok = bins.baf.notna()
    a.scatter(X(bins.chrom[ok], bins.pos[ok]), bins.baf[ok], s=4, color=BW_DOTS, alpha=0.4, lw=0, rasterized=True)
    hb = 0.25
    if "cbs_union" in J:
        seg, sm = J["cbs_union"]; rho = sm.joint_TF
        tot, mnr = effective(seg)
        b = (1 - rho + rho * mnr) / (2 * (1 - rho) + rho * tot)
        for sgn in (1, -1):
            steps(a, seg, 0.5 + sgn * seg.theta_hat, color=HMM, lw=1.8)
        steps(a, seg, b, color=JOINT, lw=1.4); steps(a, seg, 1 - b, color=JOINT, lw=1.4)
        hb = 0.25 if b.min() >= 0.26 and (0.5 - seg.theta_hat).min() >= 0.26 else 0.5
    a.set_ylim(0.5 - hb, 0.5 + hb); a.set_yticks([0.5 - hb, 0.5, 0.5 + hb])
    a.axhline(0.5, color="black", lw=0.6, ls=(0, (2, 2)), zorder=0)
    a.set_title("alleles: phased BAF per bin (grey) · segment imbalance from the allele HMM (orange) · "
                "cbs_union prediction (blue)", fontsize=9, color="black", loc="left")
    # ---- ichorCNA CN ----
    if not compact:
        a = ax[2]; frame(a, "copy number")
        for r in (iseg.itertuples() if iseg is not None else []):
            t = {"GAIN": "gain", "AMP": "gain", "HLAMP": "gain", "HETD": "loss", "HOMD": "loss"}.get(r.Corrected_Call)
            if t:
                a.axvspan(OFF[r.chrom] + r.start, OFF[r.chrom] + r.end, color=CALL[t], alpha=0.18, lw=0)
        if iseg is not None:
            steps(a, iseg, iseg.Corrected_Copy_Number.clip(upper=6), color="black", lw=1.8)
        a.set_ylim(-0.3, 6.3); a.set_yticks([0, 2, 4, 6])
        a.set_title(f"depth HMM · TF {itf:.2f} · ploidy {ipl:.2f} (total copy number only)" if iseg is not None else
                    f"depth HMM: no {s}_500k/{s}_Y.seg.txt found next to the project (copy the depth results to show it)",
                    fontsize=9, color="black", loc="left")
    # ---- joint fit modes ----
    for a, (m, (seg, sm)) in zip(ax[2 if compact else 3:], J.items()):
        frame(a, "copy number")
        tot, mnr = effective(seg)
        base = seg[seg.sub_to.isna()].groupby(["nA", "nB"]).mb.sum().idxmax(); bt = sum(base)
        for r, t_, m_ in zip(seg.itertuples(), tot, mnr):
            t = call_type(r, t_, m_, bt)
            if t:
                a.axvspan(OFF[r.chrom] + r.start, OFF[r.chrom] + r.end, color=CALL[t], alpha=0.18, lw=0)
        steps(a, seg, np.clip(tot, 0, 6) + 0.06, color="black", lw=1.8)
        steps(a, seg, np.clip(mnr, 0, 6) - 0.06, color=MINOR, lw=1.6)
        a.set_ylim(-0.3, 6.3); a.set_yticks([0, 2, 4, 6])
        det = str(sm.tumor_detected) == "True"
        ct = f" · ctDNA {sm.ctDNA_fraction:.2f}" if "ctDNA_fraction" in sm and det else ""
        amb = " · AMBIGUOUS ploidy" if str(sm.ambiguous) == "True" else ""
        con = f" · constrained: {sm.constraints} (dLL {sm.dLL_constraint:.0f})" if isinstance(sm.get("constraints"), str) and sm.constraints else ""
        head = (f"joint fit [{LABEL.get(m, m)}]{con} · TF {sm.joint_TF:.2f} · ploidy {sm.joint_ploidy:.2f}{ct} · score {sm.tumor_llr:.0f}"
                f" · other-ploidy dLL {sm.dLL_other_ploidy:.0f}{amb} · {int(sm.calls_pass)} passing calls"
                if det else f"joint fit [{LABEL.get(m, m)}] · tumor NOT detected (score {sm.tumor_llr:.1f})")
        a.set_title(head, fontsize=9, color="black", loc="left")
    ax[-1].set_xticks([OFF[c] + fai[c] / 2 for c in CHR]); ax[-1].set_xticklabels([c[3:] for c in CHR], fontsize=8)
    handles = [Line2D([], [], color="black", lw=1.8, label="total CN"), Line2D([], [], color=MINOR, lw=1.6, label="minor-allele CN")]
    handles += [Patch(color=v, alpha=0.35, label="copy-neutral LOH" if k == "CN-LOH" else k) for k, v in CALL.items()]
    fig.legend(handles=handles, loc="upper right", ncol=5, frameon=False, fontsize=8.5, labelcolor="black",
               bbox_to_anchor=(0.99, 0.995))
    fig.suptitle(f"{s}: " + ("depth, alleles and joint fit" if compact else "copy-number profiles by method")
                 + " (shading = passing calls)", x=0.01, y=0.995, ha="left", fontsize=12, fontweight="bold", color="black")
    fig.subplots_adjust(top=1 - 0.75 / fig.get_figheight())
    out = ROOT / f"work/plots/{'cna' if compact else 'cna_compare'}_{s}.png"; out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150, facecolor="white", bbox_inches="tight"); plt.close(fig)
    print(f"-> {out}")


for s in SAMPLES:
    bins = pd.read_csv(ROOT / f"work/cnv/{s}.bins.tsv", sep="\t")
    iseg, itf, ipl = ichor(s)
    J = {m: (pd.read_csv(ROOT / f"work/joint_modes/{m}/{s}.joint_segments.tsv", sep="\t"),
             pd.read_csv(ROOT / f"work/joint_modes/{m}/{s}.joint_summary.tsv", sep="\t").iloc[0])
         for m in MODES if (ROOT / f"work/joint_modes/{m}/{s}.joint_segments.tsv").exists()}
    for compact in (False, True):
        draw(s, bins, iseg, itf, ipl, J, compact)
