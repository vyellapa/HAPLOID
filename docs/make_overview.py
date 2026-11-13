"""README / wiki overview figure (docs/haploid-overview.png); usage, from the repository root: python3 docs/make_overview.py
HAPLOID overview figure from synthetic data (no sample data): reads -> phasing -> depth + haplotype BAF -> joint fit
-> results."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

OUT = sys.argv[1] if len(sys.argv) > 1 else "docs/haploid-overview.png"
rng = np.random.default_rng(7)
INK, MUTED, GRID = "#1f2933", "#6b7785", "#d9dee4"
BLUE, ORANGE, TEAL, RED, PURPLE = "#2f6fb0", "#e08a1e", "#1f9e8f", "#c8423b", "#7a5aa6"
H1, H2 = "#2f6fb0", "#e08a1e"           # haplotype colours
PANEL_BG = "#f6f8fa"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK})

fig = plt.figure(figsize=(17, 6.6), facecolor="white")
# panel frames (figure coords): a, b, c, d, e
FR = {"a": (0.012, 0.10, 0.155, 0.80), "b": (0.185, 0.10, 0.165, 0.80), "c": (0.368, 0.10, 0.255, 0.80),
      "d": (0.641, 0.10, 0.205, 0.80), "e": (0.864, 0.10, 0.124, 0.80)}
TITLE = {"a": "Low-pass cfDNA", "b": "Phase and pool", "c": "Depth and haplotype BAF", "d": "Joint fit",
         "e": "Results"}
SUB = {"a": "~1-2x: a heterozygous SNP\nhas 0-2 reads", "b": "impute + phase against the\n1000 Genomes panel; sum reads\nper haplotype",
       "c": "500 kb bins, segmented;\nCN-LOH is visible only in BAF", "d": "tumor fraction x ploidy grid;\nallele-specific copy number",
       "e": "per sample"}
for k, (x, y, w, h) in FR.items():
    fig.patches.append(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                      transform=fig.transFigure, fc=PANEL_BG, ec=GRID, lw=1.0, zorder=-10))
    fig.text(x + 0.010, y + h - 0.045, k, fontsize=15, weight="bold", color=INK)
    fig.text(x + 0.030, y + h - 0.043, TITLE[k], fontsize=11.5, weight="bold", color=INK, va="baseline")
    fig.text(x + 0.010, y + h - 0.075, SUB[k], fontsize=8.6, color=MUTED, va="top", linespacing=1.35)
for a, b in zip("abcd", "bcde"):
    x0 = FR[a][0] + FR[a][2]; x1 = FR[b][0]
    fig.patches.append(FancyArrowPatch((x0 + 0.001, 0.5), (x1 - 0.001, 0.5), transform=fig.transFigure,
                                       arrowstyle="-|>,head_length=6,head_width=4", color=MUTED, lw=1.6, zorder=5))
fig.text(0.5, 0.955, "HAPLOID: haplotype-aware LOH and allele-specific copy number from low-pass cfDNA",
         ha="center", fontsize=13.5, weight="bold", color=INK)
fig.text(0.5, 0.035, "Synthetic illustration", ha="center", fontsize=7.5, color=MUTED, style="italic")


def axes_in(k, rect):
    x, y, w, h = FR[k]; rx, ry, rw, rh = rect
    return fig.add_axes([x + rx * w, y + ry * h, rw * w, rh * h])


def clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_facecolor("none")


# ---- a: genome with sparse reads over heterozygous SNPs ----
ax = axes_in("a", (0.06, 0.08, 0.88, 0.64)); clean(ax); ax.set_xlim(0, 10); ax.set_ylim(-0.6, 7.4)
snps = np.array([0.9, 2.0, 2.6, 3.9, 4.7, 5.9, 6.6, 7.8, 8.9])
ax.plot([0.2, 9.8], [0, 0], color=INK, lw=2.2, solid_capstyle="round")
for s in snps:
    ax.plot([s, s], [-0.25, 0.25], color=INK, lw=1.2)
    ax.plot([s, s], [0.3, 7.2], color=GRID, lw=0.8, ls=(0, (2, 3)), zorder=0)
level = np.zeros(40)
for i in range(30):
    x0 = rng.uniform(0.2, 8.8); L = 0.9
    lane = next(j for j in range(1, 40) if level[j] <= x0); level[lane] = x0 + L + 0.15
    y = lane * 0.62
    if y > 7.1:
        continue
    ax.plot([x0, x0 + L], [y, y], color="#9aa5b1", lw=3.2, solid_capstyle="round", alpha=0.9)
    for s in snps[(snps > x0) & (snps < x0 + L)]:
        ax.scatter(s, y, s=26, color=H1 if rng.random() < 0.5 else H2, zorder=3, edgecolor="white", lw=0.6)
ax.text(5, -0.55, "reads over SNPs (blue / orange = the two alleles)", ha="center", va="top", fontsize=7.6, color=MUTED)

# ---- b: two haplotypes, reads pooled per haplotype per bin ----
ax = axes_in("b", (0.06, 0.30, 0.88, 0.42)); clean(ax); ax.set_xlim(0, 10); ax.set_ylim(-0.1, 2.75)
xs = np.linspace(0.6, 9.4, 12)
for y, c, lab in ((1.9, H1, "haplotype 1"), (0.6, H2, "haplotype 2")):
    ax.plot([0.3, 9.7], [y, y], color=c, lw=2.4, solid_capstyle="round")
    ax.text(9.7, y - 0.22 if y > 1 else y + 0.12, lab, fontsize=8, color=c, weight="bold", ha="right",
            va="top" if y > 1 else "bottom")
    for x in xs:
        ax.scatter(x, y, s=34, color=c, edgecolor="white", lw=0.8, zorder=3)
for x in xs:                                          # reads pooled on each haplotype (more on haplotype 1)
    for k in range(rng.integers(1, 4)):
        ax.plot([x - 0.28, x + 0.28], [1.9 + 0.17 * (k + 1)] * 2, color=H1, lw=2.4, alpha=0.55, solid_capstyle="round")
    for k in range(rng.integers(0, 3)):
        ax.plot([x - 0.28, x + 0.28], [0.6 - 0.17 * (k + 1)] * 2, color=H2, lw=2.4, alpha=0.55, solid_capstyle="round")
ax2 = axes_in("b", (0.10, 0.07, 0.80, 0.20)); ax2.set_xlim(0, 1)
ax2.barh([1], [0.62], color=H1, height=0.55); ax2.barh([1], [0.38], left=[0.62], color=H2, height=0.55)
ax2.text(0.31, 1, "h1 = 612", ha="center", va="center", color="white", fontsize=8, weight="bold")
ax2.text(0.81, 1, "h2 = 381", ha="center", va="center", color="white", fontsize=8, weight="bold")
ax2.text(0.5, 1.62, "reads per haplotype, one 500 kb bin", ha="center", fontsize=7.8, color=MUTED)
ax2.text(0.5, 0.22, "BAF = 612 / 993 = 0.62", ha="center", fontsize=8.6, color=INK, weight="bold")
clean(ax2); ax2.set_ylim(0, 2)

# ---- c: synthetic genome: depth + haplotype BAF ----
chrom_len = np.array([30, 26, 22, 20, 18, 16])          # bins per chromosome
# state per chromosome: (total, minor) at tumor fraction 0.6 -> gain, LOH (loss), CN-LOH, normal, subclonal-free
STATE = [(2, 1), (3, 1), (2, 1), (1, 0), (2, 1), (2, 0)]
NAMES = ["", "gain", "", "LOH", "", "CN-LOH"]
rho, ploidy = 0.6, 2.1
pos, logr, baf, seg = [], [], [], []
x = 0
for L, (t, m) in zip(chrom_len, STATE):
    r = np.log2((rho * t + 2 * (1 - rho)) / (rho * ploidy + 2 * (1 - rho)))
    b = (1 - rho + rho * m) / (2 * (1 - rho) + rho * t)
    px = x + np.arange(L) + 0.5
    pos.append(px); logr.append(r + rng.normal(0, 0.07, L))
    flip = rng.random(L) < 0.5
    bb = np.where(flip, b, 1 - b) + rng.normal(0, 0.03, L)
    baf.append(bb); seg.append((x, x + L, r, b)); x += L
pos, logr, baf = map(np.concatenate, (pos, logr, baf))
edges = np.cumsum(chrom_len)
axd = axes_in("c", (0.15, 0.40, 0.81, 0.27)); axb = axes_in("c", (0.15, 0.09, 0.81, 0.27))
for a_ in (axd, axb):
    a_.set_xlim(0, edges[-1]); a_.set_xticks([])
    for e in edges[:-1]:
        a_.axvline(e, color=GRID, lw=0.8)
    for s_ in ("top", "right"):
        a_.spines[s_].set_visible(False)
axd.scatter(pos, logr, s=7, color="#9aa5b1", lw=0)
axb.scatter(pos, baf, s=7, color="#9aa5b1", lw=0)
for x0, x1, r, b in seg:
    axd.plot([x0, x1], [r, r], color=BLUE, lw=2.2)
    axb.plot([x0, x1], [b, b], color=ORANGE, lw=2.0); axb.plot([x0, x1], [1 - b, 1 - b], color=ORANGE, lw=2.0)
axd.set_ylim(-0.75, 0.75); axd.set_yticks([-0.5, 0, 0.5]); axd.set_ylabel("depth\nlogR", fontsize=8.5)
axb.set_ylim(0, 1); axb.set_yticks([0, 0.5, 1]); axb.set_ylabel("haplotype\nBAF", fontsize=8.5)
for i, (nm, (x0, x1, r, b)) in enumerate(zip(NAMES, seg)):
    if nm:
        c = {"gain": ORANGE, "LOH": BLUE, "CN-LOH": TEAL}[nm]
        for a_ in (axd, axb):
            a_.axvspan(x0, x1, color=c, alpha=0.10, lw=0)
        axd.text((x0 + x1) / 2, 0.62, nm, ha="center", fontsize=8, color=c, weight="bold")
for i, (e0, e1) in enumerate(zip(np.r_[0, edges[:-1]], edges)):
    axb.text((e0 + e1) / 2, -0.12, f"chr{i + 1}", ha="center", va="top", fontsize=7.2, color=MUTED)

# ---- d: likelihood grid + allele-specific copy number ----
axg = axes_in("d", (0.18, 0.42, 0.62, 0.30))
T = np.linspace(0.05, 1, 60); P = np.linspace(1.5, 4.5, 60); TT, PP = np.meshgrid(T, P)
Z = -((TT - rho) / 0.12) ** 2 - ((PP - ploidy) / 0.35) ** 2
Z = np.maximum(Z, -((TT - 0.33) / 0.10) ** 2 - ((PP - 4.1) / 0.30) ** 2 - 1.6)   # a weaker genome-doubled mode
axg.contourf(TT, PP, Z, levels=np.linspace(-9, 0, 10), cmap="Blues", extend="min")
axg.scatter([rho], [ploidy], marker="*", s=150, color=RED, edgecolor="white", lw=0.8, zorder=3)
axg.set_xlabel("tumor fraction", fontsize=8); axg.set_ylabel("ploidy", fontsize=8)
axg.tick_params(labelsize=7); axg.set_xticks([0.2, 0.6, 1.0]); axg.set_yticks([2, 3, 4])
axg.annotate("best fit", xy=(rho, ploidy), xytext=(0.80, 3.0), fontsize=7.6, color=RED, weight="bold",
             arrowprops=dict(arrowstyle="-", color=RED, lw=0.8))
axg.annotate("genome-doubled\nalternative", xy=(0.33, 4.1), xytext=(0.62, 3.95), fontsize=7.2, color=MUTED,
             va="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
axc = axes_in("d", (0.14, 0.09, 0.82, 0.22))
axc.set_xlim(0, edges[-1]); axc.set_ylim(-0.3, 3.6); axc.set_xticks([]); axc.set_yticks([0, 1, 2, 3])
axc.tick_params(labelsize=7)
for s_ in ("top", "right"):
    axc.spines[s_].set_visible(False)
for (x0, x1, r, b), (t, m), nm in zip(seg, STATE, NAMES):
    if nm:
        axc.axvspan(x0, x1, color={"gain": ORANGE, "LOH": BLUE, "CN-LOH": TEAL}[nm], alpha=0.16, lw=0)
    axc.plot([x0, x1], [t + 0.06, t + 0.06], color=INK, lw=2.2)
    axc.plot([x0, x1], [m - 0.06, m - 0.06], color=RED, lw=2.0)
for e in edges[:-1]:
    axc.axvline(e, color=GRID, lw=0.8)
axc.set_ylabel("copies", fontsize=8)
axc.text(0, -0.55, "total", color=INK, fontsize=7.5, va="top", weight="bold")
axc.text(edges[-1] * 0.18, -0.55, "minor allele", color=RED, fontsize=7.5, va="top", weight="bold")
for nm, c, x_ in (("gain", ORANGE, 0.48), ("LOH", BLUE, 0.70), ("CN-LOH", TEAL, 0.86)):
    axc.text(edges[-1] * x_, -0.55, nm, color=c, fontsize=7.5, va="top", weight="bold")

# ---- e: results ----
ax = axes_in("e", (0.0, 0.0, 1.0, 0.70)); clean(ax); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
items = [("report.pdf", "plots + summary"), ("summary.tsv", "TF 0.60 (0.52-0.67)\nploidy 2.1, detected"),
         ("segments.tsv", "logR, BAF,\nmajor / minor CN"), ("events.tsv", "LOH, CN-LOH,\ngains, losses")]
for i, (nm, d) in enumerate(items):
    y = 0.90 - i * 0.235
    ax.add_patch(FancyBboxPatch((0.08, y - 0.17), 0.84, 0.19, boxstyle="round,pad=0.01,rounding_size=0.03",
                                fc="white", ec=GRID, lw=1.0))
    ax.add_patch(Rectangle((0.08, y - 0.17), 0.035, 0.19, color=[BLUE, PURPLE, TEAL, ORANGE][i], lw=0))
    ax.text(0.16, y - 0.015, nm, fontsize=9, weight="bold", color=INK, va="top")
    ax.text(0.16, y - 0.075, d, fontsize=7.4, color=MUTED, va="top", linespacing=1.25)

fig.savefig(OUT, dpi=170, facecolor="white")
print("->", OUT)
