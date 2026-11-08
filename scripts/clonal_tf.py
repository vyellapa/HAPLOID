"""Tumor fraction from clonal segments: the TF at which the aberrant segments fit integer copy states without subclones.

Each segment of a joint-fit segmentation (work/joint_modes/<mode>/<s>.joint_segments.tsv) carries its median logR (r_obs)
and allelic imbalance from the phased-allele HMM (theta_hat = |BAF - 0.5|). For a clonal state nA:nB at tumor fraction t:
  BAF   b = (1 - t + t nB) / (2(1 - t) + t n)              (does not depend on ploidy)
  logR  r = log2(2(1 - t) + t n) + c                          (c = normalisation, i.e. ploidy; fitted)
Grid over t and c; each segment takes its best clonal state, cost = squared z on logR and BAF, capped at CAP (a segment no
clonal state explains - subclonal, or noise - costs CAP however far off it is). Mb-weighted total + LAMBDA per copy changed
(simpler states preferred, as in PURPLE). Unlike the joint fit there are no subclonal states, so a high TF cannot be bought
by calling the losses "subclonal at 50%".
Tolerances are model-error floors, not the (tiny) statistical errors: logR 0.05, BAF 0.02 (the noise floor theta0 is added
to predicted imbalance as in joint_fit). Segments < MIN_MB, and BAF of segments without usable phasing, are left out.
Per segment the TF implied by each state from BAF alone is also written - that is the plot: clonal segments line up on
one TF.
Usage: clonal_tf.py [--mode cbs_union] [samples...]   (default: every sample in the mode's joint_summary.tsv)
Output: work/clonal_tf.tsv (per sample), work/clonal_tf_segments.tsv, work/plots/clonal_tf_<s>.png with --plot
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from allelic_lib import LINE, BW_DOTS, style_bw

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--mode", default="cbs_union")
ap.add_argument("--plot", action="store_true")
ap.add_argument("--only", choices=["logR", "BAF"], help="use one signal only (to see what each contributes)")
ap.add_argument("samples", nargs="*")
A = ap.parse_args()
J = ROOT / f"work/joint_modes/{A.mode}"
SR, ST, CAP, LAMBDA, MIN_MB = 0.05, 0.02, 9.0, 1.0, 5.0
WR, WT = (0.0, 1.0) if A.only == "BAF" else (1.0, 0.0) if A.only == "logR" else (1.0, 1.0)
STATES = [(a, b) for a in range(0, 7) for b in range(0, 4) if a >= b and 0 < a + b <= 6]
T = np.round(np.arange(0.05, 1.0001, 0.01), 2)
C = np.arange(-1.5, 0.5001, 0.01)
ABER_R, ABER_T = 0.08, 0.06          # a segment counts as aberrant beyond these (logR, theta above theta0)


def baf(t, a, b):
    return (1 - t + t * b) / (2 * (1 - t) + t * (a + b))


def fit(seg, th0):
    n = np.array([a + b for a, b in STATES]); nb = np.array([b for a, b in STATES])
    ev = np.array([abs(a - 1) + abs(b - 1) for a, b in STATES])
    w = seg.mb.to_numpy() / 10
    r = seg.r_obs.to_numpy()[:, None, None]
    th = seg.theta_hat.to_numpy()[:, None, None]; ph = seg.use_phased.to_numpy().astype(bool)
    best = (np.inf, None)
    out = np.full((len(T), len(C)), np.inf)
    for i, t in enumerate(T):
        lr = np.log2(2 * (1 - t) + t * n)                       # (state,)
        tp = np.sqrt((0.5 - (1 - t + t * nb) / (2 * (1 - t) + t * n)) ** 2 + th0 ** 2)
        zr = ((r - (lr[None, :, None] + C[None, None, :])) / SR) ** 2        # (seg, state, c)
        zt = np.where(ph[:, None, None], ((th - tp[None, :, None]) / ST) ** 2, 0.0)
        cost = np.minimum(WR * zr + WT * zt + LAMBDA * ev[None, :, None], CAP).min(axis=1)  # best state per segment
        out[i] = (w[:, None] * cost).sum(axis=0)
    i, j = np.unravel_index(np.argmin(out), out.shape)
    return out, T[i], C[j]


def assign(seg, t, c, th0):
    rows = []
    for s in seg.itertuples():
        bestc, bs = np.inf, None
        for a, b in STATES:
            n = a + b
            zr = ((s.r_obs - np.log2(2 * (1 - t) + t * n) - c) / SR) ** 2
            tp = np.hypot(0.5 - baf(t, a, b), th0)
            zt = ((s.theta_hat - tp) / ST) ** 2 if s.use_phased else 0.0
            k = WR * zr + WT * zt + LAMBDA * (abs(a - 1) + abs(b - 1))
            if k < bestc:
                bestc, bs = k, (a, b)
        rows.append((bs, bestc))
    return rows


def implied(theta, a, b):
    """TF at which state a:b gives this BAF imbalance (BAF only; None when no TF in (0, 1] does)."""
    tp = np.abs(0.5 - baf(T, a, b))
    if theta < tp.min() or theta > tp.max():
        return None
    return float(T[np.argmin(np.abs(tp - theta))])


S = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(J.glob("*.joint_summary.tsv"))]).drop_duplicates("sample", keep="last") \
    .set_index("sample")                                         # per-sample summaries (the merged file can lag behind)
samples = A.samples or S.index.tolist()
res, segrows = [], []
for s in samples:
    f = J / f"{s}.joint_segments.tsv"
    if not f.exists():
        continue
    seg = pd.read_csv(f, sep="\t")
    seg = seg[(seg.mb >= MIN_MB) & seg.r_obs.notna() & seg.theta_hat.notna()].reset_index(drop=True)
    th0 = float(S.loc[s, "theta0"]) if "theta0" in S else 0.02
    if len(seg) < 3:
        continue
    out, t, c = fit(seg, th0)
    asg = assign(seg, t, c, th0)
    n = np.array([a + b for (a, b), _ in asg])
    ploidy = float(np.average(n, weights=seg.mb))
    aber = ((seg.r_obs - np.average(seg.r_obs, weights=seg.mb)).abs() > ABER_R) | (seg.theta_hat > th0 + ABER_T)
    clonal = np.array([k < 4 for _, k in asg])
    am = seg.mb[aber].sum()
    # runner-up: best cost at a TF at least 0.15 away (another clonal reading of the same data)
    prof = out.min(axis=1)
    far = np.abs(T - t) >= 0.15
    t2 = float(T[far][np.argmin(prof[far])]) if far.any() else np.nan
    res.append(dict(sample=s, clonal_TF=t, clonal_ploidy=round(ploidy, 2),
                    ctDNA_clonal=round(t * ploidy / (t * ploidy + 2 * (1 - t)), 3),
                    aberrant_mb=round(am, 1), clonal_frac_of_aberrant=round(seg.mb[aber & clonal].sum() / am, 2) if am > 0 else np.nan,
                    runner_up_TF=t2, d_cost_runner_up=round(float(prof[far].min() - prof.min()), 1) if far.any() else np.nan,
                    joint_TF=S.loc[s, "joint_TF"], joint_ploidy=S.loc[s, "joint_ploidy"], ichor_TF=S.loc[s, "ichor_TF"],
                    tumor_detected=S.loc[s, "tumor_detected"]))
    for r, ((a, b), k), ab in zip(seg.itertuples(), asg, aber):
        imp = {f"tf_{x}{y}": implied(r.theta_hat - th0 if r.theta_hat > th0 else 0, x, y)
               for x, y in [(1, 0), (2, 0), (2, 1), (3, 1), (3, 0), (4, 1)]} if r.use_phased else {}
        segrows.append(dict(sample=s, chrom=r.chrom, start=r.start, end=r.end, mb=r.mb, r_obs=round(r.r_obs, 3),
                            theta_hat=r.theta_hat, aberrant=bool(ab), state=f"{a}:{b}", cost=round(k, 2),
                            clonal=k < 4, **imp))
    if A.plot:
        g = pd.DataFrame([x for x in segrows if x["sample"] == s and x["aberrant"]])
        fig, ax = plt.subplots(figsize=(8, 0.28 * max(len(g), 6) + 1.2), facecolor="white"); style_bw(ax)
        col = {"tf_10": LINE["blue"], "tf_20": LINE["teal"], "tf_21": LINE["red"], "tf_31": LINE["orange"],
               "tf_30": "#555555", "tf_41": "#999999"}
        lab = {"tf_10": "1:0 loss", "tf_20": "2:0 CN-LOH", "tf_21": "2:1 gain", "tf_31": "3:1 gain",
               "tf_30": "3:0", "tf_41": "4:1"}
        for y, row in enumerate(g.itertuples()):
            for k_, cl in col.items():
                v = getattr(row, k_, None)
                if v is not None and np.isfinite(v):
                    ax.scatter(v, y, s=36, color=cl, zorder=3, edgecolor="white", lw=0.5)
        ax.axvline(t, color="black", lw=1.2); ax.axvline(S.loc[s, "joint_TF"], color=LINE["blue"], lw=1, ls=(0, (3, 2)))
        if np.isfinite(S.loc[s, "ichor_TF"]):
            ax.axvline(S.loc[s, "ichor_TF"], color=LINE["red"], lw=1, ls=(0, (1, 2)))
        ax.set_yticks(range(len(g))); ax.set_yticklabels([f"{r.chrom[3:]}:{r.start/1e6:.0f}-{r.end/1e6:.0f}  {r.state}"
                                                           for r in g.itertuples()], fontsize=7)
        ax.set_xlim(0, 1.02); ax.set_xlabel("tumor fraction implied by the segment's BAF", fontsize=9)
        ax.set_title(f"{s}: clonal TF {t:.2f} (black) · joint fit {S.loc[s, 'joint_TF']:.2f} (blue dashed) · "
                     f"ichorCNA {S.loc[s, 'ichor_TF']:.2f} (red dotted)", fontsize=9, loc="left")
        from matplotlib.lines import Line2D
        ax.legend(handles=[Line2D([], [], marker="o", ls="", color=v, label=lab[k]) for k, v in col.items()],
                  fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        o = ROOT / f"work/plots/clonal_tf_{s}.png"; o.parent.mkdir(exist_ok=True, parents=True)
        fig.savefig(o, dpi=130, facecolor="white", bbox_inches="tight"); plt.close(fig)

R = pd.DataFrame(res); G = pd.DataFrame(segrows)
fo, fs = ROOT / f"work/clonal_tf{'_' + A.only if A.only else ''}.tsv", ROOT / f"work/clonal_tf{'_' + A.only if A.only else ''}_segments.tsv"
if A.samples:                               # a subset updates its own rows and keeps the rest of the table
    if fo.exists():
        R = pd.concat([pd.read_csv(fo, sep="\t").query("sample not in @samples"), R])
    if fs.exists():
        G = pd.concat([pd.read_csv(fs, sep="\t").query("sample not in @samples"), G])
R.to_csv(fo, sep="\t", index=False); G.to_csv(fs, sep="\t", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 300)
print(R.round(2).to_string(index=False))
print(f"-> {ROOT / 'work/clonal_tf.tsv'}, {ROOT / 'work/clonal_tf_segments.tsv'}")
