"""Residual per-bin depth bias measured in this cohort's own tumor-negative samples (second pass, after run_cohort.sh).

The healthy controls behind resources/pipeline_bias_500kb.tsv are blood/plasma; CSF cfDNA keeps a small extra bias of
its own (e.g. a chr19 "gain" in tumor-free CSF). Every sample whose arm-level joint fit finds no
tumor (tumor_llr < MAX_LLR) is used as a cohort normal:
  per sample : logR_precohort from work/cnv/<s>.bins.tsv (ichorCNA correctedDepth minus the control bias), median-centred
  per bin    : cohort_bias = median over samples; spread = 1.4826 * MAD over samples
  cohort_mask: spread > SPREAD_K x the genome-wide median spread (bins too variable to correct), or |cohort_bias| > BIAS_MAX
build_bins.py subtracts cohort_bias and drops cohort_mask bins whenever the output file exists, so after reviewing the
printout, rerun bins -> seg -> joint for the cohort (refit_downstream.sh; no re-phasing).
Usage: cohort_bias.py [--summary work/joint_modes/arms/joint_summary.tsv] [--max-llr 10] [--min-n 30] [--spread-k 3]
Output: resources/cohort_bias_500kb.tsv, resources/cohort_bias_review.tsv (masked blocks), work/cohort_bias.png
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from allelic_lib import INK, INK2, SURF, style

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--summary", default=str(ROOT / "work/joint_modes/arms/joint_summary.tsv"))
ap.add_argument("--max-llr", type=float, default=10.0)
ap.add_argument("--min-n", type=int, default=30)
ap.add_argument("--spread-k", type=float, default=3.0)
ap.add_argument("--bias-max", type=float, default=0.3)
A = ap.parse_args()

sm = pd.read_csv(A.summary, sep="\t", dtype={"sample": str})      # "001" stays "001" (its bins file name)
neg = sm[(sm.tumor_llr < A.max_llr) & sm.depth_source.astype(str).str.startswith("ichorCNA")]["sample"].tolist()
cols = []
for s in neg:
    f = ROOT / f"work/cnv/{s}.bins.tsv"
    if not f.exists():
        continue
    d = pd.read_csv(f, sep="\t")
    v = d.set_index(["chrom", "start"])["logR_precohort" if "logR_precohort" in d else "logR"]
    cols.append((v - v.median()).rename(s))
print(f"{len(sm)} samples in {A.summary}; {len(cols)} tumor-negative (tumor_llr < {A.max_llr:g}) used as cohort normals")
if len(cols) < A.min_n:
    raise SystemExit(f"only {len(cols)} tumor-negative samples (< --min-n {A.min_n}): not building a cohort bias")
M = pd.concat(cols, axis=1)
n = M.notna().sum(axis=1)
M = M[n >= 0.5 * M.shape[1]]
med = M.median(axis=1)
spread = 1.4826 * M.sub(med, axis=0).abs().median(axis=1)
cut = A.spread_k * spread.median()
out = pd.DataFrame({"cohort_bias": med.round(4), "spread": spread.round(4), "n": M.notna().sum(axis=1)}).reset_index()
out["cohort_mask"] = (out.spread > cut) | (out.cohort_bias.abs() > A.bias_max)
out["_c"] = out.chrom.str[3:].astype(int); out = out.sort_values(["_c", "start"]).drop(columns="_c").reset_index(drop=True)
(ROOT / "resources").mkdir(exist_ok=True)
out.to_csv(ROOT / "resources/cohort_bias_500kb.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
chrom = out.groupby("chrom", sort=False).agg(bins=("start", "size"), median_bias=("cohort_bias", "median"),
                                             masked=("cohort_mask", "sum"))
print(f"\nper-chromosome median residual bias (logR; a 1-copy gain at TF 0.10 is ~ +0.07):")
print(chrom.round(3).T.to_string())
m = out[out.cohort_mask].copy()
m["blk"] = ((m.chrom != m.chrom.shift()) | (m.start != m.start.shift() + 500_000)).cumsum()
blocks = m.groupby("blk").agg(chrom=("chrom", "first"), start=("start", "min"), end=("start", "max"), bins=("start", "size"),
                              bias=("cohort_bias", "median"), spread=("spread", "max")).reset_index(drop=True)
blocks["end"] += 499_999
blocks.to_csv(ROOT / "resources/cohort_bias_review.tsv", sep="\t", index=False)
print(f"\ncohort_mask: {len(m)} bins ({len(m) * 0.5:.1f} Mb) in {len(blocks)} blocks "
      f"(spread > {cut:.3f} = {A.spread_k:g} x median {spread.median():.3f}, or |bias| > {A.bias_max:g}):")
print(blocks.assign(start=(blocks.start / 1e6).round(1), end=(blocks.end / 1e6).round(1)).round(3).to_string(index=False))

fig, ax = plt.subplots(2, 1, figsize=(16, 5.5), facecolor=SURF, sharex=True)
for a in ax: style(a)
off, ticks, x = 0, [], np.zeros(len(out))
for c, g in out.groupby("chrom", sort=False):
    x[g.index] = off + g.start.to_numpy(); ticks.append((off + g.start.max() / 2, c[3:])); off += g.start.max() + 5e6
col = np.where(out.cohort_mask, "#d03b3b", INK2)
ax[0].scatter(x, out.cohort_bias, s=2, c=col, lw=0); ax[0].axhline(0, color=INK, lw=0.6)
ax[0].set_ylim(-0.3, 0.3); ax[0].set_ylabel("median logR", fontsize=9, color=INK2)
ax[0].set_title(f"cohort residual bias, {len(cols)} tumor-negative samples (red = masked)", fontsize=10, color=INK, loc="left")
ax[1].scatter(x, out.spread, s=2, c=col, lw=0); ax[1].axhline(cut, color="#d03b3b", lw=0.8, ls=(0, (2, 2)))
ax[1].set_ylim(0, max(0.3, 2 * cut)); ax[1].set_ylabel("spread (1.48 MAD)", fontsize=9, color=INK2)
ax[1].set_xticks([t for t, _ in ticks]); ax[1].set_xticklabels([l for _, l in ticks], fontsize=8)
fig.tight_layout(); fig.savefig(ROOT / "work/cohort_bias.png", dpi=110, facecolor=SURF); plt.close(fig)
print(f"-> resources/cohort_bias_500kb.tsv, resources/cohort_bias_review.tsv, work/cohort_bias.png")
