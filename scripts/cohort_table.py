"""One row per sample from the joint-fit modes, following the dilution results:
  detection, TF, ploidy        : arms mode (cleanest detection, TF within ~0.01 down to TF 0.10)
  events (LOH / CN-LOH / gains): cbs_union when detected with TF >= CBS_MIN_TF (most sensitive, few false calls there),
                                 otherwise arms (arm-level only; very low TF gives too little signal for breakpoints)
  ichorCNA TF / ploidy         : alongside, for comparison
ctDNA_fraction is the tumor share of cfDNA; unlike the cell fraction (TF) it does not depend on the ploidy choice, so use it
when `ambiguous` is True (a solution at a different ploidy within 10 log-likelihood units).
Allele-data QC: dup_rate / unique_reads_M from the markdup stats (work/bam/<s>.auto.markdup_stats.txt), baf_bins = bins with
>= 100 haplotype reads. Very high-duplicate (low-input) libraries keep ichorCNA depth (it counts duplicates) but have almost no
allele reads: allele_qc = "low" (baf_bins < --min-baf-bins) means the fit is close to depth-only, so "no LOH" there is not
a negative result.
fit_flag = "uniform_imbalance" when > 90% of the genome sits in one state other than 1:1 / 2:2 (e.g. all 2:1 at ploidy 3)
with flat depth: no tumor makes that, but a second person's DNA (excess het sites, smeared BAF) or noisy allele counts do.
Treat those detections as QC failures until checked.
clonal_TF / ctDNA_clonal: from clonal_tf.py (the TF at which the aberrant segments fit integer states with no subclones),
merged in when work/clonal_tf.tsv exists - run clonal_tf.py first. Where it and joint_TF differ, the joint fit has usually
explained a set of clonal events as "subclonal" at a higher TF (see work/plots/clonal_tf_<s>.png with --plot).
Refitted samples: TF, ploidy, ctDNA and events come from the refit, tf_basis says which and what it cost (dLL_constraint):
  work/joint_modes/cbs_union_user/          values you gave (refit_with.sh) - first priority
  work/joint_modes/cbs_union_ichorploidy/   ploidy held at ichorCNA's (resolve_discordant.sh)
TF_lo / TF_hi, ctDNA_lo / ctDNA_hi: segment-bootstrap 95% interval from the same fit as the reported estimate: arms,
or the constrained refit when active, detected samples only. With few informative arms, intervals can be wide because
a resample may miss the events that anchor TF. Cached refits are used only for samples in their current override table.
TF_uncertain = interval wider than --max-tf-width (default 0.2): the estimate depends on which events drive it.
Usage: cohort_table.py [--min-tf 0.10] [--min-baf-bins 2000] [--max-tf-width 0.2]
Output: work/cohort_summary.tsv, work/cohort_events.tsv (one row per called event, with the mode it came from)
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from allelic_lib import call_events
from summary_tables import read_summaries

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument("--min-tf", type=float, default=0.10)
ap.add_argument("--min-baf-bins", type=int, default=2000)
ap.add_argument("--max-tf-width", type=float, default=0.2); A = ap.parse_args()
J = ROOT / "work/joint_modes"
arms = read_summaries(J / "arms").set_index("sample")
cbs = read_summaries(J / "cbs_union").set_index("sample")


def qc(s):
    q = dict(dup_rate=None, unique_reads_M=None, baf_bins=None, hap_reads_per_bin=None)
    f = ROOT / f"work/bam/{s}.auto.markdup_stats.txt"
    if f.exists():
        kv = dict(l.split(": ", 1) for l in f.read_text().splitlines() if ": " in l)
        r, d = int(kv["READ"]), int(kv["DUPLICATE TOTAL"])
        q.update(dup_rate=round(d / r, 3), unique_reads_M=round((r - d) / 1e6, 1))
    f = ROOT / f"work/cnv/{s}.bins.tsv"
    if f.exists():
        b = pd.read_csv(f, sep="\t", usecols=["baf", "n"])
        q.update(baf_bins=int(b.baf.notna().sum()), hap_reads_per_bin=int(b.n.median()))
    q["allele_qc"] = "" if q["baf_bins"] is None else "ok" if q["baf_bins"] >= A.min_baf_bins else "low"
    return q


def fit_flag(seg, det):
    if not det:
        return ""
    st = seg.assign(k=seg.nA.astype(str) + ":" + seg.nB.astype(str) + seg.sub_to.fillna("").astype(str)).groupby("k").mb.sum()
    return "uniform_imbalance" if st.max() > 0.9 * st.sum() and st.idxmax() not in ("1:1", "2:2") else ""


REFITS = [("cbs_union_user", "your values"), ("cbs_union_ichorploidy", "ploidy held at ichorCNA's")]
active_refits = {}
for folder, f in [("cbs_union_user", ROOT / "work/user_constraints.tsv"),
                  ("cbs_union_ichorploidy", ROOT / "work/fit_overrides_ichor_ploidy.tsv")]:
    active_refits[folder] = set(pd.read_csv(f, sep="\t", dtype=str)["sample"]) if f.exists() else set()
rows, events = [], []
for s, a in arms.iterrows():
    det = str(a.tumor_detected) == "True"
    use = "cbs_union" if det and a.joint_TF >= A.min_tf and s in cbs.index else "arms"
    held = next(((d, why) for d, why in REFITS
                 if det and s in active_refits[d] and (J / d / f"{s}.joint_summary.tsv").exists()), None)
    if held:                                              # refitted: events from that fit too
        use = "cbs_union"
    folder = held[0] if held else use
    seg = pd.read_csv(J / f"{folder}/{s}.joint_segments.tsv", sep="\t")
    ev = call_events(seg) if det else pd.DataFrame(columns=["chrom", "start", "end", "type", "mb", "focal"])
    src = pd.read_csv(J / f"{folder}/{s}.joint_summary.tsv", sep="\t").iloc[0] if folder != use \
        else cbs.loc[s] if use == "cbs_union" else a
    ok = seg[seg.call_pass.astype(str) == "True"] if det else seg.iloc[0:0]
    sub = ok.sub_to.fillna("").astype(str)                # LOH also when both subclonal states lack the B allele
    loh = ok[(ok.nB == 0) & sub.map(lambda u: u == "" or u.split("@")[0].split(":")[1] == "0").astype(bool)]
    flag = fit_flag(pd.read_csv(J / f"arms/{s}.joint_segments.tsv", sep="\t"), det)
    rows.append(dict(sample=s, tumor_detected=det, fit_flag=flag, tumor_llr=a.tumor_llr,
                     TF=a.joint_TF if det else "not_detected", ctDNA_fraction=a.get("ctDNA_fraction", "") if det else "not_detected",
                     ploidy=a.joint_ploidy if det else "", ambiguous=a.ambiguous if det else "",
                     dLL_other_ploidy=a.dLL_other_ploidy if det else "",
                     ichor_TF=a.ichor_TF, ichor_ploidy=a.ichor_ploidy, events_from=use if det else "",
                     loh_mb=round(loh.mb.sum(), 1), cnloh_mb=src.cnloh_mb if det else 0.0,
                     n_loh=int((ev.type == "LOH").sum()), n_gain=int((ev.type == "gain").sum()),
                     n_loss=int((ev.type == "loss").sum()), n_focal=int(ev.focal.astype(bool).sum()) if len(ev) else 0,
                     constraints=a.get("constraints", ""), depth_source=a.depth_source, **qc(s)))
    ci = src if held else a
    lo, hi = ci.get("TF_lo", np.nan), ci.get("TF_hi", np.nan)
    rows[-1].update(TF_lo=lo if det else "", TF_hi=hi if det else "",
                    ctDNA_lo=ci.get("ctDNA_lo", np.nan) if det else "", ctDNA_hi=ci.get("ctDNA_hi", np.nan) if det else "",
                    TF_uncertain=bool(hi - lo > A.max_tf_width) if det and pd.notna(lo) and pd.notna(hi) else "",
                    tf_basis="arms" if det else "")
    if held:
        rows[-1].update(TF=src.joint_TF, ploidy=src.joint_ploidy, ctDNA_fraction=src.ctDNA_fraction,
                        tf_basis=f"cbs_union, {held[1]}: {src.constraints} (dLL_constraint {src.dLL_constraint})")
    for e in ev.itertuples():
        events.append(dict(sample=s, mode=folder, chrom=e.chrom, start=int(e.start), end=int(e.end), type=e.type,
                           mb=round(e.mb, 1), focal=e.focal))
R = pd.DataFrame(rows)
if (ROOT / "work/clonal_tf.tsv").exists():
    C = pd.read_csv(ROOT / "work/clonal_tf.tsv", sep="\t", dtype={"sample": str})[["sample", "clonal_TF", "clonal_ploidy", "ctDNA_clonal",
                                                           "clonal_frac_of_aberrant", "runner_up_TF", "d_cost_runner_up"]]
    R = R.merge(C, on="sample", how="left")
    for c in C.columns[1:]:                               # only meaningful when tumor is detected
        R[c] = R[c].where(R.tumor_detected)
R.to_csv(ROOT / "work/cohort_summary.tsv", sep="\t", index=False)
E = pd.DataFrame(events, columns=["sample", "mode", "chrom", "start", "end", "type", "mb", "focal"])
E.to_csv(ROOT / "work/cohort_events.tsv", sep="\t", index=False)
d = R[R.tumor_detected]
print(f"{len(R)} samples: {len(d)} tumor detected ({(d.ambiguous.astype(str) == 'True').sum()} ploidy-ambiguous: "
      f"use ctDNA_fraction), events from cbs_union for {(d.events_from == 'cbs_union').sum()}")
u = R[R.fit_flag == "uniform_imbalance"]
if len(u):
    print(f"fit_flag uniform_imbalance (whole genome in one non-diploid state, flat depth; check for contamination): "
          f"{', '.join(u['sample'])}")
low = R[R.allele_qc == "low"]
if len(low):
    print(f"allele_qc low (< {A.min_baf_bins} BAF bins; near depth-only, LOH not assessable): {len(low)} samples")
    print(low[["sample", "dup_rate", "unique_reads_M", "baf_bins", "tumor_detected", "TF", "ichor_TF"]].to_string(index=False))
if len(E):
    rec = E.groupby(["chrom", "type"])["sample"].nunique().sort_values(ascending=False).head(15)
    print("most frequent called events (samples per chromosome and type; a recurrent artifact shows up here too):")
    print(rec.to_string())
print(f"-> {ROOT / 'work/cohort_summary.tsv'}, {ROOT / 'work/cohort_events.tsv'}")
