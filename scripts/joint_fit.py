"""Joint purity / ploidy / allele-specific copy-number fit from three evidence terms (the "all-sites" design).

For every (rho, psi) on an ASCAT-style grid and every segment, the best integer state (nA >= nB >= 0, 1 <= nA+nB <= 8)
is chosen by log-likelihood; the grid point with the highest total wins.
  depth   : segment median logR ~ N(log2((2(1-rho) + rho n) / (2(1-rho) + rho psi)), sd_r)       [ASCAT, gamma = 1]
  phased  : allele-HMM log-likelihood curve over theta (phase switches marginalised; numbat), read at
            theta = 0.5 - b, b = (1-rho + rho nB) / (2(1-rho) + rho n); flattened so its effective SE >= SE_THETA_FLOOR.
            DROPPED for segments with het density D < D_DROPOUT: imputation removed LOH hets there and the survivors
            are a biased (more balanced) subset.
  all-sites both-allele rate : every panel site with >= 2 reads, imputation not used; observed rate (normalised
            to balanced segments) ~ N(R(b), se) with R(b) = sum 2pq(1 - b^n - (1-b)^n) / sum 2pq(1 - 2 * 0.5^n),
            se = max(phi * binomial SE, empirical window SD)  [FACETS-CH-style null]. Carries LOH where hets drop out.
            Noise floor: balanced segments measure theta ~ 0.02 rather than 0 (residual overdispersion / allelic bias), so
            the curve is read at sqrt(theta_pred^2 + theta0^2), theta0 = median theta_hat over segments balanced on every
            signal; otherwise the fit "explains" that noise with spurious low-fraction subclonal states.
  prior   : - LAMBDA * PURPLE event distance per segment (min(single-event, WGD) distance; purple README 462-488).
Segments (env HAPLOID_SEG_MODE): "ascat" (default) = ASCAT segments, ASPCF breakpoints when the ASCAT fit failed;
"arms" = chromosome arms (centromere +/- 2 Mb removed; most hets per event, for low TF); "ascat_arms" = ASCAT segments
split at centromeres; "cbs_union" = CNACS-style union of breakpoints from separate CBS runs on logR and on mirrored
phased BAF (work/cnv/<s>.cbs_segments.tsv, run_ascat.R), within arms; "cbs_ichor" = cbs_union plus the ichorCNA segment
breakpoints (ICHOR_SEG; its HMM finds weaker depth steps), an ichorCNA breakpoint dropped within 1 Mb of an existing one.
seg_source in the summary. Output folder: env HAPLOID_JOINT_OUT (default work/joint).
Switch scale nu: grid-MLE (override per sample with env HAPLOID_NU="S1=4,S2=2") on strongly imbalanced segments (phased_lib.fit_nu).
Summary: TF_report = rho (tumor CELL fraction); ctDNA_fraction = rho*ploidy / (rho*ploidy + 2(1-rho)) = tumor share of
cfDNA, which genome-doubled and diploid explanations share - report it when `ambiguous` (other ploidy within 10 LL).
Outputs: work/joint/<sample>.joint_segments.tsv, <sample>.joint_summary.tsv (merged into joint_summary.tsv),
         <sample>.joint_sunrise.png
User constraints (refit with what you know, e.g. from tissue or a clearly balanced chromosome): one row per sample in
resources/fit_overrides.tsv (env HAPLOID_OVERRIDES), columns sample, purity, ploidy, diploid (blank = free):
  purity  "0.35" (= 0.30-0.40) or "0.3-0.5";  ploidy  "2" (= 1.75-2.25) or "3.5-4.5"  -> fit restricted to that window
  diploid "chr2,chr8p,chr11q"  -> segments centred there are held at 1:1, which anchors absolute copy number
Detection (tumor_llr) still uses the free fit; dLL_constraint = constrained - free best log-likelihood (below about -10
the data disagree with the constraint).
Usage: joint_fit.py <sample> [<sample> ...]
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from allelic_lib import CEN_PAD, INK, INK2, RES, SURF, ichor_ploidy, in_mask, params, segments, site_table, style
from phased_lib import THETA, fit_nu, forward_ll, load_hapcounts

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("HAPLOID_JOINT_OUT", ROOT / "work/joint")); OUT.mkdir(parents=True, exist_ok=True)
SEG_MODE = os.environ.get("HAPLOID_SEG_MODE", "ascat")
RHO = np.round(np.arange(0.05, 1.0001, 0.01), 2)
PSI = np.round(np.arange(1.5, 6.5001, 0.05), 2)
CLONAL = [(a, b) for a in range(0, 9) for b in range(0, a + 1) if 1 <= a + b <= 8]
SUB_F = (0.3, 0.5, 0.7)
# state = (nA, nB, nA2, nB2, f): a clonal state, or a mixture with a one-copy neighbour in a fraction f of tumor cells
STATES = [(a, b, a, b, 0.0) for a, b in CLONAL]
for a, b in CLONAL:
    for a2, b2 in ((a - 1, b), (a, b - 1), (a + 1, b), (a, b + 1)):
        if a2 >= 0 and b2 >= 0 and 1 <= a2 + b2 <= 8:
            STATES += [(a, b, a2, b2, f) for f in SUB_F]
SD_R_FLOOR, SE_THETA_FLOOR, D_DROPOUT, LAMBDA, MU_SUB, MIN_SEG = 0.02, 0.01, 0.5, 2.0, 4.0, 5e6
VETO = 2.0          # a track contradicts a call when it favours the background state by more than this (log-lik)
SUB_MIN, CNLOH_MIN = 10.0, 10.0   # total evidence needed for subclonal calls / allele evidence for copy-neutral LOH
EVENT_MIN = float(os.environ.get("HAPLOID_EVENT_MIN", 10))   # summed evidence of a merged event (adjacent passing
                  # segments with the same call type); weaker events are single-segment depth wiggles at low TF
# focal events (FOCAL_MIN <= size < MIN_SEG): never used for the purity/ploidy fit; tested afterwards at the fitted
# solution and kept only with strong evidence in BOTH depth and alleles (homozygous deletion: hets must vanish)
FOCAL_MIN = 1e6
FOCAL_DEPTH = float(os.environ.get("HAPLOID_FOCAL_DEPTH", 15))
FOCAL_ALLELE = float(os.environ.get("HAPLOID_FOCAL_ALLELE", 10))
FOCAL_HOMDEL_D = 0.3
FOCAL_AVAIL = 0.5   # unmasked panel-site density (vs genome median) needed to judge alleles in a focal segment
TUMOR_LLR_MIN = float(os.environ.get("HAPLOID_TUMOR_LLR_MIN", 20))   # best fit must beat the all-1:1 null by this much
                         # (provisional; dilution nulls reached 11.6, so ~15 may do - calibrate on the cohort)
NU = {k: int(v) for k, v in (x.split("=") for x in os.environ.get("HAPLOID_NU", "").split(",") if x)}  # e.g. "sample1=4"
OVERRIDES = Path(os.environ.get("HAPLOID_OVERRIDES", ROOT / "resources/fit_overrides.tsv"))
B_OF_K = 0.5 - THETA                            # minor-allele fraction for each theta-grid index


def event_dist(a, b):
    return min(abs(a - 1) + abs(b - 1), 1 + abs(a - 2) + abs(b - 2))


def theta_se(ll, k):
    """Half-width at a 0.5 log-likelihood drop around the maximum."""
    drop = ll[k] - ll
    right = np.where(drop[k:] >= 0.5)[0]; left = np.where(drop[:k + 1][::-1] >= 0.5)[0]
    w = [THETA[k + right[0]] - THETA[k]] if len(right) else []
    w += [THETA[k] - THETA[k - left[0]]] if len(left) else []
    return min(w) if w else 0.5


def arms(bins):
    """Chromosome arms spanned by the depth bins, centromere +/- CEN_PAD removed (acrocentric p arms drop out)."""
    cb = pd.read_csv(RES / "cytoband_b38.bed", sep="\t", header=None)
    cen = cb[cb[4] == "acen"].groupby(0).agg({1: "min", 2: "max"})
    rows = []
    for c, g in bins.groupby("chrom"):
        lo, hi = int(g.start.min()), int(g.start.max() + 500_000)
        rows += [(c[3:], lo, min(hi, cen.loc[c, 1] - CEN_PAD)), (c[3:], max(lo, cen.loc[c, 2] + CEN_PAD), hi)]
    a = pd.DataFrame(rows, columns=["chr", "startpos", "endpos"])
    return a[a.endpos > a.startpos]


def build_segments(s):
    d = ROOT / f"work/cnv/{s}_ascat"
    src = "ascat" if (d / f"{s}.ascat_segments.tsv").exists() else "aspcf"
    bins = pd.read_csv(ROOT / f"work/cnv/{s}.bins.tsv", sep="\t")
    if SEG_MODE == "arms":
        seg, src = arms(bins), "arms"
    elif SEG_MODE in ("cbs_union", "cbs_ichor"):
        cbs = pd.read_csv(ROOT / f"work/cnv/{s}.cbs_segments.tsv", sep="\t")
        cbs["chrom"] = cbs.chrom.astype(str).str.replace("chr", "")
        if SEG_MODE == "cbs_ichor":
            ich = segments(s)
            ich = ich.assign(chrom=ich.chrom.astype(str).str.replace("chr", ""), start=ich.start.astype(int))
        rows = []
        for a in arms(bins).itertuples():
            x = cbs[(cbs.chrom == a.chr) & (cbs["loc.start"] > a.startpos) & (cbs["loc.start"] < a.endpos)]
            edges = [a.startpos] + sorted(set(x["loc.start"].astype(int))) + [a.endpos]
            if SEG_MODE == "cbs_ichor":
                for e in ich[(ich.chrom == a.chr) & (ich.start > a.startpos) & (ich.start < a.endpos)].start:
                    if min(abs(e - np.array(edges))) >= FOCAL_MIN:
                        edges = sorted(edges + [int(e)])
            rows += [(a.chr, e0, e1) for e0, e1 in zip(edges[:-1], edges[1:]) if e1 > e0]
        seg, src = pd.DataFrame(rows, columns=["chr", "startpos", "endpos"]), SEG_MODE
    else:
        seg = pd.read_csv(d / f"{s}.{src}_segments.tsv", sep="\t")
        seg["chr"] = seg.chr.astype(str).str.replace("chr", "")
        if SEG_MODE == "ascat_arms":            # split ASCAT segments at centromeres
            A = arms(bins)
            seg = pd.DataFrame([(a.chr, max(r.startpos, a.startpos), min(r.endpos, a.endpos))
                                for r in seg.itertuples() for a in A[A.chr == r.chr].itertuples()
                                if min(r.endpos, a.endpos) > max(r.startpos, a.startpos)],
                               columns=["chr", "startpos", "endpos"])
            src += "+arms"
    seg = seg[(seg.endpos - seg.startpos) >= FOCAL_MIN].copy()
    seg["chrom"] = "chr" + seg.chr.astype(str)
    hc = load_hapcounts(s)
    big = seg[(seg.endpos - seg.startpos) >= MIN_SEG]
    nu = NU[s] if s in NU else fit_nu(hc, list(zip(big.chrom, big.startpos, big.endpos)))[0]
    sites = {c: site_table(s, int(c[3:])) for c in sorted(set(seg.chrom), key=lambda x: int(x[3:]))}
    panel = {}                                   # unmasked 1kGP panel positions: how many sites COULD carry alleles
    for c in sites:
        p_ = pd.read_csv(ROOT / f"resources/sites_hg38/{c}.sites.tsv.gz", sep="\t", header=None, usecols=[1]).iloc[:, 0].to_numpy()
        panel[c] = np.sort(p_[~in_mask(c, p_)])
    sd_bin = 1.4826 * np.nanmedian(np.abs(np.diff(bins.logR))) / np.sqrt(2)
    rows = []
    for r in seg.itertuples():
        b = bins[(bins.chrom == r.chrom) & bins.pos.between(r.startpos, r.endpos)]
        b = b[np.isfinite(b.logR)]
        x = sites[r.chrom]; x = x[x.pos.between(r.startpos, r.endpos)]
        h = hc[(hc.chrom == r.chrom) & hc.pos.between(r.startpos, r.endpos)].sort_values("pos")
        focal = (r.endpos - r.startpos) < MIN_SEG
        if len(b) < (2 if focal else 4) or len(x) < (100 if focal else 500):
            continue
        ll = None
        if len(h) >= (50 if focal else 200):
            dcm = np.diff(h.cm.to_numpy(), prepend=h.cm.iat[0])
            ll = forward_ll(h.h1.to_numpy(), h.h2.to_numpy(), dcm, nu)
            ll = ll - ll.max()
        p0 = x.exp.to_numpy()
        pred = np.array([(2 * x.af * (1 - x.af) * (1 - bb ** x.n - (1 - bb) ** x.n)).sum() for bb in B_OF_K]) / p0.sum()
        L = r.endpos - r.startpos
        rows.append(dict(chrom=r.chrom, start=int(r.startpos), end=int(r.endpos), mb=L / 1e6, bins=len(b),
                         r_obs=float(b.logR.median()), sd_r=max(sd_bin / np.sqrt(len(b)), SD_R_FLOOR),
                         obs=float(x.obs.sum()), exp=float(p0.sum()), pvar=float(np.sum(p0 * (1 - p0))),
                         R_pred=pred, ll=ll, hets_mb=len(h) / (L / 1e6), sites=len(x), focal=focal,
                         avail_mb=(np.searchsorted(panel[r.chrom], r.endpos) - np.searchsorted(panel[r.chrom], r.startpos)) / (L / 1e6)))
    S = pd.DataFrame(rows)
    # het-density reference: median over (non-focal) segments; normalisation set = balanced, non-depleted segments
    S["D"] = S.hets_mb / S[~S.focal].hets_mb.median()
    S["A"] = S.avail_mb / S[~S.focal].avail_mb.median()     # unmasked site availability vs genome
    S["theta_hat"] = [THETA[int(np.argmax(l))] if l is not None else np.nan for l in S.ll]
    S["theta_se"] = [theta_se(l, int(np.argmax(l))) if l is not None else np.nan for l in S.ll]
    bal = S[(S.theta_hat <= 0.02) & (S.D > 0.7) & ~S.focal]
    if len(bal) < 3:                      # almost nothing balanced: normalise on the least imbalanced half
        bal = S[(S.theta_hat <= S[~S.focal].theta_hat.median()) & ~S.focal]
    norm = bal.obs.sum() / bal.exp.sum()
    rate = (bal.obs / bal.exp) / norm
    se_bin = np.sqrt(bal.pvar) / bal.exp / norm
    phi = float(np.median(np.abs(rate - 1) / se_bin) * 1.4826) if len(bal) >= 3 else 3.0
    S["R_obs"] = S.obs / S.exp / norm
    S["se_R"] = np.maximum(np.maximum(phi, 1) * np.sqrt(S.pvar) / S.exp / norm, 0.03 / np.sqrt(np.maximum(S.mb / 10, 1)))
    S["use_phased"] = [(l is not None) and (d >= D_DROPOUT) for l, d in zip(S.ll, S.D)]
    S.attrs["depth_source"] = str(bins.logR_source.iat[0]) if "logR_source" in bins else "unknown"
    calm = S[S.use_phased & S.R_obs.between(0.95, 1.05) & S.D.between(0.8, 1.25) & (S.theta_hat <= 0.05) & ~S.focal]
    S.attrs["theta0"] = float(calm.theta_hat.median()) if len(calm) >= 3 else 0.02
    return S, norm, phi, nu, src


def score(S):
    """Total log-likelihood over the (rho, psi) grid plus the best state per segment."""
    R_, P_ = np.meshgrid(RHO, PSI, indexing="ij")
    total = np.zeros_like(R_); best = []
    th0 = S.attrs.get("theta0", 0.0)
    for s in S.itertuples():
        best_ll = np.full(R_.shape, -np.inf); best_st = np.zeros(R_.shape, int)
        for i, (a0, b0, a2, b2, f) in enumerate(STATES):
            if getattr(s, "force11", False) and (a0, b0, f) != (1, 1, 0.0):
                continue
            a, bmin = a0 + f * (a2 - a0), b0 + f * (b2 - b0)     # effective (possibly fractional) copies
            if bmin > a:
                a, bmin = bmin, a
            n = a + bmin
            r_pred = np.log2((2 * (1 - R_) + R_ * n) / (2 * (1 - R_) + R_ * P_))
            ll = -0.5 * ((s.r_obs - r_pred) / s.sd_r) ** 2
            b = (1 - R_ + R_ * bmin) / (2 * (1 - R_) + R_ * n)
            k = np.clip(np.rint((0.5 - b) / 0.005).astype(int), 0, len(THETA) - 1)
            ll += -0.5 * ((s.R_obs - s.R_pred[k]) / s.se_R) ** 2
            if s.use_phased:
                tau = max(1.0, (SE_THETA_FLOOR / max(s.theta_se, 1e-4)) ** 2)
                kp = np.clip(np.rint(np.sqrt((0.5 - b) ** 2 + th0 ** 2) / 0.005).astype(int), 0, len(THETA) - 1)
                ll += s.ll[kp] / tau
            ll -= LAMBDA * event_dist(max(a0, b0), min(a0, b0)) + (MU_SUB if f > 0 else 0.0)
            upd = ll > best_ll
            best_ll[upd] = ll[upd]; best_st[upd] = i
        total += best_ll; best.append(best_st)
    return total, best


def null_ll(S):
    """Log-likelihood of 'no tumor': every segment 1:1 (theta 0, both-allele rate 1), best common logR offset.
    Purity is unidentifiable then, so a sample whose best fit does not beat this clearly gets TF = not detected."""
    th0 = S.attrs.get("theta0", 0.0)
    kp = int(np.rint(th0 / 0.005))
    base = sum(-0.5 * ((s.R_obs - s.R_pred[0]) / s.se_R) ** 2 + (s.ll[kp] / max(1.0, (SE_THETA_FLOOR / max(s.theta_se, 1e-4)) ** 2)
               if s.use_phased else 0.0) for s in S.itertuples())
    return base + max(sum(-0.5 * ((s.r_obs - c) / s.sd_r) ** 2 for s in S.itertuples()) for c in np.arange(-0.3, 0.301, 0.005))


def call_evidence(S, st, rho, psi, base):
    """Per segment, log-likelihood of the called state minus the background (most common clonal) state, per track:
    depth, phased allele HMM, all-sites both-allele rate. Then the call filter:
      a track vetoes a call it contradicts (llr < -VETO): only possible when the call predicts a signal on that
      track that it should have seen, so copy-neutral LOH is never vetoed by flat depth;
      subclonal calls need both tracks >= 0 and total >= SUB_MIN;
      copy-neutral LOH (nB = 0, same total as background) needs allele evidence >= CNLOH_MIN, neither allele signal
      against it, and more hets than a germline run of homozygosity would leave (ROH flag)."""
    th0 = S.attrs.get("theta0", 0.0)

    def pred(n, bm):
        r = np.log2((2 * (1 - rho) + rho * n) / (2 * (1 - rho) + rho * psi))
        b = (1 - rho + rho * bm) / (2 * (1 - rho) + rho * n)
        k = int(np.clip(np.rint((0.5 - b) / 0.005), 0, len(THETA) - 1))
        kp = int(np.clip(np.rint(np.sqrt((0.5 - b) ** 2 + th0 ** 2) / 0.005), 0, len(THETA) - 1))
        return r, k, kp
    rb, kb, kpb = pred(base[0] + base[1], base[1])
    rows = []
    for s, x in zip(S.itertuples(), st):
        a0, b0, a2, b2, f = x
        a, bm = a0 + f * (a2 - a0), b0 + f * (b2 - b0)
        if bm > a:
            a, bm = bm, a
        r, k, kp = pred(a + bm, bm)
        d = -0.5 * ((s.r_obs - r) ** 2 - (s.r_obs - rb) ** 2) / s.sd_r ** 2
        al = -0.5 * ((s.R_obs - s.R_pred[k]) ** 2 - (s.R_obs - s.R_pred[kb]) ** 2) / s.se_R ** 2
        ph = ((s.ll[kp] - s.ll[kpb]) / max(1.0, (SE_THETA_FLOOR / max(s.theta_se, 1e-4)) ** 2)) if s.use_phased else 0.0
        is_call = (max(a0, b0), min(a0, b0)) != base or f > 0
        reason = []
        if is_call and s.focal:
            if s.A < FOCAL_AVAIL:
                reason.append("focal_masked_region")      # too few unmasked sites to judge alleles
            if f > 0:
                reason.append("focal_subclonal")
            if d < FOCAL_DEPTH:
                reason.append("focal_depth_weak")
            if a0 + b0 == 0:
                if not s.D < FOCAL_HOMDEL_D:
                    reason.append("focal_homdel_hets_present")
            elif ph + al < FOCAL_ALLELE:
                reason.append("focal_allele_weak")
        if is_call:
            if d < -VETO:
                reason.append("depth_contradicts")
            if ph + al < -VETO:
                reason.append("allele_contradicts")
            if f > 0 and (d < 0 or ph + al < 0 or d + ph + al < SUB_MIN):
                reason.append("subclonal_weak")
            if min(a0, b0) == 0 and f == 0 and a0 + b0 == sum(base):
                if ph + al < CNLOH_MIN or ph < -1 or al < -1:
                    reason.append("cnloh_allele_weak")
                if s.R_obs < s.R_pred[k] - 3 * s.se_R:
                    reason.append("roh_like")
        rows.append(dict(llr_depth=round(d, 2), llr_phased=round(ph, 2), llr_allsites=round(al, 2), is_call=is_call,
                         call_pass=is_call and not reason, filter_reason=",".join(reason)))
    return pd.DataFrame(rows, index=S.index)


def window(v, half):
    """ "0.35" -> (0.35 - half, 0.35 + half); "0.3-0.5" -> (0.3, 0.5)."""
    v = str(v).strip()
    if "-" in v:
        a, b = v.split("-", 1)
        return float(a), float(b)
    return float(v) - half, float(v) + half


def constraints(s):
    """(allowed (rho, psi) grid mask, [(chrom, start, end)] held at 1:1, description) from OVERRIDES for sample s."""
    allow = np.ones((len(RHO), len(PSI)), bool); force, desc = [], []
    if not OVERRIDES.exists():
        return allow, force, ""
    o = pd.read_csv(OVERRIDES, sep="\t", dtype=str).fillna("")
    o = o[o["sample"] == s]
    if not len(o):
        return allow, force, ""
    r = o.iloc[0]
    if r.get("purity", ""):
        lo, hi = window(r.purity, 0.05); allow &= ((RHO >= lo - 1e-9) & (RHO <= hi + 1e-9))[:, None]
        desc.append(f"purity {lo:.2f}-{hi:.2f}")
    if r.get("ploidy", ""):         # psi = the grid ploidy; the reported ploidy is the segment-weighted one (close)
        lo, hi = window(r.ploidy, 0.25); allow &= ((PSI >= lo - 1e-9) & (PSI <= hi + 1e-9))[None, :]
        desc.append(f"ploidy {lo:.2f}-{hi:.2f}")
    if r.get("diploid", ""):
        cb = pd.read_csv(RES / "cytoband_b38.bed", sep="\t", header=None)
        cen = cb[cb[4] == "acen"].groupby(0).agg({1: "min", 2: "max"}); end = cb.groupby(0)[2].max()
        for x in [x.strip() for x in r.diploid.split(",") if x.strip()]:
            c = x if x.startswith("chr") else "chr" + x
            arm = c[-1] if c[-1] in "pq" else ""
            c = c.rstrip("pq")
            if c not in cen.index:
                raise SystemExit(f"{OVERRIDES}: unknown region {x!r} for {s}")
            a, b = {"p": (0, cen.loc[c, 1]), "q": (cen.loc[c, 2], end[c]), "": (0, end[c])}[arm]
            force.append((c, a, b))
        desc.append(f"diploid {r.diploid}")
    if not allow.any():
        raise SystemExit(f"{OVERRIDES}: purity/ploidy window for {s} leaves no grid point")
    return allow, force, "; ".join(desc)


def event_filter(S, base):
    """Merge adjacent passing non-focal calls of the same type (gain / loss / LOH / imbalance vs the background state)
    into events and fail every segment of an event whose summed evidence (depth + phased + all-sites) < EVENT_MIN.
    Focal calls keep their own, stricter test. In the dilutions this removed 9 of 11 false CBS-union events at
    TF >= 0.10 and kept 70 of 76 true ones (the 6 lost were weak pieces at TF 0.10)."""
    S = S.copy()
    S["_c"] = S.chrom.str[3:].astype(int); S = S.sort_values(["_c", "start"]).drop(columns="_c").reset_index(drop=True)
    bt = sum(base)

    def typ(r):
        a, b = r.nA, r.nB
        if r.sub_to:
            x = r.sub_to.split("@")[0].split(":"); a, b = int(x[0]), int(x[1])
        return "LOH" if b == 0 else ("gain" if a + b > bt else ("loss" if a + b < bt else "imbalance"))
    ev, cur, prev = [], None, None
    for i, r in S.iterrows():
        if not (r.call_pass and not r.focal):
            cur = None; continue
        ty = typ(r)
        if cur is not None and prev.chrom == r.chrom and cur[0] == ty:
            cur[1].append(i)
        else:
            cur = (ty, [i]); ev.append(cur)
        prev = r
    for ty, idx in ev:
        tot = float((S.loc[idx, "llr_depth"] + S.loc[idx, "llr_phased"] + S.loc[idx, "llr_allsites"]).sum())
        S.loc[idx, "event_llr"] = round(tot, 1)
        if tot < EVENT_MIN:
            S.loc[idx, "call_pass"] = False
            S.loc[idx, "filter_reason"] = "event_weak"
    return S


def local_maxima(T, n=3):
    out = []
    for i in range(3, T.shape[0] - 3):
        for j in range(3, T.shape[1] - 3):
            w = T[i - 3:i + 4, j - 3:j + 4]
            if T[i, j] >= w.max():
                out.append((T[i, j], i, j))
    i0, j0 = np.unravel_index(np.argmax(T), T.shape)
    out.append((T[i0, j0], i0, j0))
    out = sorted(o for o in set(out) if np.isfinite(o[0]))[::-1]
    return out[:n]


summary = []
for smp in sys.argv[1:] or sys.exit("usage: joint_fit.py <samples...>"):
    S_all, norm, phi, nu, src = build_segments(smp)
    allow, force, cons = constraints(smp)
    mid = (S_all.start + S_all.end) / 2
    S_all["force11"] = [any(c == fc and a <= m_ < b for fc, a, b in force) for c, m_ in zip(S_all.chrom, mid)]
    attrs0 = dict(S_all.attrs)
    S = S_all[~S_all.focal].reset_index(drop=True); S.attrs.update(attrs0)        # purity/ploidy fit: >= MIN_SEG only
    F = S_all[S_all.focal].reset_index(drop=True); F.attrs.update(attrs0)
    T, best = score(S)
    L = S.mb.to_numpy()
    if S.force11.any():                                   # free fit (nothing held at 1:1) for detection and dLL
        Sf = S.assign(force11=False); Sf.attrs.update(attrs0); llfree = float(score(Sf)[0].max())
    else:
        llfree = float(T.max())
    opts = local_maxima(np.where(allow, T, -np.inf))
    llmax = opts[0][0]
    (_, i0, j0) = opts[0]
    rho, psi = RHO[i0], PSI[j0]
    st = [STATES[b[i0, j0]] for b in best]
    S["nA"], S["nB"] = [max(x[0], x[1]) for x in st], [min(x[0], x[1]) for x in st]
    S["sub_to"] = [f"{max(x[2], x[3])}:{min(x[2], x[3])}@{x[4]:.1f}" if x[4] > 0 else "" for x in st]
    eff = [(x[0] + x[4] * (x[2] - x[0])) + (x[1] + x[4] * (x[3] - x[1])) for x in st]
    ploidy = float(np.sum(np.array(eff) * L) / L.sum())
    try:
        tf_i, pl_i = params(smp)[0], ichor_ploidy(smp)
        ii, jj = int(np.argmin(np.abs(RHO - tf_i))), int(np.argmin(np.abs(PSI - pl_i)))
        dll_i = round(float(T[ii, jj] - llmax), 1)
    except (FileNotFoundError, KeyError, ValueError):
        tf_i = pl_i = dll_i = np.nan
    alts, best_alt = [], np.nan
    for v, i, j in opts[1:]:
        stj = [STATES[b[i, j]] for b in best]
        pj = np.sum([(x[0] + x[4] * (x[2] - x[0]) + x[1] + x[4] * (x[3] - x[1])) * l for x, l in zip(stj, L)]) / L.sum()
        alts.append(f"rho={RHO[i]:.2f} ploidy={pj:.2f} dLL={v - llmax:.1f}")
    # best solution anywhere on the grid at a clearly different ploidy (e.g. genome doubling either way)
    far = np.abs(PSI - psi) >= 1.0
    if far.any():
        ia, ja = np.unravel_index(np.argmax(np.where(far[None, :], T, -np.inf)), T.shape)
        best_alt = float(T[ia, ja] - llmax)
        stj = [STATES[b[ia, ja]] for b in best]
        pa = np.sum([(x[0] + x[4] * (x[2] - x[0]) + x[1] + x[4] * (x[3] - x[1])) * l for x, l in zip(stj, L)]) / L.sum()
        alts.append(f"other-ploidy best: rho={RHO[ia]:.2f} ploidy={pa:.2f} dLL={best_alt:.1f}")
    base = S[S.sub_to == ""].groupby(["nA", "nB"]).mb.sum().idxmax()
    tumor_llr = float(llfree - null_ll(S))                # detection from the >= MIN_SEG segments, free fit
    attrs = dict(S.attrs)
    S = S.join(call_evidence(S, st, rho, psi, base)); S.attrs.update(attrs)      # join drops attrs
    if len(F):                                            # focal segments: best state at the fitted (rho, psi)
        _, bestF = score(F)
        stF = [STATES[b[i0, j0]] for b in bestF]
        F["nA"], F["nB"] = [max(x[0], x[1]) for x in stF], [min(x[0], x[1]) for x in stF]
        F["sub_to"] = [f"{max(x[2], x[3])}:{min(x[2], x[3])}@{x[4]:.1f}" if x[4] > 0 else "" for x in stF]
        F = F.join(call_evidence(F, stF, rho, psi, base))
        S = pd.concat([S, F], ignore_index=True)
        S["_c"] = S.chrom.str[3:].astype(int); S = S.sort_values(["_c", "start"]).drop(columns="_c").reset_index(drop=True)
        S.attrs.update(attrs)
    S = event_filter(S, base); S.attrs.update(attrs)
    out = S.drop(columns=["R_pred", "ll"]).assign(rho=rho)
    out.to_csv(OUT / f"{smp}.joint_segments.tsv", sep="\t", index=False)
    # LOH = no B allele in any tumor cell: clonal nB = 0, or a subclonal mix of two nB = 0 states (e.g. 1:0 with 2:0)
    sub_minor0 = np.array([u == "" or u.split("@")[0].split(":")[1] == "0" for u in S.sub_to])
    loh = S[(S.nB == 0) & sub_minor0 & S.call_pass]
    detected = tumor_llr > TUMOR_LLR_MIN
    # tumor share of cfDNA: rho counts tumor CELLS, ploidy-weighted here. Unlike rho it is the same for the
    # genome-doubled and the diploid explanation of the same data, and it is what dilution scales.
    ctdna = round(float(rho * ploidy / (rho * ploidy + 2 * (1 - rho))), 3)
    row = dict(sample=smp, tumor_llr=round(tumor_llr, 1), tumor_detected=bool(detected),
               TF_report=rho if detected else "not_detected", ctDNA_fraction=ctdna,
               ctDNA_report=ctdna if detected else "not_detected", ichor_TF=tf_i, ichor_ploidy=pl_i, joint_TF=rho, joint_ploidy=round(ploidy, 2),
               dLL_ichor=dll_i, dLL_other_ploidy=round(float(best_alt), 1), ambiguous=bool(best_alt > -10),
               loh_mb=round(float(loh.mb.sum()), 1),
               cnloh_mb=round(float(loh[loh.nA == int(round(ploidy))].mb.sum()), 1),   # copy-neutral vs sample ploidy
               calls=int(S.is_call.sum()), calls_pass=int(S.call_pass.sum()), segments=len(S),
               focal_calls_pass=int((S.call_pass & S.focal).sum()), dropout_segs=int((~S.use_phased).sum()), allsites_phi=round(phi, 2), nu=nu, theta0=S.attrs["theta0"],
               constraints=cons, dLL_constraint=round(llmax - llfree, 1) if cons else 0.0, seg_source=src, depth_source=S.attrs["depth_source"], alternatives="; ".join(alts))
    tmp1 = OUT / f".{smp}.{os.getpid()}.tsv"                # atomic: concurrent samples read these files
    pd.DataFrame([row]).to_csv(tmp1, sep="\t", index=False); os.replace(tmp1, OUT / f"{smp}.joint_summary.tsv")
    summary.append(row)
    print(f"\n=== {smp}: tumor LLR vs all-normal {tumor_llr:.1f} ({'detected' if detected else 'NOT detected'})")
    print(f"=== {smp}: joint TF {rho:.2f}, ploidy {ploidy:.2f}, ctDNA fraction {ctdna:.2f}  (ichorCNA {tf_i:.2f} / {pl_i:.2f}: dLL {dll_i})")
    print("   alternatives:", "; ".join(alts) or "none")
    if cons:
        print(f"   constraints: {cons}  (dLL vs free fit {llmax - llfree:.1f})")
    base = S.groupby(["nA", "nB"]).mb.sum().idxmax()          # most common state = background
    show = S[(S.nA != base[0]) | (S.nB != base[1]) | (S.sub_to != "") | (S.theta_hat > 0.03) | (S.D < 0.7)]
    for r in show.itertuples():
        print(f"   {r.chrom}:{r.start/1e6:6.1f}-{r.end/1e6:6.1f}  logR {r.r_obs:+.2f}  theta_hat {r.theta_hat:.3f}  "
              f"both-allele {r.R_obs:.2f}  hetD {r.D:.2f}{' (phased dropped)' if not r.use_phased else ''}  -> {r.nA}:{r.nB}{(' subclonal ' + r.sub_to) if r.sub_to else ''}"
              f"{'  [PASS]' if r.call_pass else ('  [filtered: ' + r.filter_reason + ']' if r.is_call else '')}")

    fig, ax = plt.subplots(figsize=(7, 5.2), facecolor=SURF); style(ax); ax.grid(False)
    dl = np.clip(T - llmax, -200, 0)
    im = ax.imshow(dl.T, origin="lower", aspect="auto", cmap="Blues", vmin=-200, vmax=0,
                   extent=[RHO[0], RHO[-1], PSI[0], PSI[-1]])
    cb = fig.colorbar(im, ax=ax, pad=0.02); cb.set_label("log-likelihood vs best (clipped at -200)", color=INK2, fontsize=9)
    cb.ax.tick_params(labelsize=8, colors=INK2); cb.outline.set_visible(False)
    if np.isfinite(tf_i):
        ax.scatter([tf_i], [pl_i], s=90, facecolor="none", edgecolor="#eb6834", lw=2.2, zorder=4)
        ax.annotate(f"ichorCNA {tf_i:.2f}/{pl_i:.2f}", (tf_i, pl_i), xytext=(8, 8) if tf_i < 0.6 else (-10, 8),
                textcoords="offset points", ha="left" if tf_i < 0.6 else "right", fontsize=9, color=INK, fontweight="bold")
    ax.scatter([rho], [psi], s=110, marker="X", color=INK, edgecolor=SURF, lw=1.2, zorder=5)
    ax.annotate(f"joint fit {rho:.2f}/{ploidy:.2f}", (rho, psi), xytext=(-10, -16) if rho > 0.6 else (8, -16),
                textcoords="offset points", ha="right" if rho > 0.6 else "left", fontsize=9, color=INK, fontweight="bold")
    ax.set_xlabel("tumor fraction (rho)", fontsize=9, color=INK2); ax.set_ylabel("ploidy grid (psi)", fontsize=9, color=INK2)
    ax.set_title(f"{smp}: joint fit (depth + phased BAF + all-sites allele rate)", fontsize=10.5, color=INK, loc="left",
                 fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / f"{smp}.joint_sunrise.png", dpi=130, facecolor=SURF); plt.close(fig)

# merged table rebuilt from the per-sample files; written via a temp file so concurrent runs never leave it half-written
allsum = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(OUT.glob("*.joint_summary.tsv"))])
tmp = OUT / f".joint_summary.{os.getpid()}.tsv"
allsum.to_csv(tmp, sep="\t", index=False); os.replace(tmp, OUT / "joint_summary.tsv")
