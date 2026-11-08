"""Shared phased-allele helpers: numbat-style 2-state allele HMM over GLIMPSE2-phased hets.

hidden state = which phased haplotype carries the major allele (up/down); switch probability between consecutive
hets p_s = (1 - exp(-2 nu d_cM)) / 2 (numbat R/utils.R:463-475); H1 reads ~ Binomial(n, 0.5 +/- theta)
(numbat R/hmm.R:14-61). The forward likelihood is evaluated on a theta grid, so phase is marginalised.
"""
import numpy as np, pandas as pd
from allelic_lib import ROOT, in_mask

THETA = np.round(np.arange(0.0, 0.4951, 0.005), 3)
NU_GRID = (0.25, 0.5, 1, 2, 4, 8, 16)
LLR_IMBALANCED = 10.0     # 2*LLR threshold for calling a region imbalanced


def forward_ll(h1, h2, dcm, nu):
    """Log-likelihood over the THETA grid for one run of hets (sorted by position)."""
    ps = np.clip((1 - np.exp(-2 * nu * dcm)) / 2, 1e-10, 0.5)
    lu, ld = np.log1p(2 * THETA), np.log1p(-2 * THETA)
    a_up = np.full(len(THETA), 0.5); a_dn = np.full(len(THETA), 0.5); ll = np.zeros(len(THETA))
    for i in range(len(h1)):
        if i:
            p = ps[i]
            a_up, a_dn = a_up * (1 - p) + a_dn * p, a_up * p + a_dn * (1 - p)
        eu = np.exp(h1[i] * lu + h2[i] * ld); ed = np.exp(h1[i] * ld + h2[i] * lu)
        a_up, a_dn = a_up * eu, a_dn * ed
        sc = a_up + a_dn; ll += np.log(sc); a_up /= sc; a_dn /= sc
    return ll


def theta_fit(x, nu):
    """theta_hat, curvature SE, 2*LLR vs theta=0 for a DataFrame of hets (needs h1, h2, cm)."""
    x = x.sort_values("pos")
    dcm = np.diff(x.cm.to_numpy(), prepend=x.cm.iat[0])
    ll = forward_ll(x.h1.to_numpy(), x.h2.to_numpy(), dcm, nu)
    k = int(np.argmax(ll)); se = np.nan
    if 0 < k < len(THETA) - 1:
        curv = (ll[k - 1] - 2 * ll[k] + ll[k + 1]) / (THETA[1] - THETA[0]) ** 2
        se = 1 / np.sqrt(-curv) if curv < 0 else np.nan
    return float(THETA[k]), se, float(2 * (ll[k] - ll[0]))


def load_hapcounts(s):
    """Phased het counts, re-masked on load (mask files may be newer than the counts)."""
    out = []
    for f in sorted((ROOT / "work/phased").glob(f"{s}.chr*.hapcounts.tsv.gz")):
        d = pd.read_csv(f, sep="\t")
        if len(d):
            out.append(d[~in_mask(d.chrom.iat[0], d.pos.to_numpy())])
    return pd.concat(out, ignore_index=True)


def fit_nu(d, regions, min_llr=50):
    """Grid-MLE of nu on regions [(chrom, start, end)] that are strongly imbalanced at nu = 1."""
    runs = []
    for c, a, b in regions:
        x = d[(d.chrom == c) & d.pos.between(a, b)].sort_values("pos")
        if len(x) < 200:
            continue
        dcm = np.diff(x.cm.to_numpy(), prepend=x.cm.iat[0])
        ll = forward_ll(x.h1.to_numpy(), x.h2.to_numpy(), dcm, 1.0)
        if 2 * (ll.max() - ll[0]) > min_llr:
            runs.append((x.h1.to_numpy(), x.h2.to_numpy(), dcm))
    if not runs:
        return 1.0, {}
    tot = {nu: sum(forward_ll(h1, h2, dcm, nu).max() for h1, h2, dcm in runs) for nu in NU_GRID}
    return max(tot, key=tot.get), tot
