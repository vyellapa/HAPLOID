"""Fragment counts per 500 kb bin, stratified by fragment-length class x fragment GC (CNACS-style input).

For each sample: read-1 of proper pairs, MAPQ >= 20, not duplicate / secondary / supplementary / QC-fail,
fragment length 50-600 bp. Fragment GC is computed over the fragment's own reference span (CNACS bait_gc.pl
idea) from a per-chromosome GC prefix sum. Also builds (once, sample-independent) the 'possible fragments'
table: GC of every reference window of each class's representative length, sampled every STEP bp, so
expected depth per bin can be predicted from GC (CNACS correct_gc.sh / predicted_depth).

Usage: fragment_counts.py <sample> [<sample> ...]     (BAMs: work/bam/<sample>.auto.md.bam)
Output: work/depth/<sample>.frag_counts.npz (obs[bin, len_class, gc_pct]), work/depth/possible.npz
"""
import os, sys
from pathlib import Path
import numpy as np, pysam

ROOT = Path(__file__).resolve().parents[1]
REF = Path(os.environ.get("REF", ROOT / "resources/hg38.fa"))
BIN = 500_000
LEN_EDGES = np.array([50, 130, 160, 190, 601])          # classes: <130 (tumor-enriched), 130-160, 160-190, >=190
LEN_REP = np.array([115, 145, 175, 215])                # representative length per class (possible table)
STEP = 50
CHROMS = [f"chr{i}" for i in range(1, 23)]
OUT = ROOT / "work/depth"; OUT.mkdir(parents=True, exist_ok=True)

fa = pysam.FastaFile(str(REF))
nbins = {c: fa.get_reference_length(c) // BIN + 1 for c in CHROMS}


LUT_GC = np.zeros(256, np.int8); LUT_GC[[ord(x) for x in "GCgc"]] = 1
LUT_N = np.zeros(256, np.int8); LUT_N[[ord(x) for x in "Nn"]] = 1


def gc_prefix(c):
    """Prefix sums of G/C and N counts along a chromosome (int32 lookup tables keep peak memory ~3 GB)."""
    s = np.frombuffer(fa.fetch(c).encode("ascii"), dtype=np.uint8)
    out = []
    for lut in (LUT_GC, LUT_N):
        p = np.empty(len(s) + 1, np.int32); p[0] = 0
        np.cumsum(lut[s], dtype=np.int32, out=p[1:])
        out.append(p)
    return out[0], out[1]


def possible_table():
    f = OUT / "possible.npz"
    if f.exists():
        return
    tabs = {}
    for c in CHROMS:
        cg, cn = gc_prefix(c)
        L = len(cg) - 1
        t = np.zeros((nbins[c], len(LEN_REP), 101), np.int64)
        for k, ln in enumerate(LEN_REP):
            st = np.arange(0, L - ln, STEP)
            ok = (cn[st + ln] - cn[st]) == 0
            st = st[ok]
            g = np.rint(100 * (cg[st + ln] - cg[st]) / ln).astype(int)
            b = (st + ln // 2) // BIN
            np.add.at(t, (b, k, g), 1)
        tabs[c] = t
        print(f"possible: {c}", flush=True)
    np.savez_compressed(f, **tabs)


def sample_table(s):
    bam = pysam.AlignmentFile(str(ROOT / f"work/bam/{s}.auto.md.bam"))
    tabs, n_tot = {}, 0
    for c in CHROMS:
        cg, cn = gc_prefix(c)
        starts, lens = [], []
        for r in bam.fetch(c):
            if (not r.is_read1 or not r.is_proper_pair or r.is_duplicate or r.is_secondary or r.is_supplementary
                    or r.is_qcfail or r.mapping_quality < 20):
                continue
            tl = abs(r.template_length)
            if 50 <= tl <= 600:
                starts.append(min(r.reference_start, r.next_reference_start)); lens.append(tl)
        st, ln = np.array(starts, np.int64), np.array(lens, np.int64)
        en = np.minimum(st + ln, len(cg) - 1)
        g = np.rint(100 * (cg[en] - cg[st]) / np.maximum(en - st, 1)).astype(int).clip(0, 100)
        k = np.searchsorted(LEN_EDGES, ln, side="right") - 1
        b = (st + ln // 2) // BIN
        t = np.zeros((nbins[c], len(LEN_REP), 101), np.int64)
        np.add.at(t, (b, k, g), 1)
        tabs[c] = t; n_tot += len(st)
        print(f"{s} {c}: {len(st)} fragments", flush=True)
    np.savez_compressed(OUT / f"{s}.frag_counts.npz", **tabs)
    print(f"{s}: total fragments {n_tot}")


if __name__ == "__main__":
    if sys.argv[1] == "--possible":
        possible_table()
    else:
        for smp in sys.argv[1:]:
            sample_table(smp)
