#!/usr/bin/env bash
# Samples where the joint fit (cbs_union) and ichorCNA disagree on ploidy: refit cbs_union with ploidy held at ichorCNA's
# (+/- PLOIDY_WIN), tumor fraction free. The joint fit drifts to ploidy ~4 (or below 2) with subclonal states where
# ichorCNA's depth-only ploidy is right, so ichorCNA anchors the ploidy and the joint fit
# still sets TF and the allele-specific states. Second pass, after cohort_table.py:
#   1. work/fit_overrides_ichor_ploidy.tsv: detected samples with |cbs_union ploidy - ichorCNA ploidy| > MIN_PLOIDY_DIFF,
#      minus QC-flagged ones (fit_flag, allele_qc low in work/cohort_summary.tsv) and the samples given as arguments
#   2. refit (joint fit only, no re-phasing) -> work/joint_modes/cbs_union_ichorploidy/ (the free fit is kept)
#   3. cohort_table.py: TF / ploidy / ctDNA / events of those samples from the constrained fit (tf_basis column);
#      dLL_constraint = what holding the ploidy costs the joint model
# Usage: resolve_discordant.sh [samples to leave out, e.g. suspected germline-CNV samples]
#   env: MIN_PLOIDY_DIFF (0.5), PLOIDY_WIN (0.3), JOBS (8), PY (python3), BAM_DIR (ichorCNA results: <BAM_DIR>/<s>_500k/)
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd); SC=$ROOT/scripts
PY=${PY:-python3}
OV=$ROOT/work/fit_overrides_ichor_ploidy.tsv
(cd "$SC" && "$PY" cohort_table.py > /dev/null)        # current QC flags for every finished sample
MIN_PLOIDY_DIFF=${MIN_PLOIDY_DIFF:-0.5} PLOIDY_WIN=${PLOIDY_WIN:-0.3} BAM_DIR=${BAM_DIR:-} OV=$OV ROOT=$ROOT \
"$PY" - "$@" << 'PYEOF'
import os, sys, glob
from pathlib import Path
import pandas as pd
ROOT = Path(os.environ["ROOT"]); leave = set(sys.argv[1:])
dmin, win, bam = float(os.environ["MIN_PLOIDY_DIFF"]), float(os.environ["PLOIDY_WIN"]), os.environ["BAM_DIR"]
S = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(str(ROOT / "work/joint_modes/cbs_union/*.joint_summary.tsv"))]) \
    .drop_duplicates("sample", keep="last").set_index("sample")
Q = pd.read_csv(ROOT / "work/cohort_summary.tsv", sep="\t").set_index("sample") if (ROOT / "work/cohort_summary.tsv").exists() else None
def ichor_ploidy(s):
    for f in ([Path(bam) / f"{s}_500k/{s}_Y.params.txt"] if bam else []) + [ROOT / f"ichor/{s}_Y.params.txt", ROOT / f"{s}_500k/{s}_Y.params.txt"]:
        if f.exists():
            kv = dict(l.split(":\t", 1) for l in f.read_text().splitlines() if ":\t" in l)
            return float(kv["Ploidy"])
    v = S.loc[s, "ichor_ploidy"] if "ichor_ploidy" in S else float("nan")
    return float(v)
rows = []
for s, r in S.iterrows():
    if str(r.tumor_detected) != "True" or s in leave:
        continue
    q = Q.loc[s] if Q is not None and s in Q.index else {}
    if (isinstance(q.get("fit_flag", ""), str) and q.get("fit_flag", "")) or q.get("allele_qc", "") == "low":
        continue
    p = ichor_ploidy(s)
    if p == p and abs(r.joint_ploidy - p) > dmin:
        rows.append(dict(sample=s, purity="", ploidy=f"{max(p - win, 1.0):.2f}-{p + win:.2f}", diploid="",
                         ichor_p=p, cbs_p=r.joint_ploidy))
P = pd.DataFrame(rows, columns=["sample", "purity", "ploidy", "diploid", "ichor_p", "cbs_p"])
for r in P.itertuples():
    print(f"  {r.sample}: cbs_union ploidy {r.cbs_p:.2f} -> held at ichorCNA {r.ichor_p:.2f} ({r.ploidy})")
P.drop(columns=["ichor_p", "cbs_p"]).to_csv(os.environ["OV"], sep="\t", index=False)
print(f"{len(P)} samples -> {os.environ['OV']}")
PYEOF
samples=$(tail -n +2 "$OV" | cut -f1)
[ -n "$samples" ] || { echo "nothing to refit"; exit 0; }
HAPLOID_OVERRIDES=$OV STEPS=joint MODES=cbs_union OUT_SUFFIX=_ichorploidy JOBS=${JOBS:-8} BAM_DIR=${BAM_DIR:-} PY=$PY \
  bash "$SC/refit_downstream.sh" $samples
cd "$SC" && "$PY" cohort_table.py
