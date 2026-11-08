#!/usr/bin/env bash
# Refit one sample with values you know (tissue purity / ploidy, or chromosomes known to be unchanged).
#   refit_with.sh <sample> [purity=0.35 | purity=0.3-0.5] [ploidy=2 | ploidy=2.6-3.1] [diploid=chr2,chr8p] ...
#   refit_with.sh <sample> clear          # remove that sample's values (the free fit is used again)
# purity "0.35" = 0.30-0.40, ploidy "2" = 1.75-2.25, diploid = chromosomes / arms held at 1:1 (anchors copy number).
# The values are kept in work/user_constraints.tsv (one row per sample, later calls replace that sample's row); the joint
# fit (cbs_union) is rerun for that sample only into work/joint_modes/cbs_union_user/ - the free fit is kept - and
# cohort_table.py then reports TF / ploidy / ctDNA / events from it (user values take priority over the ichorCNA ploidy hold).
# Printed: free fit vs your values, and dLL_constraint (how much worse the fit is with them: ~ -10 or better = the data
# agree; hundreds below 0 = the data disagree, check the values). Tumor detection always uses the free fit.
# env: PY (python3), BAM_DIR (ichorCNA results: <BAM_DIR>/<s>_500k/)
set -euo pipefail
[ $# -ge 2 ] || { sed -n '2,10p' "$0"; exit 1; }
ROOT=$(cd "$(dirname "$0")/.." && pwd); SC=$ROOT/scripts; PY=${PY:-python3}
S=$1; shift
UC=$ROOT/work/user_constraints.tsv
UC=$UC S=$S "$PY" - "$@" << 'PYEOF'
import os, sys
from pathlib import Path
import pandas as pd
f, s = Path(os.environ["UC"]), os.environ["S"]
T = pd.read_csv(f, sep="\t", dtype=str).fillna("") if f.exists() else pd.DataFrame(columns=["sample", "purity", "ploidy", "diploid"])
T = T[T["sample"] != s]
if sys.argv[1:] != ["clear"]:
    row = {"sample": s, "purity": "", "ploidy": "", "diploid": ""}
    for a in sys.argv[1:]:
        k, _, v = a.partition("=")
        if k not in ("purity", "ploidy", "diploid") or not v:
            sys.exit(f"unknown value '{a}': use purity=..., ploidy=..., diploid=...")
        row[k] = v
    T = pd.concat([T, pd.DataFrame([row])])
    print(f"{s}: " + ", ".join(f"{k} {v}" for k, v in row.items() if k != "sample" and v))
else:
    print(f"{s}: values cleared")
T.to_csv(f, sep="\t", index=False)
PYEOF
if [ "$*" = "clear" ]; then
  rm -f "$ROOT/work/joint_modes/cbs_union_user/$S".joint_*
else
  HAPLOID_OVERRIDES=$UC STEPS=joint MODES=cbs_union OUT_SUFFIX=_user JOBS=1 BAM_DIR=${BAM_DIR:-} PY=$PY \
    bash "$SC/refit_downstream.sh" "$S"
  for d in cbs_union cbs_union_user; do
    awk -F'\t' -v d="$d" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==2{printf "  %-15s TF %s  ploidy %s  ctDNA %s  calls %s  dLL_constraint %s\n", (d=="cbs_union"?"free fit":"your values"), $h["joint_TF"], $h["joint_ploidy"], $h["ctDNA_fraction"], $h["calls_pass"], (d=="cbs_union"?"-":$h["dLL_constraint"])}' \
      "$ROOT/work/joint_modes/$d/$S.joint_summary.tsv"
  done
fi
cd "$SC" && HAPLOID_OVERRIDES=$UC "$PY" cohort_table.py > /dev/null && echo "  cohort table updated (tf_basis shows the source)"
