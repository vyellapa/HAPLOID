"""Read per-sample fit summaries and serialize updates of their merged table."""
import fcntl
import os
from pathlib import Path

import pandas as pd


def read_summaries(folder):
    """Prefer complete per-sample files; accept a merged-only table for older runs."""
    folder = Path(folder)
    files = sorted(folder.glob("*.joint_summary.tsv"))
    if not files:
        merged = folder / "joint_summary.tsv"
        files = [merged] if merged.exists() else []
    if not files:
        return pd.DataFrame(columns=["sample"])
    return pd.concat([pd.read_csv(f, sep="\t", dtype={"sample": str}) for f in files], ignore_index=True)


def merge_summaries(folder):
    """Lock before taking the snapshot so an older merge cannot replace a newer one."""
    folder = Path(folder)
    with (folder / ".joint_summary.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        table = read_summaries(folder)
        tmp = folder / f".joint_summary.{os.getpid()}.tsv"
        try:
            table.to_csv(tmp, sep="\t", index=False)
            os.replace(tmp, folder / "joint_summary.tsv")
        finally:
            tmp.unlink(missing_ok=True)
