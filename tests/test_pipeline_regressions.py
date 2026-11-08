"""Pipeline regressions using synthetic summaries and stubbed sequencing tools.

Run with: python3 -m unittest discover -s tests -v
"""
import multiprocessing
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import summary_tables


def delayed_merge(folder, snapshot, release):
    read = summary_tables.read_summaries

    def paused_read(path):
        table = read(path)
        snapshot.set()
        if not release.wait(10):
            raise TimeoutError("merge test was not released")
        return table

    summary_tables.read_summaries = paused_read
    summary_tables.merge_summaries(folder)


def next_merge(folder, started, done):
    started.set()
    summary_tables.merge_summaries(folder)
    done.set()


class PipelineRegressions(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="haploid-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        shutil.copytree(SCRIPTS, self.root / "scripts", ignore=shutil.ignore_patterns("__pycache__"))

    def run_script(self, name, *args, env=None):
        command = ["bash"] if name.endswith(".sh") else [sys.executable]
        script_env = {**os.environ, "PY": sys.executable}
        script_env.pop("HAPLOID_OVERRIDES", None)  # Tests select their own constraint tables.
        return subprocess.run(
            command + [str(self.root / "scripts" / name), *args],
            env={**script_env, **(env or {})},
            capture_output=True, text=True, timeout=20,
        )

    def write_fit(self, mode, sample="S", **values):
        row = dict(sample=sample, tumor_detected=True, tumor_llr=100., joint_TF=0.4, joint_ploidy=2.,
                   ctDNA_fraction=0.4, ambiguous=False, dLL_other_ploidy=-20., ichor_TF=0.4,
                   ichor_ploidy=2., cnloh_mb=0., constraints="", depth_source="ichorCNA",
                   dLL_constraint=0., TF_lo=0.38, TF_hi=0.42, ctDNA_lo=0.38, ctDNA_hi=0.42)
        row.update(values)
        folder = self.root / "work/joint_modes" / mode
        folder.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([row]).to_csv(folder / f"{sample}.joint_summary.tsv", sep="\t", index=False)
        pd.DataFrame([dict(chrom="chr1", start=1, end=10000001, nA=1, nB=1, sub_to="",
                           mb=10., call_pass=False, focal=False)]).to_csv(
            folder / f"{sample}.joint_segments.tsv", sep="\t", index=False)

    def write_overrides(self, name, samples):
        pd.DataFrame([dict(sample=s, purity="", ploidy="1.7-2.3", diploid="") for s in samples],
                     columns=["sample", "purity", "ploidy", "diploid"]).to_csv(
            self.root / "work" / name, sep="\t", index=False)

    def cohort(self, env=None):
        result = self.run_script("cohort_table.py", env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        return pd.read_csv(self.root / "work/cohort_summary.tsv", sep="\t", dtype={"sample": str}).set_index("sample")

    def test_intervals_follow_reported_estimates(self):
        self.write_fit("arms")
        self.write_fit("cbs_union", joint_TF=0.8, ctDNA_fraction=0.8,
                       TF_lo=0.79, TF_hi=0.81, ctDNA_lo=0.79, ctDNA_hi=0.81)
        row = self.cohort().loc["S"]
        self.assertEqual((row.TF, row.TF_lo, row.TF_hi), (0.4, 0.38, 0.42))
        self.assertEqual((row.ctDNA_fraction, row.ctDNA_lo, row.ctDNA_hi), (0.4, 0.38, 0.42))
        self.assertEqual(row.tf_basis, "arms")

    def test_cohort_reads_samples_missing_from_merged_table(self):
        self.write_fit("arms", sample="001")
        self.write_fit("arms", sample="002")
        folder = self.root / "work/joint_modes/arms"
        shutil.copy(folder / "001.joint_summary.tsv", folder / "joint_summary.tsv")
        self.assertEqual(set(self.cohort().index), {"001", "002"})

    def test_unavailable_intervals_and_negative_samples_stay_empty(self):
        self.write_fit("arms", TF_lo=float("nan"), TF_hi=float("nan"))
        self.write_fit("cbs_union", TF_lo=0.79, TF_hi=0.81)
        self.write_fit("arms", sample="negative", tumor_detected=False)
        table = self.cohort()
        self.assertTrue(pd.isna(table.loc["S"].TF_lo))
        self.assertTrue(pd.isna(table.loc["S"].TF_uncertain))
        self.assertEqual(table.loc["negative"].TF, "not_detected")
        self.assertTrue(pd.isna(table.loc["negative"].TF_lo))

    def test_clonal_results_merge_with_numeric_sample_names(self):
        self.write_fit("arms", sample="123")
        pd.DataFrame([dict(sample="123", clonal_TF=0.3, clonal_ploidy=2., ctDNA_clonal=0.3,
                           clonal_frac_of_aberrant=1., runner_up_TF=0.7, d_cost_runner_up=2.)]).to_csv(
            self.root / "work/clonal_tf.tsv", sep="\t", index=False)
        self.assertEqual(self.cohort().loc["123"].clonal_TF, 0.3)

    def test_active_refits_use_their_own_estimates_and_intervals(self):
        self.write_fit("arms")
        self.write_fit("cbs_union")
        self.write_fit("cbs_union_ichorploidy", joint_TF=0.6, ctDNA_fraction=0.6,
                       TF_lo=0.58, TF_hi=0.62, ctDNA_lo=0.58, ctDNA_hi=0.62)
        self.write_overrides("fit_overrides_ichor_ploidy.tsv", ["S"])
        row = self.cohort().loc["S"]
        self.assertEqual((row.TF, row.TF_lo, row.TF_hi), (0.6, 0.58, 0.62))
        self.assertEqual((row.ctDNA_fraction, row.ctDNA_lo, row.ctDNA_hi), (0.6, 0.58, 0.62))
        self.write_fit("cbs_union_user", joint_TF=0.7, TF_lo=0.68, TF_hi=0.72)
        self.write_overrides("user_constraints.tsv", ["S"])
        self.assertEqual(self.cohort().loc["S"].TF, 0.7)
        self.write_overrides("user_constraints.tsv", [])
        self.assertEqual(self.cohort().loc["S"].TF, 0.6)
        self.write_overrides("fit_overrides_ichor_ploidy.tsv", [])
        self.assertEqual(self.cohort().loc["S"].TF, 0.4)

    def test_custom_override_file_selects_refit_and_intervals(self):
        self.write_fit("arms", sample="001")
        self.write_fit("cbs_union", sample="001")
        self.write_fit("cbs_union_user", sample="001", joint_TF=0.7, ctDNA_fraction=0.7,
                       TF_lo=0.68, TF_hi=0.72, ctDNA_lo=0.68, ctDNA_hi=0.72)
        self.write_overrides("batch_constraints.tsv", ["001"])
        env = {"HAPLOID_OVERRIDES": str(self.root / "work/batch_constraints.tsv")}
        row = self.cohort(env=env).loc["001"]
        self.assertEqual((row.TF, row.TF_lo, row.TF_hi), (0.7, 0.68, 0.72))
        self.assertEqual((row.ctDNA_fraction, row.ctDNA_lo, row.ctDNA_hi), (0.7, 0.68, 0.72))
        self.assertIn("your values", row.tf_basis)

    def test_custom_override_membership_controls_cached_user_refits(self):
        for sample in ["001", "002"]:
            self.write_fit("arms", sample=sample)
            self.write_fit("cbs_union", sample=sample)
            self.write_fit("cbs_union_user", sample=sample, joint_TF=0.7)
            self.write_fit("cbs_union_ichorploidy", sample=sample, joint_TF=0.6)
        self.write_overrides("user_constraints.tsv", ["001", "002"])
        self.write_overrides("fit_overrides_ichor_ploidy.tsv", ["001", "002"])
        self.write_overrides("batch_constraints.tsv", ["001"])
        env = {"HAPLOID_OVERRIDES": str(self.root / "work/batch_constraints.tsv")}
        table = self.cohort(env=env)
        self.assertEqual(table.loc["001"].TF, 0.7)
        self.assertEqual(table.loc["002"].TF, 0.6)
        self.write_overrides("batch_constraints.tsv", [])
        table = self.cohort(env=env)
        self.assertTrue((table.TF == 0.6).all())
        self.assertTrue(table.tf_basis.str.contains("ichorCNA").all())
        self.write_overrides("fit_overrides_ichor_ploidy.tsv", [])
        table = self.cohort(env=env)
        self.assertTrue((table.TF == 0.4).all())
        self.assertTrue((table.tf_basis == "arms").all())

    def test_exclusion_updates_cohort_without_refitting(self):
        self.write_fit("arms")
        self.write_fit("cbs_union", joint_ploidy=4.)
        self.write_fit("cbs_union_ichorploidy", joint_TF=0.6)
        self.write_overrides("fit_overrides_ichor_ploidy.tsv", ["S"])
        self.assertEqual(self.cohort().loc["S"].TF, 0.6)
        result = self.run_script("resolve_discordant.sh", "S")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("nothing to refit", result.stdout)
        row = pd.read_csv(self.root / "work/cohort_summary.tsv", sep="\t").iloc[0]
        self.assertEqual(row.TF, 0.4)
        self.assertEqual(row.tf_basis, "arms")
        self.assertTrue((self.root / "work/joint_modes/cbs_union_ichorploidy/S.joint_summary.tsv").exists())

    def test_manual_refits_use_their_own_override_table(self):
        self.write_fit("arms")
        self.write_fit("cbs_union", calls_pass=0)
        self.write_fit("cbs_union_user", joint_TF=0.7, calls_pass=0)
        self.write_overrides("batch_constraints.tsv", [])
        # Keep the prepared fit so this test exercises table selection without sequencing tools.
        (self.root / "scripts/refit_downstream.sh").write_text("#!/usr/bin/env bash\nexit 0\n")
        env = {"HAPLOID_OVERRIDES": str(self.root / "work/batch_constraints.tsv")}
        result = self.run_script("refit_with.sh", "S", "ploidy=2", env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        row = pd.read_csv(self.root / "work/cohort_summary.tsv", sep="\t").iloc[0]
        self.assertEqual(row.TF, 0.7)
        self.assertIn("your values", row.tf_basis)

    def test_refit_failures_stop_sample_and_propagate_across_batches(self):
        runner = self.root / "runner"
        runner.write_text(
            f"#!{sys.executable}\nimport os, sys\n"
            "with open(os.environ['CALLS'], 'a') as f:\n    f.write(' '.join(sys.argv[1:3]) + '\\n')\n"
            "if sys.argv[1] == os.environ['FAIL_STEP'] and sys.argv[2] == 'bad':\n    sys.exit(7)\n")
        runner.chmod(0o755)
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        shutil.copy(runner, bin_dir / "Rscript")
        calls = self.root / "calls.txt"
        steps = ["build_bins.py", "run_ascat.R", "joint_fit.py"]
        for fail_step in steps:
            for samples in [("bad", "good", "later"), ("good", "later", "bad")]:
                with self.subTest(step=fail_step, samples=samples):
                    calls.write_text("")
                    env = dict(PY=str(runner), PATH=str(bin_dir) + os.pathsep + os.environ["PATH"],
                               CALLS=str(calls), FAIL_STEP=fail_step, STEPS="bins seg joint",
                               MODES="arms cbs_union", JOBS="2")
                    result = self.run_script("refit_downstream.sh", *samples, env=env)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("bad refit FAILED", result.stderr)
                    self.assertNotIn("bad refit done", result.stdout)
                    for sample in ["good", "later"]:
                        self.assertIn(f"{sample} refit done", result.stdout)
                    bad_calls = [line.split()[0] for line in calls.read_text().splitlines() if line.endswith(" bad")]
                    self.assertEqual(bad_calls, steps[:steps.index(fail_step) + 1])
        result = self.run_script("refit_downstream.sh", "good", env=env)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_concurrent_merges_keep_all_completed_samples(self):
        folder = self.root / "summaries"
        folder.mkdir()
        pd.DataFrame([dict(sample="A")]).to_csv(folder / "A.joint_summary.tsv", sep="\t", index=False)
        context = multiprocessing.get_context("spawn")
        snapshot, release, started, done = [context.Event() for _ in range(4)]
        first = context.Process(target=delayed_merge, args=(folder, snapshot, release))
        second = context.Process(target=next_merge, args=(folder, started, done))
        first.start()
        try:
            self.assertTrue(snapshot.wait(5))
            pd.DataFrame([dict(sample="B")]).to_csv(folder / "B.joint_summary.tsv", sep="\t", index=False)
            second.start()
            self.assertTrue(started.wait(5))
            self.assertFalse(done.wait(0.2), "a second merge bypassed the snapshot lock")
        finally:
            release.set()
            for process in [first, second]:
                if process.pid is not None:
                    process.join(10)
                    if process.is_alive():
                        process.terminate()
                        process.join()
        self.assertEqual((first.exitcode, second.exitcode), (0, 0))
        table = pd.read_csv(folder / "joint_summary.tsv", sep="\t")
        self.assertEqual(set(table["sample"]), {"A", "B"})


if __name__ == "__main__":
    unittest.main()
