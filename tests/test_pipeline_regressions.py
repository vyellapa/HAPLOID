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

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import summary_tables

HAS_ASCAT = shutil.which("Rscript") is not None and subprocess.run(
    ["Rscript", "-e", "suppressMessages({library(ASCAT); library(DNAcopy)})"], capture_output=True).returncode == 0
STEPS = ["markdup", "gl", "phase", "hap", "depth", "ascat", "joint", "plot", "cleanup"]


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

    def test_detected_samples_report_cbs_union_with_its_interval(self):
        self.write_fit("arms")
        self.write_fit("cbs_union", joint_TF=0.8, ctDNA_fraction=0.8, joint_ploidy=2.1,
                       TF_lo=0.79, TF_hi=0.81, ctDNA_lo=0.79, ctDNA_hi=0.81)
        row = self.cohort().loc["S"]
        self.assertEqual((row.TF, row.TF_lo, row.TF_hi, row.ploidy), (0.8, 0.79, 0.81, 2.1))
        self.assertEqual((row.ctDNA_fraction, row.ctDNA_lo, row.ctDNA_hi), (0.8, 0.79, 0.81))
        self.assertEqual(row.tf_basis, "cbs_union")

    def test_low_tf_samples_report_arms_with_its_interval(self):
        self.write_fit("arms", joint_TF=0.06, ctDNA_fraction=0.06, TF_lo=0.05, TF_hi=0.09, ctDNA_lo=0.05, ctDNA_hi=0.09)
        self.write_fit("cbs_union", joint_TF=0.9, TF_lo=0.85, TF_hi=0.95)
        row = self.cohort().loc["S"]
        self.assertEqual((row.TF, row.TF_lo, row.TF_hi), (0.06, 0.05, 0.09))
        self.assertEqual(row.tf_basis, "arms")

    def test_cohort_reads_samples_missing_from_merged_table(self):
        self.write_fit("arms", sample="001")
        self.write_fit("arms", sample="002")
        folder = self.root / "work/joint_modes/arms"
        shutil.copy(folder / "001.joint_summary.tsv", folder / "joint_summary.tsv")
        self.assertEqual(set(self.cohort().index), {"001", "002"})

    def test_unavailable_intervals_and_negative_samples_stay_empty(self):
        self.write_fit("arms", TF_lo=0.38, TF_hi=0.42)              # not the reported fit
        self.write_fit("cbs_union", TF_lo=float("nan"), TF_hi=float("nan"))
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
        self.assertTrue(table.tf_basis.str.contains("depth HMM").all())
        self.write_overrides("fit_overrides_ichor_ploidy.tsv", [])
        table = self.cohort(env=env)
        self.assertTrue((table.TF == 0.4).all())
        self.assertTrue((table.tf_basis == "cbs_union").all())

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
        self.assertEqual(row.tf_basis, "cbs_union")
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

    def executable(self, path, text):
        path = self.root / path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        path.chmod(0o755)
        return path

    def stub_pipeline(self):
        """Tools, reference and resources run_cohort.sh checks for; run_sample.sh logs the steps it would run."""
        sh = "#!/bin/sh\n"
        self.executable("bin/samtools", sh + 'for c in $(seq 1 22); do printf "@SQ\\tSN:chr%s\\tLN:%s\\n" $c $c; done\n')
        self.executable("bin/bcftools", sh)
        self.executable("bin/Rscript", sh)
        self.executable("bin/python", sh)
        self.executable("glimpse/GLIMPSE2_phase", sh + "echo --bam-file\n")
        (self.root / "ref.fa").write_text(">chr1\nN\n")
        (self.root / "ref.fa.fai").write_text("".join(f"chr{c}\t{c}\n" for c in range(1, 23)))
        for f in ["hg38-blacklist.v2.bed", "cytoband_b38.bed", "exclusion/sites_hg38.bed.gz", "exclusion/bins_500kb.tsv"] + [
                f"{d}/chr{c}.{x}" for c in range(1, 23)
                for d, x in [("sites_hg38", "sites.vcf.gz"), ("sites_hg38", "sites.tsv.gz"), ("maps_b38", "b38.gmap.gz")]]:
            (self.root / "resources" / f).parent.mkdir(parents=True, exist_ok=True)
            (self.root / "resources" / f).write_text("placeholder\n")
        for c in range(1, 23):
            self.executable(f"glimpse_ref/.done.chr{c}", "")
        self.executable("scripts/run_sample.sh",
                        '#!/usr/bin/env bash\nROOT=$(cd "$(dirname "$0")/.." && pwd)\n'
                        'printf "%s\\n" "$2" >> "$ROOT/sample_bams.txt"\n'
                        f'for step in {" ".join(STEPS)}; do\n'
                        '  [ -e "$ROOT/work/done/$1.$step" ] || { echo $step >> "$ROOT/steps.txt"; touch "$ROOT/work/done/$1.$step"; }\n'
                        'done\n')
        self.executable("scripts/run_ichor.sh", '#!/usr/bin/env bash\n'
                        'for f in correctedDepth params seg; do echo computed > "$3/$1_Y.$f.txt"; done\n')
        self.executable("scripts/resolve_discordant.sh", sh)
        return dict(PATH=f"{self.root / 'bin'}{os.pathsep}{os.environ['PATH']}", PY=str(self.root / "bin/python"),
                    REF=str(self.root / "ref.fa"), GLIMPSE_RUNNER="native", GLIMPSE_BIN=str(self.root / "glimpse"),
                    HAPLOID_GLIMPSE_REF=str(self.root / "glimpse_ref"), HAPLOID_CONTROLS="",
                    HAPLOID_RES=str(self.root / "resources"), HAPLOID_FINAL=str(self.root / "work/final"))

    def bam(self, path):
        bam = self.root / path
        bam.parent.mkdir(parents=True, exist_ok=True)
        bam.write_bytes(b"synthetic BAM placeholder")
        Path(f"{bam}.bai").touch()
        return str(bam)

    def depth_results(self, folder, name):
        folder = self.root / folder
        folder.mkdir(parents=True, exist_ok=True)
        (folder / f"{name}.correctedDepth.txt").write_text("chr\tstart\tend\tlog2_TNratio_corrected\nchr1\t1\t500000\t0\n")
        for f in ["params", "seg"]:
            (folder / f"{name}.{f}.txt").write_text(f"{name}\n")
        return str(folder)

    def test_new_depth_results_rerun_the_fit(self):
        env = self.stub_pipeline()
        bam, supplied = self.bam("bams/S.bam"), self.depth_results("supplied", "S_Y")

        def steps_run(*args):
            (self.root / "steps.txt").write_text("")
            result = self.run_script("run_bam.sh", "--bam", bam, *args, env=env)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            return (self.root / "steps.txt").read_text().split()
        self.assertEqual(steps_run("--ichor", supplied), STEPS)
        self.assertEqual(steps_run("--ichor", supplied), [])
        self.assertEqual(steps_run(), STEPS[4:])                  # computed here instead of supplied
        self.assertEqual(steps_run(), [])
        self.assertEqual(steps_run("--ichor", supplied), STEPS[4:])

    def test_shared_depth_folder_matches_whole_sample_names(self):
        env = self.stub_pipeline()
        for name in ["S1_Y", "S10_Y"]:
            shared = self.depth_results("shared", name)
        result = self.run_script("run_bam.sh", "--bam", self.bam("bams/S1.bam"), "--ichor", shared, env=env)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual((self.root / "work/inputs/S1_500k/S1_Y.params.txt").read_text(), "S1_Y\n")

    def test_cohort_passes_bam_paths_literally(self):
        env = self.stub_pipeline()
        bam = self.bam("bams/a b $(touch $INJECTED);c/S.bam")
        self.depth_results(Path(bam).parent, "S_Y")
        injected = self.root / "injected"
        result = self.run_script("run_cohort.sh", bam, env={**env, "BAM_DIR": str(self.root / "bams"),
                                                            "INJECTED": str(injected)})
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse(injected.exists())
        self.assertEqual((self.root / "sample_bams.txt").read_text().splitlines(), [bam])
        self.assertTrue((self.root / "work/done/S.cleanup").exists())

    def test_joint_fit_without_phased_segments_still_fits(self):
        self.executable("bin/bcftools", '#!/bin/sh\nexec cat "$4"\n')      # query -f <format> <file>: inputs kept as text
        rng = np.random.default_rng(7)
        chroms = {"chr1": (248_956_422, 123_400_000), "chr2": (242_193_529, 93_900_000)}   # length, centromere
        res, work = self.root / "resources", self.root / "work"
        (res / "sites_hg38").mkdir(parents=True)
        (res / "cytoband_b38.bed").write_text("".join(
            f"{c}\t0\t{m - 10**6}\tp\tgneg\n{c}\t{m - 10**6}\t{m + 10**6}\tcen\tacen\n{c}\t{m + 10**6}\t{n}\tq\tgneg\n"
            for c, (n, m) in chroms.items()))
        (res / "hg38-blacklist.v2.bed").write_text("chr1\t0\t10000\tHigh Signal Region\n")
        for d in ["cnv", "gl", "phased"]:
            (work / d).mkdir(parents=True)
        bins = []
        for c, (n, _) in chroms.items():
            start = np.arange(0, n - 500_000, 500_000)
            bins.append(pd.DataFrame(dict(chrom=c, start=start + 1, pos=start + 250_000,
                                          logR=rng.normal(0, 0.05, len(start)), logR_source="ichorCNA")))
            pos = np.arange(10**6, n - 10**6, 20_000)
            af, reads = rng.uniform(0.1, 0.9, len(pos)), rng.integers(2, 6, len(pos))
            het = rng.random(len(pos)) < 2 * af * (1 - af)
            alt = np.where(het, rng.binomial(reads, 0.5), np.where(rng.random(len(pos)) < af, reads, 0))
            pd.DataFrame(dict(pos=pos, ad=[f"{r - a},{a}" for r, a in zip(reads, alt)])).to_csv(
                work / f"gl/S.{c}.gl.bcf", sep="\t", header=False, index=False)
            pd.DataFrame(dict(pos=pos, af=af.round(3))).to_csv(
                res / f"sites_hg38/{c}.sites.vcf.gz", sep="\t", header=False, index=False, compression=None)
            pd.DataFrame(dict(chrom=c, pos=pos)).to_csv(res / f"sites_hg38/{c}.sites.tsv.gz", sep="\t", header=False, index=False)
        pd.concat(bins).to_csv(work / "cnv/S.bins.tsv", sep="\t", index=False)
        hets = np.arange(50) * 10**6 + 2 * 10**6                  # 50 phased hets: far below 200 per arm
        pd.DataFrame(dict(chrom="chr1", pos=hets, h1=rng.integers(0, 3, 50), h2=rng.integers(0, 3, 50), cm=hets / 1e6)).to_csv(
            work / "phased/S.chr1.hapcounts.tsv.gz", sep="\t", index=False)
        out = work / "joint_modes/arms"
        result = self.run_script("joint_fit.py", "S", env=dict(
            PATH=f"{self.root / 'bin'}{os.pathsep}{os.environ['PATH']}", HAPLOID_SEG_MODE="arms", HAPLOID_JOINT_OUT=str(out), HAPLOID_RES=str(res),
            HAPLOID_CONTROLS=str(res), ICHOR_PARAMS=str(self.root / "missing.params.txt")))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        row = pd.read_csv(out / "S.joint_summary.tsv", sep="\t").iloc[0]
        self.assertEqual((row.segments, row.dropout_segs), (4, 4))
        self.assertTrue(np.isfinite([row.joint_TF, row.joint_ploidy, row.allsites_phi]).all())

    @unittest.skipUnless(HAS_ASCAT, "needs R with ASCAT and DNAcopy")
    def test_failed_ascat_refit_drops_earlier_segments(self):
        rng = np.random.default_rng(1)
        (self.root / "work/cnv").mkdir(parents=True)
        pd.DataFrame([(f"chr{c}", i * 500_000 + 1, i * 500_000 + 250_000, rng.normal(0, 0.05), 0.5 + rng.normal(0, 0.01))
                      for c in range(1, 23) for i in range(30)], columns=["chrom", "start", "pos", "logR", "baf"]).to_csv(
            self.root / "work/cnv/S.bins.tsv", sep="\t", index=False)
        fail = self.root / "fail.Rprofile"
        fail.write_text('setHook(packageEvent("ASCAT", "onLoad"), function(...)\n'
                        '  assignInNamespace("ascat.runAscat", function(...) stop("simulated failure"), "ASCAT"))\n')
        out = self.root / "work/cnv/S_ascat"
        for profile, fitted in [({}, True), ({"R_PROFILE_USER": str(fail)}, False)]:
            result = subprocess.run(["Rscript", str(self.root / "scripts/run_ascat.R"), "S", "XX", "logR", "70"],
                                    env={**os.environ, **profile}, capture_output=True, text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((out / "S.ascat_segments.tsv").exists(), fitted, result.stderr)
            self.assertTrue((out / "S.aspcf_segments.tsv").exists())

    def test_cohort_bias_keeps_leading_zero_sample_names(self):
        (self.root / "work/cnv").mkdir(parents=True)
        for s in ["001", "002"]:
            self.write_fit("arms", sample=s, tumor_llr=0.)
            pd.DataFrame(dict(chrom="chr1", start=np.arange(20) * 500_000 + 1, logR=np.linspace(-0.1, 0.1, 20))).to_csv(
                self.root / f"work/cnv/{s}.bins.tsv", sep="\t", index=False)
        summary_tables.merge_summaries(self.root / "work/joint_modes/arms")
        result = self.run_script("cohort_bias.py", "--min-n", "2")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((pd.read_csv(self.root / "resources/cohort_bias_500kb.tsv", sep="\t").n == 2).all())


if __name__ == "__main__":
    unittest.main()
