"""Exercise container orchestration without downloading sequencing reference data."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


PROJECT = Path(__file__).resolve().parents[1]


class ContainerEntrypointTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="haploid-container-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "container").mkdir()
        shutil.copy(PROJECT / "container/entrypoint.sh", self.root / "container")
        scripts = self.root / "scripts"
        scripts.mkdir()
        (scripts / "run_cohort.sh").write_text(
            '#!/usr/bin/env bash\nset -euo pipefail\n'
            'printf "%s\\n" "$@" > "$ROOT/cohort_bams.txt"\n'
            'exit_code=${PIPELINE_EXIT:-0}\n'
            '[ "$exit_code" -eq 0 ] || exit "$exit_code"\n'
            'mkdir -p "$ROOT/work/done"\n'
            'for bam in "$@"; do\n'
            '  sample=$(basename "$bam" .bam)\n'
            '  if [ "$sample" = "${FAIL_SAMPLE:-}" ]; then\n'
            '    echo phase > "$ROOT/work/done/$sample.FAILED"\n'
            '  else\n'
            '    touch "$ROOT/work/done/$sample.cleanup"\n'
            '  fi\n'
            'done\n'
        )
        for name in ("clonal_tf.py", "cohort_table.py"):
            (scripts / name).write_text(
                'import os\nfrom pathlib import Path\n'
                'root = Path(__file__).resolve().parents[1]\n'
                'if os.environ.get("REPORT_FAIL"):\n    raise SystemExit(9)\n'
                'with (root / "reports.txt").open("a") as out:\n'
                '    out.write(Path(__file__).name + "\\n")\n'
            )

    def run_container(self, *args, **env):
        return subprocess.run(
            ["bash", str(self.root / "container/entrypoint.sh"), *args],
            env={**os.environ, "ROOT": str(self.root), **env},
            capture_output=True, text=True, timeout=20,
        )

    def bam(self, sample):
        path = self.root / f"bams/{sample}_500k/{sample}.bam"
        path.parent.mkdir(parents=True)
        path.write_bytes(b"synthetic BAM placeholder")
        return str(path)

    def test_help_works_without_data(self):
        result = self.run_container()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("cohort", result.stdout)
        self.assertFalse((self.root / "cohort_bams.txt").exists())

    def test_sample_includes_shared_preparation_and_reports(self):
        bam = self.bam("001")
        result = self.run_container("sample", bam)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.root / "cohort_bams.txt").read_text().splitlines(), [bam])
        self.assertTrue((self.root / "work/done/001.cleanup").exists())
        self.assertEqual((self.root / "reports.txt").read_text().splitlines(),
                         ["clonal_tf.py", "cohort_table.py"])

    def test_cohort_discovers_bams_and_keeps_numeric_sample_names(self):
        bams = [self.bam("001"), self.bam("002")]
        result = self.run_container("cohort", BAM_DIR=str(self.root / "bams"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.root / "cohort_bams.txt").read_text().splitlines(), bams)

    def test_missing_bams_stop_before_pipeline(self):
        for args in [("sample", str(self.root / "missing.bam")), ("cohort",)]:
            with self.subTest(args=args):
                result = self.run_container(*args, BAM_DIR=str(self.root / "missing"))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("ERROR:", result.stderr)
                self.assertFalse((self.root / "cohort_bams.txt").exists())

    def test_failed_samples_cause_nonzero_exit_despite_runner_success(self):
        result = self.run_container("cohort", self.bam("good"), self.bam("bad"), FAIL_SAMPLE="bad")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("bad did not finish", result.stderr)
        self.assertTrue((self.root / "work/done/good.cleanup").exists())
        self.assertFalse((self.root / "reports.txt").exists())

    def test_pipeline_and_report_failures_propagate(self):
        bam = self.bam("S")
        result = self.run_container("sample", bam, PIPELINE_EXIT="7")
        self.assertEqual(result.returncode, 7, result.stderr)
        self.assertFalse((self.root / "reports.txt").exists())
        result = self.run_container("sample", bam, REPORT_FAIL="1")
        self.assertEqual(result.returncode, 9, result.stderr)


if __name__ == "__main__":
    unittest.main()
