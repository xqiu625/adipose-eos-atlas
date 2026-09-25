"""Shell setup tests with mocked modules; real HPCC batch validation still required."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class ModuleSetupTests(unittest.TestCase):
    def run_runner(self, stage="decontx", fail_purge=False):
        runner = Path(__file__).resolve().parents[1] / "scripts/run_preprocessing_task.sh"
        with tempfile.TemporaryDirectory() as base:
            env = dict(os.environ, SLURM_CPUS_PER_TASK="2", SLURM_ARRAY_TASK_ID="1",
                       RERUN_R_LIB=base, CMDSTAN=base, FAIL_PURGE=str(int(fail_purge)))
            driver = '''
module() {
  echo "MODULE $*"
  if [[ "$1" == purge && "$FAIL_PURGE" == 1 ]]; then return 1; fi
  return 0
}
Rscript() { echo "RSCRIPT $*"; }
export -f module Rscript
exec bash "$@"
'''
            return subprocess.run(["bash", "-c", driver, "module-test", str(runner),
                                   stage, base, "worker", "/bin/echo"],
                                  env=env, capture_output=True, text=True)

    def test_purge_precedes_pinned_modules_and_version_check(self):
        result = self.run_runner()
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(lines[:4], ["MODULE purge", "MODULE load gcc/12.2.0",
                                    "MODULE load R/4.5.0", "MODULE list"])
        self.assertTrue(lines[4].startswith("RSCRIPT --vanilla -e stopifnot"))
        self.assertIn("4.5.0", lines[4])

    def test_failed_purge_stops_worker(self):
        result = self.run_runner(fail_purge=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("MODULE load", result.stdout)
        self.assertNotIn("parallel_preprocessing.py", result.stdout)

    def test_python_stage_does_not_purge_modules(self):
        result = self.run_runner(stage="scrublet")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("MODULE", result.stdout)
        self.assertNotIn("RSCRIPT", result.stdout)
