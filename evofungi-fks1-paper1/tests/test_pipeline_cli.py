import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from paper1.stages import digest


class PipelineCLITests(unittest.TestCase):
    def test_real_source_ledger_is_hashed_and_resumable(self):
        repo = Path(__file__).resolve().parents[1]
        source = repo / 'reference_audit/evidence/05_results/paper1_hierarchical/P1_NEW2_MODULE_RESULT_SUMMARY_v1.0.csv'
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            local = root / 'NEW2.csv'
            local.write_bytes(source.read_bytes())
            cfg = root / 'config.json'
            config = {'project_root': str(root), 'output_relative': 'new-results', 'inputs': {'new2': {'relative_path': local.name, 'sha256': digest(local)}}, 'integration_roles': ['new2']}
            cfg.write_text(json.dumps(config))
            command = [sys.executable, '-m', 'paper1', '--config', str(cfg), '--stage', 'integrate-records']
            a = subprocess.run(command, cwd=repo, capture_output=True, text=True, check=True)
            record = json.loads(a.stdout)
            ledger = json.loads((Path(record['directory']) / 'claim_source_ledger.json').read_text())
            self.assertEqual(len(ledger), 25)
            self.assertEqual(ledger[0]['CSV_record_1based_including_header'], 2)
            self.assertEqual(ledger[0]['sha256'], digest(source))
            b = subprocess.run(command, cwd=repo, capture_output=True, text=True, check=True)
            self.assertEqual(json.loads(b.stdout)['status'], 'REUSED_VALIDATED')
            local.write_text('changed source')
            c = subprocess.run(command, cwd=repo, capture_output=True, text=True)
            self.assertNotEqual(c.returncode, 0)
            self.assertIn('Hash/path identity failed', c.stderr)


if __name__ == '__main__':
    unittest.main()
