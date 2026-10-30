import csv
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path

from tools.audit_handoff import ranks
from tools.remote_intake import inspect_file
from tools.run_remote_readiness import build_wrapper, confirmed, redact


class ReadinessTests(unittest.TestCase):
    def test_proxy_and_oauth_query_secrets_are_redacted_without_losing_error(self):
        original='SSLEOFError https://example.test/file?authuser=0&colab-runtime-proxy-token=FAKE.PROXY.VALUE&content=1 Authorization: Bearer FAKE_BEARER https://example.test/?state=FAKE_STATE&code=FAKE_CODE'
        result=redact(original)
        for secret in ['FAKE.PROXY.VALUE','FAKE_BEARER','FAKE_STATE','FAKE_CODE']:
            self.assertNotIn(secret,result)
        self.assertIn('SSLEOFError',result)
        self.assertIn('content=1',result)

    def test_remote_wrapper_uses_fresh_main_and_rejects_child_error(self):
        payload = b'import __main__; print(hasattr(__main__, "_WRAPPER_SENTINEL"))\n'
        code = '_WRAPPER_SENTINEL=True\n' + build_wrapper(payload, 'probe.py', 'marker ', 'hash')
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertTrue(result.stdout.startswith('False\n'))
        self.assertTrue(confirmed(result.stdout, 'marker ', 'hash'))
        bad = build_wrapper(b'raise RuntimeError("expected child failure")\n', 'fail.py', 'marker ', 'hash')
        result = subprocess.run([sys.executable, '-c', bad], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('expected child failure', result.stderr)
        self.assertFalse(confirmed(result.stdout, 'marker ', 'hash'))

    def test_cli_zero_exit_without_remote_success_is_failure(self):
        self.assertFalse(confirmed('Traceback: remote kernel exception\n', 'marker ', 'hash'))
        self.assertFalse(confirmed('marker {"success": true, "code_sha256": "old-code"}\n', 'marker ', 'hash'))
        self.assertTrue(confirmed('marker {"success": true, "code_sha256": "hash"}\n', 'marker ', 'hash'))

    def test_tied_ranks(self):
        self.assertEqual(ranks([20, 10, 20, 40]), [2.5, 1, 2.5, 4])

    def test_real_aa_schema_and_unit(self):
        root = Path(__file__).resolve().parents[1] / 'reference_audit/evidence'
        path = root / '02_processed_data/MASTER_AA_TABLE_v0.1.csv'
        record = inspect_file(path, root)
        self.assertEqual(record['status'], 'HASHED')
        self.assertFalse(record['duplicate_columns'])
        self.assertEqual(record['malformed_row_examples'], [])
        with path.open() as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(record['rows'], len({r['aa_variant_id'] for r in rows}))
        self.assertIn('fitness_none', record['columns'])

    def test_duplicate_columns_and_bad_width_are_visible(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path = root / 'bad.csv'
            path.write_text('aa,aa\nx,y,z\n')
            record = inspect_file(path, root)
            self.assertTrue(record['duplicate_columns'])
            self.assertEqual(record['malformed_row_examples'], [2])

    def test_invalid_json_is_error_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            path = root / 'bad.json'
            path.write_text('{')
            self.assertEqual(inspect_file(path, root)['status'], 'SCHEMA_ERROR')


if __name__ == '__main__':
    unittest.main()
