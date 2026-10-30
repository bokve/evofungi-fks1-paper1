import json
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from paper1.stages import run_stage, writer_lock


class StageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'input.txt'
        self.source.write_text('original')
        self.calls = 0

    def run_example(self, params=None, mutate=False):
        def action(attempt):
            self.calls += 1
            (attempt / 'output.txt').write_text(self.source.read_text())
            if mutate:
                self.source.write_text('changed during execution')
            return ['output.txt']
        def validate(attempt, outputs):
            if not (attempt / outputs[0]).read_text():
                raise ValueError('Empty output')
        return run_stage('example', self.root / 'results', {'source': self.source}, {'test': Path(__file__)}, params or {}, action, validate)

    def test_resume_checks_outputs_and_preserves_previous_attempt(self):
        first = self.run_example()
        self.assertEqual(self.run_example()['status'], 'REUSED_VALIDATED')
        self.assertEqual(self.calls, 1)
        (Path(first['directory']) / 'output.txt').write_text('corrupted')
        second = self.run_example()
        self.assertEqual(second['status'], 'EXECUTED_VALIDATED')
        self.assertNotEqual(first['directory'], second['directory'])
        self.assertEqual((Path(first['directory']) / 'output.txt').read_text(), 'corrupted')

    def test_changed_input_and_params_invalidate_resume(self):
        a = self.run_example()
        self.source.write_text('new input')
        b = self.run_example()
        c = self.run_example({'seed': 42})
        self.assertEqual(len({a['fingerprint'], b['fingerprint'], c['fingerprint']}), 3)

    def test_mutated_input_cannot_be_completed(self):
        with self.assertRaisesRegex(RuntimeError, 'changed during execution'):
            self.run_example(mutate=True)
        self.assertEqual(list((self.root / 'results').rglob('completion.json')), [])
        self.assertEqual(len(list((self.root / 'results').rglob('failure.json'))), 1)

    def test_second_writer_is_rejected(self):
        with writer_lock(self.root / 'results'):
            with self.assertRaisesRegex(RuntimeError, 'Writer lock exists'):
                with writer_lock(self.root / 'results'):
                    self.fail('Concurrent writer admitted')
        self.assertFalse((self.root / 'results/.paper1-writer.lock').exists())

    def test_relocation_reuses_content_but_software_change_invalidates(self):
        with patch('paper1.stages.software',return_value={'python':'fixed','packages':{'fixture':'1'}}):
            a=self.run_example()
            moved=self.root/'relocated.txt';moved.write_bytes(self.source.read_bytes());self.source=moved
            b=self.run_example()
        self.assertEqual(b['status'],'REUSED_VALIDATED')
        self.assertEqual(a['fingerprint'],b['fingerprint'])
        with patch('paper1.stages.software',return_value={'python':'fixed','packages':{'fixture':'2'}}):
            c=self.run_example()
        self.assertNotEqual(c['fingerprint'],a['fingerprint'])


if __name__ == '__main__':
    unittest.main()
