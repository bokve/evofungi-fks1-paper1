import unittest

from paper1.esm_contract import targets, window


class ESMContractTests(unittest.TestCase):
    def test_long_protein_window_preserves_target_and_length_at_edges(self):
        sequence = 'ACDEFGHIKLMNPQRSTVWY' * 100
        for p in [1, 2, 511, 512, 639, 1353, 1999, 2000]:
            w = window(sequence, p)
            self.assertEqual(len(w['sequence']), 1022)
            self.assertEqual(w['sequence'][w['local_index_0based']], sequence[p - 1])
            self.assertEqual(w['start_0based'] + w['local_index_0based'], p - 1)
        self.assertEqual(window(sequence, 1)['start_0based'], 0)
        self.assertEqual(window(sequence, 2000)['end_0based_exclusive'], 2000)

    def test_short_sequence_not_padded_or_truncated(self):
        self.assertEqual(window('ACD', 2)['sequence'], 'ACD')
        with self.assertRaises(ValueError):
            window('ACD', 0)

    def test_wrong_ref_and_duplicate_variant_stop_preparation(self):
        row = {'aa_variant_id': 'FKS1:A1C', 'aa_position': '1', 'wt_aa': 'A', 'mut_aa': 'C'}
        with self.assertRaisesRegex(ValueError, 'WT mismatch'):
            targets('CC', [row])
        with self.assertRaisesRegex(ValueError, 'Duplicate AA variant'):
            targets('AC', [row, row])

    def test_targets_are_phenotype_blind_and_use_all_five_checkpoints(self):
        row = {'aa_variant_id': 'FKS1:A1C', 'aa_position': '1', 'wt_aa': 'A', 'mut_aa': 'C', 'fitness_none': 'unreadable sentinel'}
        a = targets('AC', [row])
        b = targets('AC', [dict(row, fitness_none='different')])
        self.assertEqual(a, b)
        self.assertEqual(len({r['checkpoint'] for r in a}), 5)


if __name__ == '__main__':
    unittest.main()
