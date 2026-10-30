import unittest
from pathlib import Path

from paper1.identity import read_csv, revcomp, validate_master


class IdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1] / 'reference_audit/evidence/02_processed_data'
        cls.aa = read_csv(root / 'MASTER_AA_TABLE_v0.1.csv')
        cls.dna = read_csv(root / 'MASTER_VARIANT_TABLE_v0.1_exact_codons.csv')

    def test_real_master_identity_strand_codon_and_phenotype_units(self):
        result = validate_master(self.aa, self.dna)
        self.assertEqual(result['AA_variants'], 435)
        self.assertEqual(result['exact_codon_rows'], 589)
        self.assertFalse(result['genome_REF_checked'])

    def test_wrong_strand_and_copied_fitness_fail(self):
        for field, value in [('strand', '+'), ('fitness_none', '123.456')]:
            rows = [dict(r) for r in self.dna]
            rows[0][field] = value
            with self.assertRaises(ValueError):
                validate_master(self.aa, rows)

    def test_duplicate_variant_cannot_be_independent_observation(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate exact-codon'):
            validate_master(self.aa, [*self.dna, self.dna[0]])

    def test_reverse_complement(self):
        self.assertEqual(revcomp('CTT'), 'AAG')


if __name__ == '__main__':
    unittest.main()
