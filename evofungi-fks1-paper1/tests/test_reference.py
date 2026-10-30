import hashlib
import unittest

from paper1.reference import derive_verified_reference


class ReferenceTests(unittest.TestCase):
    def test_minus_strand_full_genomic_translation_and_stop(self):
        rows = [{'strand': '-', 'contig': 'ctg', 'aa_position': '1', 'genomic_codon_start_1based': '7', 'genomic_codon_end_1based': '9', 'coding_wt_codon': 'ATG', 'wt_aa': 'M'}]
        h = hashlib.sha256(b'MK').hexdigest()
        cds, protein, record = derive_verified_reference({'ctg': 'TTATTTCAT'}, rows, 'MK', h)
        self.assertEqual((cds, protein), ('ATGAAA', 'MK'))
        self.assertTrue(record['full_query_identity'])
        with self.assertRaisesRegex(ValueError, 'adjacent in-frame stop'):
            derive_verified_reference({'ctg': 'CCCTTTCAT'}, rows, 'MK', h)
        with self.assertRaisesRegex(ValueError, 'differs from full authoritative query'):
            derive_verified_reference({'ctg': 'TTAATTCAT'}, rows, 'MK', h)


if __name__ == '__main__':
    unittest.main()
