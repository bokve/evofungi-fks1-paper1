"""Strict checks of the already-inspected MASTER AA/exact-codon schemas."""
import csv
import itertools
import math
from pathlib import Path

CODON_TABLE = dict(zip((''.join(c) for c in itertools.product('TCAG', repeat=3)), 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
FITNESS = tuple('fitness_' + c for c in ('none', 'anidulafungin', 'caspofungin', 'micafungin'))


def read_csv(path, required=()):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)) or not set(required).issubset(columns):
            raise ValueError(f'Duplicate or missing required columns in {path}')
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f'Inconsistent row width in {path}')
    return rows


def revcomp(sequence):
    if set(sequence) - set('ACGT'):
        raise ValueError('Ambiguous DNA bases')
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def validate_master(aa_rows, dna_rows, genome=None):
    by_aa = {}
    for r in aa_rows:
        identity = r['aa_variant_id']
        if identity in by_aa:
            raise ValueError(f'Duplicate AA phenotype identity: {identity}')
        if r['gene'] != 'FKS1' or identity != f"FKS1:{r['wt_aa']}{int(r['aa_position'])}{r['mut_aa']}":
            raise ValueError(f'AA identity fields disagree: {identity}')
        if int(r['aa_position']) < 1 or r['wt_aa'] == r['mut_aa']:
            raise ValueError(f'Not a valid missense AA identity: {identity}')
        by_aa[identity] = r
    seen = set()
    codons = {key: [] for key in by_aa}
    for r in dna_rows:
        identity = r['aa_variant_id']
        key = r['codon_variant_id']
        if key in seen:
            raise ValueError(f'Duplicate exact-codon identity: {key}')
        seen.add(key)
        if identity not in by_aa:
            raise ValueError(f'DNA identity absent from AA phenotype table: {identity}')
        a = by_aa[identity]
        if any(r[c] != a[c] for c in ('gene', 'hotspot', 'aa_position', 'wt_aa', 'mut_aa')):
            raise ValueError(f'DNA/AA identity mismatch: {identity}')
        wt, alt = r['coding_wt_codon'], r['coding_mut_codon']
        if CODON_TABLE.get(wt) != a['wt_aa'] or CODON_TABLE.get(alt) != a['mut_aa']:
            raise ValueError(f'Codon translation mismatch: {key}')
        if key != f'{identity}:{wt}>{alt}':
            raise ValueError(f'Codon identifier disagrees with fields: {key}')
        if r['phenotype_unit'] != 'amino_acid_substitution' or r['back_translated'].lower() != 'false':
            raise ValueError(f'Phenotype unit/back-translation mismatch: {key}')
        convert = revcomp if r['strand'] == '-' else (lambda x: x)
        if r['strand'] not in {'+', '-'} or convert(wt) != r['genomic_ref_codon'] or convert(alt) != r['genomic_alt_codon']:
            raise ValueError(f'Genomic/coding strand mismatch: {key}')
        start, end = int(r['genomic_codon_start_1based']), int(r['genomic_codon_end_1based'])
        if end-start != 2 or start < 1:
            raise ValueError(f'Invalid genomic codon coordinates: {key}')
        if genome is not None:
            if r['contig'] not in genome or genome[r['contig']][start-1:end] != r['genomic_ref_codon']:
                raise ValueError(f'Frozen genome REF mismatch: {key}')
        for c in FITNESS:
            x, y = float(r[c]), float(a[c])
            if not math.isfinite(x) or not math.isfinite(y) or not math.isclose(x,y,rel_tol=0,abs_tol=1e-12):
                raise ValueError(f'Copied AA phenotype mismatch: {key}/{c}')
        codons[identity].append(alt)
    for identity, values in codons.items():
        expected = by_aa[identity]['experimental_mut_codons'].split('|')
        # Verified on actual multi-codon MASTER rows (e.g. FKS1:D1353L).
        if sorted(values) != sorted(expected) or len(values) != int(by_aa[identity]['n_experimental_codons']):
            raise ValueError(f'AA exact-codon membership/count mismatch: {identity}')
    return {'AA_variants': len(by_aa), 'physical_positions': len({r['aa_position'] for r in aa_rows}), 'exact_codon_rows': len(seen), 'phenotype_unit': 'AA_substitution', 'genome_REF_checked': genome is not None}
