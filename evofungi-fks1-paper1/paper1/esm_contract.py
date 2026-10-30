"""Phenotype-blind ESM-1v preparation. No model downloads or inference here."""
import hashlib
from pathlib import Path

CANONICAL_AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
CHECKPOINTS = tuple(f'esm1v_t33_650M_UR90S_{i}' for i in range(1, 6))


def read_protein(path, record_id):
    records, current, parts = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if current is not None:
                if current in records:
                    raise ValueError('Duplicate FASTA record')
                records[current] = ''.join(parts)
            current, parts = line[1:].split()[0], []
        elif line.strip():
            if current is None:
                raise ValueError('FASTA sequence without header')
            parts.append(line.strip())
    if current is not None:
        if current in records:
            raise ValueError('Duplicate FASTA record')
        records[current] = ''.join(parts)
    if record_id not in records:
        raise ValueError('Explicit frozen FASTA record ID missing; no automatic substitute')
    sequence = records[record_id]
    if not sequence or set(sequence) - CANONICAL_AA:
        raise ValueError('Protein sequence must contain canonical amino acids only')
    return sequence


def window(sequence, position_1based):
    p = position_1based - 1
    if not 0 <= p < len(sequence):
        raise ValueError('Mutation position outside protein')
    start = 0 if len(sequence) <= 1022 else max(0, min(p - 511, len(sequence) - 1022))
    subsequence = sequence[start:start + 1022]
    return {'start_0based': start, 'end_0based_exclusive': start + len(subsequence), 'local_index_0based': p - start, 'sequence': subsequence}


def targets(sequence, aa_rows):
    """Use sequence/variant identities only; never inspect fitness columns."""
    seen, positions = set(), {}
    for row in aa_rows:
        identity = row['aa_variant_id']
        if identity in seen:
            raise ValueError(f'Duplicate AA variant: {identity}')
        seen.add(identity)
        position = int(row['aa_position'])
        wt, alt = row['wt_aa'], row['mut_aa']
        if wt not in CANONICAL_AA or alt not in CANONICAL_AA or wt == alt:
            raise ValueError(f'Invalid missense identity: {identity}')
        if identity != f'FKS1:{wt}{position}{alt}':
            raise ValueError(f'AA identifier disagrees with fields: {identity}')
        w = window(sequence, position)
        if sequence[position - 1] != wt:
            raise ValueError(f'Frozen protein WT mismatch: {identity}')
        if position in positions and positions[position]['wt_aa'] != wt:
            raise ValueError('Inconsistent WT at one physical position')
        positions[position] = {'aa_position': position, 'wt_aa': wt, **w}
    return [{'checkpoint': checkpoint, **positions[p], 'window_sha256': hashlib.sha256(positions[p]['sequence'].encode()).hexdigest()} for checkpoint in CHECKPOINTS for p in sorted(positions)]
