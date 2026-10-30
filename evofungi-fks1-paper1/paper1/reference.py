"""Derive a frozen BY4741 CDS only when full-protein identity is demonstrated."""
import hashlib

from .identity import CODON_TABLE, revcomp


def derive_verified_reference(genome, dna_rows, query, expected_query_sha256):
    if hashlib.sha256(query.encode()).hexdigest() != expected_query_sha256:
        raise ValueError('Full query protein bytes disagree with authoritative source log')
    anchors = set()
    for row in dna_rows:
        strand = row['strand']
        p = int(row['aa_position'])
        if strand == '-':
            anchor = int(row['genomic_codon_end_1based']) - 1 + 3*(p-1)
        elif strand == '+':
            anchor = int(row['genomic_codon_start_1based']) - 1 - 3*(p-1)
        else:
            raise ValueError('Unknown source strand')
        anchors.add((row['contig'], strand, anchor))
    if len(anchors) != 1:
        raise ValueError('Master coordinates do not define one consistent coding origin')
    contig, strand, anchor = next(iter(anchors))
    sequence = genome[contig]
    length = 3*len(query)
    if strand == '-':
        start, end = anchor-length+1, anchor+1
        if start < 3 or end > len(sequence):
            raise ValueError('CDS/stop extends outside frozen contig')
        cds = revcomp(sequence[start:end])
        stop = revcomp(sequence[start-3:start])
    else:
        start, end = anchor, anchor+length
        if start < 0 or end+3 > len(sequence):
            raise ValueError('CDS/stop extends outside frozen contig')
        cds, stop = sequence[start:end], sequence[end:end+3]
    translated = ''.join(CODON_TABLE[cds[i:i+3]] for i in range(0, len(cds), 3))
    if not cds.startswith('ATG') or CODON_TABLE.get(stop) != '*':
        raise ValueError('Source-derived candidate lacks canonical start/adjacent in-frame stop')
    if translated != query:
        raise ValueError('Frozen BY4741 genomic translation differs from full authoritative query; do not substitute query protein')
    for row in dna_rows:
        p = int(row['aa_position'])
        if cds[3*(p-1):3*p] != row['coding_wt_codon'] or translated[p-1] != row['wt_aa']:
            raise ValueError('Source-derived full CDS does not match a frozen WT codon')
    return cds, translated, {'contig': contig, 'strand': strand, 'CDS_genomic_start_1based': start+1, 'CDS_genomic_end_1based': end, 'adjacent_stop_codon_coding_orientation': stop, 'protein_length': len(translated), 'full_protein_sha256': hashlib.sha256(translated.encode()).hexdigest(), 'full_query_identity': True, 'all_source_WT_codons_checked': len(dna_rows), 'derivation': 'all_MASTER_codon_coordinates_agree_on_origin; full_contiguous_genomic_translation_equals_authoritatively_verified_query; ATG_and_adjacent_stop_verified', 'not_annotation_claim': True}
