"""Run preparation or bounded Paper 1 execution from a versioned config.

Historical preparation remains separate. Execution permits only the fixed
remote ESM comparator, with hash gates and no historical model refits.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from .esm_contract import read_protein, targets
from .identity import FITNESS, read_csv, validate_master
from .reference import derive_verified_reference
from .stages import atomic_json, digest, run_stage


def fasta(path):
    records, current, parts = {}, None, []
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if current is not None:
                    if current in records:
                        raise ValueError('Duplicate genome FASTA record')
                    records[current] = ''.join(parts).upper()
                current, parts = line[1:].split()[0], []
            elif line.strip():
                if current is None:
                    raise ValueError('Genome FASTA sequence has no header')
                parts.append(line.strip())
    if current is not None:
        if current in records:
            raise ValueError('Duplicate genome FASTA record')
        records[current] = ''.join(parts).upper()
    return records


def indexed(rows, label):
    result = {r['aa_variant_id']: r for r in rows}
    if len(rows) != len(result):
        raise ValueError(f'Duplicate AA key in {label}')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--stage', required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    if config.get('execution_run'):
        from .execution import execute
        execute(args.config, args.stage)
        return
    root = Path(config['project_root']).resolve(strict=True)
    output = (root / config['output_relative']).resolve()
    if not output.is_relative_to(root) or output == root:
        raise ValueError('New output directory must lie inside the verified project root')
    required = {
        'identity': ['aa', 'dna', 'phenotype', 'evo2_aa', 'genome', 'reference_provenance', 'new2_rows', 'alt490', 'ref28', 'completion'],
        'prepare-esm': ['aa', 'protein', 'reference_derivation'],
        'derive-reference': ['dna', 'genome', 'reference_provenance', 'query_msa', 'query_qc', 'query_ref_log'],
        'integrate-records': config.get('integration_roles', []),
    }[args.stage]
    if not required:
        raise ValueError('Explicit integration roles required')
    inputs = {'configuration': args.config.resolve()}
    for role in required:
        if role not in config['inputs']:
            raise FileNotFoundError(f'BLOCKED_INPUT: explicit input role required: {role}')
        spec = config['inputs'][role]
        p = (root / spec['relative_path']).resolve(strict=True)
        if not p.is_relative_to(root) or not p.is_file() or digest(p) != spec['sha256']:
            raise ValueError(f'Hash/path identity failed for {role}')
        inputs[role] = p

    def identity_action(attempt):
        aa = read_csv(inputs['aa'], ('aa_variant_id', *FITNESS))
        dna = read_csv(inputs['dna'], ('codon_variant_id', 'phenotype_unit'))
        reference = json.loads(inputs['reference_provenance'].read_text())
        if digest(inputs['genome']) != reference['genome_fasta_sha256']:
            raise ValueError('Frozen genome does not match original reference provenance')
        checked = validate_master(aa, dna, fasta(inputs['genome']))
        for key, expected in {'AA_variants': 435, 'physical_positions': 28, 'exact_codon_rows': 589}.items():
            if checked[key] != expected:
                raise ValueError(f'Frozen universe mismatch: {key}')
        a = indexed(aa, 'MASTER')
        phen = indexed(read_csv(inputs['phenotype'], ('aa_variant_id', *FITNESS)), 'phenotype')
        evo = indexed(read_csv(inputs['evo2_aa'], ('aa_variant_id', 'evo2_primary_AA_score')), 'Evo2 AA')
        if set(a) != set(phen) or set(a) != set(evo):
            raise ValueError('AA input masks differ; no implicit row loss allowed')
        for key, row in a.items():
            if not math.isfinite(float(evo[key]['evo2_primary_AA_score'])):
                raise ValueError(f'Nonfinite Evo2 AA score: {key}')
            for col in FITNESS:
                if not math.isclose(float(row[col]), float(phen[key][col]), rel_tol=0, abs_tol=1e-12):
                    raise ValueError(f'MASTER/phenotype disagreement: {key}/{col}')
            for drug in ('anidulafungin', 'caspofungin', 'micafungin'):
                difference = float(row['fitness_' + drug])-float(row['fitness_none'])
                if not math.isclose(difference, float(phen[key]['drug_specific_' + drug]), rel_tol=0, abs_tol=1e-12):
                    raise ValueError(f'Historical D-C identity failed: {key}/{drug}')
        repeated = read_csv(inputs['new2_rows'], ('aa_variant_id', 'compound', 'Replicate'))
        if len(repeated) != 3884 or len({(r['aa_variant_id'],r['compound'],r['Replicate']) for r in repeated}) != len(repeated):
            raise ValueError('NEW-2 measurement identity/count mismatch')
        if {r['aa_variant_id'] for r in repeated} != set(a) or len({(r['aa_variant_id'],r['compound']) for r in repeated}) != 1740:
            raise ValueError('NEW-2 AA/condition coverage mismatch')
        completion = json.loads(inputs['completion'].read_text())
        if digest(inputs['alt490']) != completion['alt_sha256'] or digest(inputs['ref28']) != completion['ref_sha256']:
            raise ValueError('Hosted ALT/REF bytes differ from original completion record')
        checked.update(status_scope='input_identity_only_not_P1_provenance_closure', matched_AA_mask=435, NEW2_measurements=3884, NEW2_AA_condition_groups=1740, hosted_checkpoint_identity='UNVERIFIED', source_aggregation_and_corrections='REQUIRE_SEPARATE_ORIGINAL_SOURCE_AUDIT')
        atomic_json(attempt / 'identity_preflight.json', checked)
        return ['identity_preflight.json']

    def esm_action(attempt):
        if not config.get('protein_record_id') or not config.get('protein_derivation_evidence'):
            raise ValueError('Explicit frozen protein record and derivation evidence required')
        evidence = json.loads(inputs['reference_derivation'].read_text())
        if evidence.get('full_query_identity') is not True:
            raise ValueError('Full frozen-genome protein identity has not been demonstrated')
        sequence = read_protein(inputs['protein'], config['protein_record_id'])
        jobs = targets(sequence, read_csv(inputs['aa']))
        if len(jobs) != 140:
            raise ValueError('Expected five checkpoints by 28 source-verified positions')
        atomic_json(attempt / 'esm_target_jobs.json', jobs)
        atomic_json(attempt / 'esm_contract.json', {'mode': 'masked_marginal_ALT_minus_WT', 'ensemble': 'arithmetic_mean_of_all_five_checkpoints', 'window_max_AA': 1022, 'association_not_executed': True, 'GPU_smoke_and_token_limit_verification': 'PENDING_REQUIRED', 'protein_derivation_evidence': config['protein_derivation_evidence']})
        return ['esm_target_jobs.json', 'esm_contract.json']

    def reference_action(attempt):
        reference = json.loads(inputs['reference_provenance'].read_text())
        if digest(inputs['genome']) != reference['genome_fasta_sha256']:
            raise ValueError('Genome/reference provenance mismatch')
        qc = read_csv(inputs['query_qc'])
        if len(qc) != 1:
            raise ValueError('Expected one explicit query QC record')
        records = fasta(inputs['query_msa'])
        query_id = qc[0]['query_id']
        if query_id not in records or '-' in records[query_id]:
            raise ValueError('Explicit source query absent/gapped in query-anchored MSA')
        verified = json.loads(inputs['query_ref_log'].read_text())
        if qc[0]['query_sha256'] != verified['candidate_sha256']:
            raise ValueError('Historical MSA query identity disagrees with authoritative verification log')
        cds, protein, evidence = derive_verified_reference(fasta(inputs['genome']), read_csv(inputs['dna']), records[query_id], verified['candidate_sha256'])
        (attempt / 'BY4741_FKS1_verified.fasta').write_text('>BY4741_FKS1_derived_verified\n' + protein + '\n')
        (attempt / 'BY4741_FKS1_verified_CDS.fna').write_text('>BY4741_FKS1_derived_verified_CDS\n' + cds + '\n')
        atomic_json(attempt / 'reference_derivation.json', evidence)
        return ['BY4741_FKS1_verified.fasta', 'BY4741_FKS1_verified_CDS.fna', 'reference_derivation.json']

    def integration_action(attempt):
        ledger = []
        for role in required:
            source = inputs[role]
            source_hash = digest(source)
            if source.suffix.lower() == '.csv':
                for line, row in enumerate(read_csv(source), 2):
                    ledger.append({'role': role, 'source_relative_path': str(source.relative_to(root)), 'sha256': source_hash, 'CSV_record_1based_including_header': line, 'record': row, 'status': 'REUSED_SOURCE_RECORD_NOT_REFITTED'})
            elif source.suffix.lower() == '.json':
                ledger.append({'role': role, 'source_relative_path': str(source.relative_to(root)), 'sha256': source_hash, 'JSON_pointer': '', 'record': json.loads(source.read_text()), 'status': 'REUSED_SOURCE_RECORD_NOT_REFITTED'})
            else:
                raise ValueError('Integration requires explicit CSV or JSON evidence')
        atomic_json(attempt / 'claim_source_ledger.json', ledger)
        return ['claim_source_ledger.json']

    def validate(attempt, outputs):
        for relative in outputs:
            value = json.loads((attempt / relative).read_text()) if relative.endswith('.json') else (attempt / relative).read_text()
            if not value:
                raise ValueError('Empty preparation output')
    actions = {'identity': identity_action, 'derive-reference': reference_action, 'prepare-esm': esm_action, 'integrate-records': integration_action}
    code = {p.name: p for p in Path(__file__).parent.glob('*.py')}
    print(json.dumps(run_stage(args.stage, output, inputs, code, {'stage': args.stage, 'seed': 0}, actions[args.stage], validate), indent=2))


if __name__ == '__main__':
    main()
