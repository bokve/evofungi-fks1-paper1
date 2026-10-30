"""Small, standard-library checks of supplied evidence; never run a notebook."""
import csv
import hashlib
import json
import statistics
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def ranks(values):
    ordered = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        for i in ordered[start:end]:
            result[i] = (start + 1 + end) / 2
        start = end
    return result


def main():
    root = Path(__file__).resolve().parents[1] / 'reference_audit'
    manifest = list(csv.DictReader((root / 'PACKAGE_MANIFEST.csv').open()))
    failures = []
    for row in manifest:
        p = root / row['package_path']
        if not p.is_file() or p.stat().st_size != int(row['bytes']) or sha(p) != row['sha256']:
            failures.append(row['package_path'])
    sources = {}

    def load(name):
        paths = list((root / 'evidence').rglob(name))
        if len(paths) != 1:
            raise ValueError(f'Ambiguous/missing source: {name}: {paths}')
        p = paths[0]
        with p.open(newline='') as f:
            rows = list(csv.DictReader(f))
        sources[name] = {'path': str(p.relative_to(root)), 'sha256': sha(p), 'rows': len(rows)}
        return rows

    aa = load('MASTER_AA_TABLE_v0.1.csv')
    dna = load('MASTER_VARIANT_TABLE_v0.1_exact_codons.csv')
    theoretical = load('P1_NEW3D3_THEORETICAL814_EXACT_GENOTYPE_FEATURE_MANIFEST_v1.0.csv')
    alt = load('ALT490_CURRENT_HOSTED7B_v0.1.csv')
    ref = load('REF28_CURRENT_HOSTED7B_v0.1.csv')
    repeated = load('P1_NEW2B0A_MODEL_READY_TABLE_3884_v1.0.csv')
    y = load('FKS1_PHENOTYPE_COMPONENTS_v1.0.csv')
    xrows = load('EVO2_AA_SCORES_v1.0.csv')
    x = {r['aa_variant_id']: float(r['evo2_primary_AA_score']) for r in xrows}
    if len(x) != len(xrows) or len({r['aa_variant_id'] for r in y}) != len(y):
        raise ValueError('Duplicate AA join key')
    if set(x) != {r['aa_variant_id'] for r in y}:
        raise ValueError('AA join set mismatch')
    correlations = []
    for endpoint in ['fitness_none', 'drug_specific_anidulafungin', 'drug_specific_caspofungin', 'drug_specific_micafungin']:
        def rho(rows):
            return statistics.correlation(ranks([x[r['aa_variant_id']] for r in rows]), ranks([float(r[endpoint]) for r in rows]))
        positions = sorted({r['aa_position'] for r in y})
        correlations.append({'endpoint': endpoint, 'pooled': rho(y), 'median_position': statistics.median(rho([r for r in y if r['aa_position'] == p]) for p in positions), 'source_rows': 'all 435 AA rows; grouped by aa_position for median'})
    summary = load('P1_NEW2_MODULE_RESULT_SUMMARY_v1.0.csv')
    completion_path = root / 'evidence/04_evo2_scores/COMPLETION_v0.1.json'
    completion = json.loads(completion_path.read_text())
    report = {
        'purpose': 'handoff_validation_only_no_new_scientific_analysis',
        'package': {'entries': len(manifest), 'matched': len(manifest)-len(failures), 'failures': failures},
        'sources': sources,
        'universes': {'AA_rows': len(aa), 'AA_unique': len({r['aa_variant_id'] for r in aa}), 'positions': len({r['aa_position'] for r in aa}), 'observed_DNA_rows': len(dna), 'theoretical_DNA_rows': len(theoretical), 'theoretical_AA_groups': len({r['aa_group_key'] for r in theoretical}), 'scored_ALT_rows': len(alt), 'scored_ALT_AA_groups': len({r['aa_group_key'] for r in alt}), 'scored_REF_rows': len(ref), 'NEW2_measurement_rows': len(repeated), 'NEW2_AA_condition_groups': len({(r['aa_variant_id'], r['compound']) for r in repeated})},
        'completion_hash_checks': {'ALT': sources['ALT490_CURRENT_HOSTED7B_v0.1.csv']['sha256'] == completion['alt_sha256'], 'REF': sources['REF28_CURRENT_HOSTED7B_v0.1.csv']['sha256'] == completion['ref_sha256']},
        'historical_correlations_recomputed': correlations,
        'NEW2_recorded_LRT_not_refitted': [dict(row_number=i+2, **r) for i,r in enumerate(summary) if r['result_type'] == 'primary_omnibus_LRT'],
        'unresolved': ['Drive source identities and pinned-source aggregation/correction lineage', 'source465 vs project435 reconciliation', 'HS2 continuous-frequency source comparability', 'hosted immutable checkpoint SHA'],
        'not_executed': ['H2 refit', 'H3 resampling', 'new inference', 'phenotype correction', 'NEW4'],
    }
    print(json.dumps(report, indent=2))
    if failures or not all(report['completion_hash_checks'].values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
