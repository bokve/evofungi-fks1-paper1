"""Source-keyed aggregate panels; no phenotype is manufactured for codons."""
import json
import numpy as np
import pandas as pd

from .esm_contract import read_protein
from .stages import atomic_json, digest

DRUGS = ('anidulafungin', 'caspofungin', 'micafungin')
FITNESS = ['fitness_none', *['fitness_' + d for d in DRUGS]]
ENDPOINTS = ['fitness_none', *['drug_specific_' + d for d in DRUGS]]


def keyed(path):
    frame = pd.read_csv(path)
    if frame.aa_variant_id.duplicated().any():
        raise ValueError('Duplicate canonical AA identity: ' + str(path))
    return frame.set_index('aa_variant_id', drop=False).sort_index()


def replace_panel(raw, replacements):
    expected = {(aa, condition) for aa in ['FKS1:V641W', 'FKS1:L642K']
                for condition in ['none', *DRUGS]}
    keys = list(zip(replacements.aa_variant_id, replacements.compound))
    if len(keys) != 8 or set(keys) != expected:
        raise ValueError('Exactly two specified AA by four conditions required')
    result = raw.copy(deep=True)
    for row in replacements.itertuples():
        result.loc[row.aa_variant_id, 'fitness_' + row.compound] = float(row.s)
    for drug in DRUGS:
        result['drug_specific_' + drug] = result['fitness_' + drug] - result.fitness_none
    untouched = ~raw.aa_variant_id.isin(['FKS1:V641W', 'FKS1:L642K'])
    if not np.array_equal(raw.loc[untouched, FITNESS].to_numpy(), result.loc[untouched, FITNESS].to_numpy()):
        raise ValueError('Unexpected phenotype changes outside replacements')
    return result


def prepare(inputs, attempt, config):
    raw = keyed(inputs['aa'])
    phenotype, evo, gemme = [keyed(inputs[x]) for x in ['phenotype', 'evo2_aa', 'gemme']]
    if len(raw) != 435 or raw.aa_position.nunique() != 28:
        raise ValueError('Primary AA/position universe changed')
    for label, data in [('phenotype', phenotype), ('Evo2', evo), ('GEMME', gemme)]:
        if set(data.index) != set(raw.index):
            raise ValueError('Primary identity mask mismatch: ' + label)
        for col in ['aa_position', 'hotspot', 'wt_aa', 'mut_aa']:
            if not np.array_equal(raw[col], data.loc[raw.index, col]):
                raise ValueError('Canonical mapping mismatch: ' + label + '/' + col)
    for col in FITNESS:
        np.testing.assert_allclose(raw[col], phenotype.loc[raw.index, col], rtol=0, atol=1e-12)
    for drug in DRUGS:
        col = 'drug_specific_' + drug
        raw[col] = phenotype.loc[raw.index, col]
        np.testing.assert_allclose(raw[col], raw['fitness_' + drug] - raw.fitness_none, rtol=0, atol=1e-12)
    raw['evo2_primary_AA_score'] = evo.loc[raw.index, 'evo2_primary_AA_score']
    raw['gemme_primary'] = gemme.loc[raw.index, 'gemme_primary']
    np.testing.assert_allclose(raw.evo2_primary_AA_score, gemme.loc[raw.index, 'evo2_local53'], rtol=0, atol=1e-12)
    np.testing.assert_allclose(raw.gemme_primary, gemme.loc[raw.index, 'gemme_norm_combi'], rtol=0, atol=0)
    if set(gemme.gemme_sign_transform) != {'NONE_RAW_NORMALIZED_SCORE'}:
        raise ValueError('Source comparator direction changed')
    sequence = read_protein(inputs['protein'], config['protein_record_id'])
    for row in raw.itertuples():
        if sequence[row.aa_position - 1] != row.wt_aa or row.aa_variant_id != f'FKS1:{row.wt_aa}{row.aa_position}{row.mut_aa}':
            raise ValueError('Protein/AA identity mismatch')
    classes, source_sets = [], []
    for hs, group in raw.groupby('hotspot'):
        start, end = int(group.aa_position.min()), int(group.aa_position.max())
        reference = sequence[start - 1:end]
        source = inputs['classified_' + hs]
        data = pd.read_csv(source)
        for rownum, row in data.iterrows():
            seq = row.aa_seq
            if len(seq) != len(reference):
                raise ValueError('Classified source hotspot sequence length mismatch')
            changed = [i for i, (a, b) in enumerate(zip(reference, seq)) if a != b]
            if len(changed) > 1:
                raise ValueError('Classified source contains multiple substitutions')
            if not changed:
                identity, kind = hs + ':WT', 'WT'
            else:
                i = changed[0]
                identity = f'FKS1:{reference[i]}{start+i}{seq[i]}'
                kind = 'stop' if seq[i] == '*' else 'missense'
            classes.append(dict(row, aa_variant_id=identity, hotspot=hs, identity_type=kind,
                                source_path=str(source), source_sha256=digest(source), source_CSV_row=int(rownum)+2,
                                phenotype_source='AUTHOR_CLASSIFIED_REPLACEMENT_VERSION', primary_member=identity in raw.index))
        source_sets.append({'hotspot':hs,'unique_sequences_including_WT':data.aa_seq.nunique()})
    classes = pd.DataFrame(classes)
    if classes.duplicated(['aa_variant_id', 'compound']).any():
        raise ValueError('Ambiguous classified AA/condition join')
    repl = classes[classes.aa_variant_id.isin(['FKS1:V641W','FKS1:L642K'])].copy()
    validation = pd.read_csv(inputs['author_validation_DMS_missing_estimates'])
    validation['aa_variant_id'] = [f'FKS1:{r.wt_aa}{int(r.aa_pos)}{r.alt_aa}' for r in validation.itertuples()]
    if validation.duplicated(['aa_variant_id', 'compound']).any():
        raise ValueError('Ambiguous original author validation table')
    checked = validation.merge(classes[['aa_variant_id','compound','aa_seq','s']],
                               on=['aa_variant_id','compound'],validate='one_to_one',suffixes=('_validation','_classified'))
    if len(checked)!=len(validation) or len(checked)!=28 or not (checked.aa_seq_validation==checked.aa_seq_classified).all():
        raise ValueError('Author rescue sequence/condition lineage mismatch')
    np.testing.assert_allclose(checked.s_validation,checked.s_classified,rtol=0,atol=1e-12)
    checked['validation_source_sha256'] = digest(inputs['author_validation_DMS_missing_estimates'])
    checked['rescue_script_sha256'] = digest(inputs['author_rescue_missing_mutants'])
    checked.to_csv(attempt/'AUTHOR_VALIDATION_LINEAGE.csv',index=False)
    replaced = replace_panel(raw, repl)
    repl['original_raw_value'] = [raw.loc[r.aa_variant_id, 'fitness_'+r.compound] for r in repl.itertuples()]
    repl['replacement_value'] = repl.s
    repl['validation_source'] = str(inputs['author_validation_DMS_missing_estimates'])
    repl['validation_sha256'] = digest(inputs['author_validation_DMS_missing_estimates'])
    repl['join_key'] = 'canonical hotspot WT sequence with exactly one AA substitution + compound'
    for name, data in [('RAW435',raw),('AUTHOR_REPLACEMENT435',replaced)]:
        data = data.copy();data['panel'] = name
        for role in ['aa','phenotype','evo2_aa','gemme']:
            data[role+'_source_path'] = str(inputs[role]);data[role+'_sha256'] = digest(inputs[role])
        data.to_csv(attempt/(name+'.csv'),index=False)
        long = data.melt(id_vars=['aa_variant_id','aa_position','hotspot','panel'],value_vars=FITNESS,var_name='condition',value_name='fitness')
        long.to_csv(attempt/(name+'_LONG.csv'),index=False)
    repl.to_csv(attempt/'CORRECTION_LEDGER.csv',index=False)
    classes['raw_fitness'] = [raw.loc[r.aa_variant_id,'fitness_'+r.compound] if r.aa_variant_id in raw.index else np.nan for r in classes.itertuples()]
    classes['source_fitness_matches_raw'] = np.isclose(classes.s, classes.raw_fitness, atol=1e-12,rtol=0)
    classes.to_csv(attempt/'SOURCE_LABEL_MAPPING.csv',index=False)
    miss = set(classes.loc[classes.identity_type=='missense','aa_variant_id'])
    extras = sorted(miss-set(raw.index))
    if len(miss)!=440 or len(extras)!=5 or not set(raw.index).issubset(miss):
        raise ValueError('Source/project set definitions changed')
    classes.loc[classes.aa_variant_id.isin(extras)].to_csv(attempt/'VALIDATION_ONLY_EXTRAS.csv',index=False)
    note={'source_rows':source_sets,'missense':len(miss),'stop':classes.loc[classes.identity_type=='stop','aa_variant_id'].nunique(),
          'WT_hotspot_references':3,'project_subset':True,'extra_missense':extras,
          'source465_wording':'Operational table has 465 hotspot-AA identities including WT references; article wording is not silently reinterpreted as 465 missense.',
          'source_labels':'Annotations from author replacement source; not independent experimental validation or new clinical labels.',
          'two_AA_replacement_all_four_conditions':True,'historical_replicate_models_refit':False}
    atomic_json(attempt/'SOURCE465_NOTE.json',note)
    pd.DataFrame([{'aa_variant_id':k,'included':True,'reason':'validated historical complete-case mask','unit':'AA substitution'} for k in raw.index]).to_csv(attempt/'MASK_LEDGER.csv',index=False)
    return [p.name for p in attempt.iterdir() if p.is_file() and p.name!='started.json']
