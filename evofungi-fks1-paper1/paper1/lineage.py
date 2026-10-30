"""Resolve Method 3 source-unit ambiguity without refitting its statistics."""
import numpy as np
import pandas as pd
from .stages import atomic_json,digest

def run(inputs,attempt,config):
    keys=['strain','paralog','hotspot','Mutated_seq','compound','Replicate','nt_seq']
    source=pd.read_csv(inputs['raw_nt_scores'],low_memory=False,usecols=keys+['confidence_score','mutation_type','Nham_aa','Nham_codons','s_T2_T0'])
    used=pd.read_csv(inputs['method3_long'])
    for column in ['confidence_score','Nham_aa','Nham_codons']:
        source[column]=pd.to_numeric(source[column],errors='coerce')
    source=source[(source.strain=='BY4741')&(source.paralog=='FKS1')&source.hotspot.isin(['HS1','HS2','HS3'])&
                  (source.confidence_score==1)&(source.mutation_type.str.lower()=='missense')&(source.Nham_aa==1)&(source.Nham_codons==1)].copy()
    for frame in [source,used]:
        frame['Replicate']=pd.to_numeric(frame.Replicate,errors='raise').astype(int)
        if frame.duplicated(keys).any():raise ValueError('Ambiguous exact nucleotide/condition/replicate source key')
    if len(used)!=5212 or len(source)!=len(used):raise ValueError('Historical Method3 source population differs')
    joined=used.merge(source[keys+['s_T2_T0']],on=keys,how='outer',validate='one_to_one',indicator=True,suffixes=('_used','_source'))
    if not (joined['_merge']=='both').all():raise ValueError('Historical Method3 row identities differ from NT source')
    np.testing.assert_allclose(joined.s_T2_T0_used,joined.s_T2_T0_source,rtol=0,atol=1e-12)
    dna=pd.read_csv(inputs['dna'])
    mapping=dna[['hotspot','experimental_hotspot_nt_seq','aa_variant_id','aa_position','codon_variant_id']].rename(columns={'experimental_hotspot_nt_seq':'nt_seq'})
    if mapping.duplicated(['hotspot','nt_seq']).any():raise ValueError('Ambiguous frozen genotype map')
    mapped=used.merge(mapping,on=['hotspot','nt_seq'],how='left',validate='many_to_one',suffixes=('_used','_master'))
    for col in ['aa_variant_id','aa_position','codon_variant_id']:
        if not (mapped[col+'_used']==mapped[col+'_master']).all():raise ValueError('Method3 canonical genotype/AA map mismatch: '+col)
    median=used.groupby(['aa_variant_id','compound','Replicate']).s_T2_T0.transform('median')
    np.testing.assert_allclose(median,used.experimental_AA_median_s_T2_T0,rtol=0,atol=1e-12)
    np.testing.assert_allclose(used.s_T2_T0-median,used.experimental_within_AA_genotype_residual,rtol=0,atol=1e-12)
    groups=used.groupby(['aa_variant_id','compound','Replicate']).agg(n_codons=('nt_seq','nunique'),n_distinct_source_scores=('s_T2_T0','nunique'))
    checks=joined[keys+['aa_variant_id','s_T2_T0_used','s_T2_T0_source']].copy();checks['source_matches']=True
    checks.to_csv(attempt/'METHOD3_SOURCE_ROW_CHECKS.csv',index=False)
    atomic_json(attempt/'METHOD3_UNIT_AUDIT.json',{'status':'EXECUTED_VALIDATED','scope':'source-unit/row validation, not Method3 refit or new correlation',
       'rows_checked':len(used),'exact_codons':int(used.codon_variant_id.nunique()),'AA':int(used.aa_variant_id.nunique()),
       'all_exported_values_match_source_NT_coefficients':True,'residual_medians_reproduced':True,
       'canonical_genotype_AA_map_checked':True,
       'multi_codon_AA_condition_replicate_groups':int((groups.n_codons>=2).sum()),
       'groups_with_distinct_NT_coefficients':int(((groups.n_codons>=2)&(groups.n_distinct_source_scores>=2)).sum()),
       'not_AA_aggregate_copied_to_codons':True,'not_independent_biological_ground_truth':True,
       'measurement_limit':'Distinct nucleotide-level source coefficients retain shared assay/normalization/replicate dependence; historical weak reliability remains.',
       'sources':{role:{'path':str(path),'sha256':digest(path)} for role,path in inputs.items() if role in ['raw_nt_scores','method3_long','dna']}})
    return ['METHOD3_SOURCE_ROW_CHECKS.csv','METHOD3_UNIT_AUDIT.json']
