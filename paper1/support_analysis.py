"""Independent, exploratory H208 score-versus-categorical-support branch."""
import json
import numpy as np
import pandas as pd
from .esm_contract import read_protein
from .metrics import position_records,paired_positions,draw_positions,rank_summaries
from .stages import atomic_json,digest
from datetime import datetime,timezone


def run(inputs,attempt,config):
    membership=pd.read_csv(inputs['h_membership'])
    membership=membership.loc[membership.in_theoretical490==1].copy()
    alt=pd.read_csv(inputs['alt490'])
    if len(membership)!=490 or set(membership.exact_genotype_key)!=set(alt.exact_genotype_key):
        raise ValueError('Theoretical H score membership mismatch')
    if alt.exact_genotype_key.duplicated().any():raise ValueError('Duplicate exact H score')
    for col,old in [('aa_group_key','aa_group_key'),('canonical_aa_position','aa_position'),('wt_codon','wt_codon'),('mutant_codon','mutant_codon')]:
        a=membership.set_index('exact_genotype_key')[col];b=alt.set_index('exact_genotype_key').loc[a.index,old]
        if not np.array_equal(a,b):raise ValueError('H score identity disagreement: '+col)
    mult=pd.read_csv(inputs['multiplicity']).set_index('aa_group_key')
    ensemble=pd.read_csv(inputs['esm-score/ESM_ENSEMBLE.csv']).set_index('aa_variant_id')
    reference=read_protein(inputs['protein'],config['protein_record_id'])
    matrix=pd.read_csv(inputs['gemme_matrix'],sep=r'\s+',index_col=0)
    # Same source-label normalization as frozen M4-R1's corrected ingestion.
    matrix.index=[str(x).strip().replace('"','').upper() for x in matrix.index]
    matrix.columns=[str(x).strip().replace('"','') for x in matrix.columns]
    if matrix.index.duplicated().any() or matrix.columns.duplicated().any():raise ValueError('Ambiguous normalized GEMME matrix labels')
    matrix=matrix.apply(pd.to_numeric,errors='raise')
    validated=pd.read_csv(inputs['gemme'])
    # Reproduce every existing matrix-derived score before using previously
    # unmeasured AA groups. This is source validation, not new phenotype fit.
    for row in validated.itertuples():
        if reference[row.aa_position-1]!=row.wt_aa:raise ValueError('GEMME WT identity mismatch')
        value=float(matrix.loc[row.mut_aa,'V'+str(row.aa_position)])
        if not np.isclose(value,row.gemme_primary,rtol=0,atol=1e-12):raise ValueError('Existing GEMME matrix lookup not reproduced')
    full=membership.merge(alt[['exact_genotype_key','evo2_primary_delta_local53']],on='exact_genotype_key',validate='one_to_one')
    rows=[]
    for key,g in full.groupby('aa_group_key',sort=True):
        p=int(g.canonical_aa_position.iloc[0]);wt=g.wt_aa.iloc[0];mut=g.mutant_aa.iloc[0]
        identity=f'FKS1:{wt}{p}{mut}'
        if reference[p-1]!=wt:raise ValueError('Support AA/reference identity mismatch')
        m=mult.loc[key]
        if len(g)!=m.n_theoretical_codons:raise ValueError('Designed H codon denominator mismatch')
        value=float(matrix.loc[mut,'V'+str(p)])
        rows.append({'aa_variant_id':identity,'aa_group_key':key,'aa_position':p,'hotspot':g.hotspot.iloc[0],
                     'designed_codons':len(g),'PRE_HC_codons':int(m.n_PRE_HC_codons),'T0_HC_codons':int(m.n_T0_HC_codons),
                     'support_fraction':float(m.n_T0_HC_codons)/len(g),
                     'Evo2_theoretical_H_median':float(g.evo2_primary_delta_local53.median()),
                     'GEMME':value,'ESM':float(ensemble.loc[identity,'esm1v_ensemble']) if identity in ensemble.index else np.nan})
    panel=pd.DataFrame(rows)
    if len(panel)!=208:raise ValueError('Expected theoretical multi-codon 208-AA universe')
    feature_cols=['Evo2_theoretical_H_median','GEMME','ESM']
    missing=~np.isfinite(panel[feature_cols].to_numpy()).all(axis=1)
    gate={'candidate_AA':208,'exact_codons':490,'zero_T0_support_AA':int((panel.T0_HC_codons==0).sum()),
          'existing_GEMME_lookup_verified_on_RAW435':True,'missing_AA':panel.loc[missing,'aa_variant_id'].tolist(),
          'continuous_readcount_gate_used':False,'primary_Evo2_score_replaced':False,'fitness_added':False,
          'bridge_caveat':'H490 hosted bridge and 339 historical overlap support scoring consistency; immutable hosted checkpoint bytes remain UNVERIFIED.',
          'exploratory_protocol_sha256':config['protocol_sha256']}
    panel.to_csv(attempt/'EXPLORATORY_SUPPORT_AA_PANEL.csv',index=False)
    if missing.any():
        gate.update(status='BLOCKED_FEATURE_UNIVERSE',reason='Incomplete three-score coverage including unsupported AA; no survivor-only rescue or new 814 scoring')
        atomic_json(attempt/'HS2_SUPPORT_SCORE_GATE.json',gate)
        return ['EXPLORATORY_SUPPORT_AA_PANEL.csv','HS2_SUPPORT_SCORE_GATE.json']
    gate.update(status='EXECUTED_VALIDATED',interpretation='Exploratory categorical representation association, not lethality/causality/binomial biological replication')
    atomic_json(attempt/'SUPPORT_PANEL_LOCK.json',{'locked_utc':datetime.now(timezone.utc).isoformat(),
      'panel_sha256':digest(attempt/'EXPLORATORY_SUPPORT_AA_PANEL.csv'),'protocol_sha256':config['protocol_sha256'],
      'feature_gate_passed_before_association':True,'AA_groups':208,'exact_codons':490,
      'sources':{role:{'path':str(path),'sha256':digest(path)} for role,path in inputs.items() if role in ['h_membership','alt490','multiplicity','gemme_matrix','gemme','protein'] or role=='esm-score/ESM_ENSEMBLE.csv'},
      'Evo2_summary':'new median over existing theoretical H codons; historical primary AA score unchanged'})
    predictors={name:name for name in feature_cols}
    rec=position_records(panel,'support_fraction',predictors)
    positions=paired_positions(rec,list(predictors));hs=panel.groupby('aa_position').hotspot.first().to_dict()
    draws=draw_positions(positions,hs)
    summaries,reps=rank_summaries(rec,list(predictors),[('ESM','GEMME'),('ESM','Evo2_theoretical_H_median'),('Evo2_theoretical_H_median','GEMME')],draws)
    pd.DataFrame(summaries).to_csv(attempt/'SUPPORT_SCORE_ASSOCIATIONS.csv',index=False)
    rec.to_csv(attempt/'SUPPORT_POSITION_RESULTS.csv',index=False)
    panel.groupby(['hotspot','designed_codons']).agg(AA=('aa_variant_id','size'),mean_support_fraction=('support_fraction','mean'),zero_support_AA=('T0_HC_codons',lambda x:int((x==0).sum()))).reset_index().to_csv(attempt/'SUPPORT_MULTIPLICITY_DESCRIPTION.csv',index=False)
    np.savez_compressed(attempt/'SUPPORT_CONDITIONAL_DRAWS.npz',positions=draws,**reps)
    atomic_json(attempt/'HS2_SUPPORT_SCORE_GATE.json',gate)
    return [p.name for p in attempt.iterdir() if p.name!='started.json']
