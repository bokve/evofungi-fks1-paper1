"""Integrate existing evidence without refitting historical analyses."""
import json
import re
import pandas as pd
from .stages import atomic_json, digest


def hs2_core(inputs, attempt, config):
    panel = pd.read_csv(inputs['ascertainment'])
    multiplicity = pd.read_csv(inputs['multiplicity'])
    if len(panel)!=814 or panel.exact_genotype_key.duplicated().any():
        raise ValueError('NNK exact-genotype universe mismatch')
    for col in ['PRE_HC_support','T0_HC_support','FINAL589_inclusion']:
        if not set(panel[col]).issubset({0,1}): raise ValueError('Invalid support indicator')
    if not (panel.T0_HC_support==panel.FINAL589_inclusion).all():
        raise ValueError('T0-HC/FINAL identity changed')
    if panel.PRE_HC_support.sum()!=811 or panel.T0_HC_support.sum()!=589 or len(multiplicity)!=532:
        raise ValueError('Frozen support universe changed')
    starts=panel.groupby('hotspot').canonical_aa_position.min().to_dict()
    if not (panel.canonical_aa_position == panel.hotspot.map(starts)+panel.local_codon_index_0based).all():
        raise ValueError('Canonical sequence-derived coordinate mismatch')
    grouped=panel.groupby('aa_group_key').agg(n=('exact_genotype_key','size'),pre=('PRE_HC_support','sum'),t0=('T0_HC_support','sum'))
    m=multiplicity.set_index('aa_group_key')
    if set(m.index)!=set(grouped.index):raise ValueError('AA support mask mismatch')
    for a,b in [('n','n_theoretical_codons'),('pre','n_PRE_HC_codons'),('t0','n_T0_HC_codons')]:
        if not (grouped.loc[m.index,a]==m[b]).all():raise ValueError('Multiplicity lineage mismatch')
    panel['source_path']=str(inputs['ascertainment']);panel['source_sha256']=digest(inputs['ascertainment'])
    panel.to_csv(attempt/'HS2_CORE_MEMBERSHIP.csv',index=False)
    multiplicity['source_path']=str(inputs['multiplicity']);multiplicity['source_sha256']=digest(inputs['multiplicity'])
    multiplicity.to_csv(attempt/'HS2_AA_MULTIPLICITY.csv',index=False)
    summary=panel.groupby(['hotspot','canonical_aa_position']).agg(designed_codons=('exact_genotype_key','size'),PRE_HC=('PRE_HC_support','sum'),T0_HC_FINAL=('T0_HC_support','sum'),AA_designed=('aa_group_key','nunique')).reset_index()
    for col,name in [('PRE_HC_support','AA_PRE_HC'),('T0_HC_support','AA_T0_HC')]:
        counts=panel.loc[panel[col]==1].groupby(['hotspot','canonical_aa_position']).aa_group_key.nunique()
        summary[name]=[int(counts.get((r.hotspot,r.canonical_aa_position),0)) for r in summary.itertuples()]
    informative=multiplicity.loc[multiplicity.n_T0_HC_codons>=2].groupby(['hotspot','canonical_aa_position']).size()
    summary['AA_T0_multi_codon']=[int(informative.get((r.hotspot,r.canonical_aa_position),0)) for r in summary.itertuples()]
    summary.to_csv(attempt/'HS2_CORE_SOURCE_TABLE.csv',index=False)
    atomic_json(attempt/'L1360_UNIVERSE_NOTE.json',{'position':1360,'theoretical_design_present':bool((panel.canonical_aa_position==1360).any()),
        'observed_informative_same_AA_panel':'absent by group multiplicity, not by phenotype/model-performance exclusion',
        'theoretical_H_panel':'present; distinguish 28-position H panel from observed 27-position same-AA panel',
        'support_interpretation':'High-confidence source representation; not death or absolute survival. FINAL589 is T0-HC, not an additional time stage.'})
    requirements=[]
    for role,path in inputs.items():
        if not role.startswith('frequency_'):continue
        # Bounded metadata inspection only: no frequency values, ratios or associations.
        fields=pd.read_csv(path,nrows=0).columns.tolist()
        unit=[x for x in ['Sample attributes','Sample_name','Replicate','Timepoint'] if x in fields]
        meta=pd.read_csv(path,usecols=unit).drop_duplicates()
        counts=[x for x in fields if re.search(r'count|denominator|total_reads',x,re.I)]
        requirements.append({'role':role,'source_path':str(path),'sha256':digest(path),'metadata':meta.to_dict('records'),'count_denominator_columns':counts})
    atomic_json(attempt/'HS2_FEASIBILITY_GATE.json',{'status':'BLOCKED_INPUT','scope':'one bounded read of already located exports and pinned source configuration; no ratio computed',
       'requirements':{'material_transformation_relationship':'not documented by existing source metadata',
        'replicate_pairing':'PRE/mainrun labels do not establish paired biological material',
        'denominator_count_provenance':'raw count and denominator fields absent from selected frequency exports',
        'platform_library_normalization':'cannot establish comparable units; PRE-HS2 metadata anomaly retained',
        'zero_support':'source confidence labels do not identify absolute survival'},
       'sources':requirements,'configuration_source':str(inputs['author_config']),
       'configuration_sha256':digest(inputs['author_config']),
       'no_additional_pseudocount':True,'categorical_score_support_independent_gate':True})
    return [p.name for p in attempt.iterdir() if p.name!='started.json']


def history(inputs, attempt, config):
    records=[]
    for role,path in inputs.items():
        sha=digest(path)
        if path.suffix=='.csv':
            data=pd.read_csv(path).fillna('')
            for i,row in data.iterrows():
                records.append({'role':role,'source_path':str(path),'sha256':sha,'row_or_pointer':str(i+2),'record':row.to_dict(),'evidence_status':'REUSED_VALIDATED_RECORD_NOT_REFIT'})
        elif path.suffix=='.json':
            records.append({'role':role,'source_path':str(path),'sha256':sha,'row_or_pointer':'','record':json.loads(path.read_text()),'evidence_status':'REUSED_VALIDATED_RECORD_NOT_REFIT'})
    atomic_json(attempt/'HISTORICAL_SOURCE_LEDGER.json',records)
    flat=pd.DataFrame([dict(r,record=json.dumps(r['record'],sort_keys=True)) for r in records])
    flat.to_csv(attempt/'NEW3_PAPER1_INTEGRATION_TABLE.csv',index=False)
    for role in ['NEW2_summary','H2_points','H3_deltas','H3_Q2','h3_paired','H4B_claims','M5_robustness','NEW1_primary','benchmark_pooled','benchmark_positions','method1_summary','method1_components','method2_summary','method3_summary','method3_reliability','structure_descriptors']:
        pd.read_csv(inputs[role]).to_csv(attempt/(role+'.csv'),index=False)
    (attempt/'NEW3_SCOPE_AMENDMENT.md').write_text('''# NEW-3 Paper 1 integration — 18 October 2026

Original A–H freezes, all M0–M3 and H panels are reused as source records, not newly refitted results. H4B supersedes H4A\'s invalid outcome-dependent success assertion and closes H only. F2 ordinary resampling retained structural failures and valid-only intervals; F2B/H3 support-preserving Dirichlet intervals are sensitivity intervals, not sampling confidence intervals, posterior intervals or p-values. All signs and unselected branches remain in the linked tables.

NEW-3I full-standard-code and separate WT-synonymous synthetic panels: DEFERRED_SCOPE. They answer additional synthetic questions outside the chosen Paper 1 experimental NNK scope; they were not executed. Alternative codons encoding the same mutant AA are not a WT-synonymous DMS panel.

The original H3 summary bytes and audits were validated; this stage does not reconstruct all historical NPZ/10,000-fit arrays or 518 hosted API requests. Hosted immutable Evo2 weight SHA remains UNVERIFIED. Method 3 reliability is reused with its assay/unit limitations; source NT-to-AA aggregation does not establish independent clinical or codon fitness truth.
''')
    atomic_json(attempt/'LIMITED_REPRODUCTION_INVENTORY.json',{'records_validated':len(records),'fresh_historical_refits':False,'H3_arrays_freshly_reconstructed':False,'hosted_API_recalled':False,'hosted_immutable_checkpoint_SHA':'UNVERIFIED',
       'preserved_tests':['NEW1 negative/uncertain','NEW2 primary and rank/aggregate sensitivity','M5 including varying GEMME estimability','Methods1-3','NEW3 M0-M3/F2 failures/F2B/H2/H3/H4B'],
       'status_scope':'chosen Paper1 NEW3 integration complete; synthetic future scope deferred; full historical compute reproduction not claimed'})
    atomic_json(attempt/'OPTIONAL_SOURCE_DECISION.json',{'status':'DEFERRED_SCOPE','scope':'optional separate same-Fks1 validation comparison',
       'reason':'Available author validation-derived records are used for the required replacement source sensitivity. A separate validation comparison is not included in this bounded extension and is not represented as an independent test of those replacements. No claim of global source absence; no source reconstruction campaign launched.'})
    return [p.name for p in attempt.iterdir() if p.name!='started.json']
