"""All planned signed comparisons, fixed LOPO and bounded sensitivities."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .comparison import matched_lopo
from .metrics import signed_rho,position_records,paired_positions,draw_positions,rank_summaries,prediction_metrics
from .panels import ENDPOINTS,DRUGS
from .stages import atomic_json,digest

PREDICTORS={'Evo2':'evo2_primary_AA_score','GEMME':'gemme_primary','ESM':'esm1v_ensemble'}
OOF={'GEMME_ESM':'GEMME_ESM_OOF','GEMME_ESM_Evo2':'GEMME_ESM_Evo2_OOF'}
CONTRASTS=[('ESM','GEMME'),('ESM','Evo2'),('Evo2','GEMME')]


def analyse(frame,endpoint,uncertainty=True):
    cols=[*PREDICTORS.values(),endpoint]
    keep=np.isfinite(frame[cols].to_numpy(dtype=float)).all(axis=1)
    data=frame.loc[keep].copy().sort_values('aa_variant_id').reset_index(drop=True)
    if data.aa_variant_id.duplicated().any():raise ValueError('AA duplicates in comparison')
    mask=data[['aa_variant_id','aa_position','hotspot']].copy()
    mask['endpoint']=endpoint
    mask['shared_mask_sha256']=__import__('hashlib').sha256('\n'.join(data.aa_variant_id).encode()).hexdigest()
    oof,audit,dropped=matched_lopo(data,endpoint)
    if dropped or set(oof.aa_variant_id)!=set(data.aa_variant_id):raise ValueError('Nested OOF comparison mask changed')
    # matched_lopo deliberately returns only required columns; keyed metadata join.
    oof=oof.merge(data[['aa_variant_id','hotspot']],on='aa_variant_id',validate='one_to_one')
    for record in audit:
        if set(record['train_AA_ids'])&set(record['test_AA_ids']) or record['heldout_position'] in record['train_positions']:
            raise ValueError('Cross-fold leakage')
    zero=position_records(data,endpoint,PREDICTORS)
    supervised=position_records(oof,endpoint,OOF)
    hs=data.groupby('aa_position').hotspot.first().to_dict()
    draws={};summaries=[];replicates={}
    for kind,records,models,contrasts in [('zero-shot',zero,list(PREDICTORS),CONTRASTS),('supervised-LOPO',supervised,list(OOF),[('GEMME_ESM_Evo2','GEMME_ESM')])]:
        eligible=paired_positions(records,models)
        d=draw_positions(eligible,hs,B=10000 if uncertainty else 0)
        draws[kind]=d
        # Influence calls use only point summaries; no fake B=0 percentile.
        if uncertainty:
            s,reps=rank_summaries(records,models,contrasts,d)
            replicates.update({kind+'/'+k:v for k,v in reps.items()})
        else:
            pivot=records.pivot(index='aa_position',columns='model',values='rho').loc[eligible]
            s=[{'contrast_or_model':m,'point':float(pivot[m].median()) if len(pivot) else np.nan,'n_positions':len(eligible)} for m in models]
            s += [{'contrast_or_model':a+' minus '+b,'point':float((pivot[a]-pivot[b]).median()) if len(pivot) else np.nan,'n_positions':len(eligible)} for a,b in contrasts]
        for r in s:r.update(endpoint=endpoint,analysis=kind,n_AA=len(data),metric='median_signed_position_rho' if ' minus ' not in r['contrast_or_model'] else 'median_paired_position_rho_difference')
        summaries.extend(s)
    for model in OOF:
        metrics=prediction_metrics(oof,endpoint,model)
        for metric,value in metrics.items():
            summaries.append({'endpoint':endpoint,'analysis':'supervised-LOPO','contrast_or_model':model,'metric':metric,'point':value,'n_AA':len(data),'n_positions':data.aa_position.nunique()})
    for metric in ['pooled_Q2','equal_position_Q2']:
        a=prediction_metrics(oof,endpoint,'GEMME_ESM_Evo2')[metric];b=prediction_metrics(oof,endpoint,'GEMME_ESM')[metric]
        summaries.append({'endpoint':endpoint,'analysis':'supervised-LOPO','contrast_or_model':'GEMME_ESM_Evo2 minus GEMME_ESM','metric':'delta_'+metric,'point':a-b,'n_AA':len(data),'n_positions':data.aa_position.nunique()})
    for model,col in PREDICTORS.items():
        pooled,reason=signed_rho(data[col],data[endpoint])
        summaries.append({'endpoint':endpoint,'analysis':'pooled_descriptive','contrast_or_model':model,'metric':'pooled_Spearman','point':pooled,'n_AA':len(data),'n_positions':data.aa_position.nunique(),'reason':reason})
    return {'summaries':pd.DataFrame(summaries),'zero_positions':zero,'OOF_positions':supervised,'OOF':oof,
            'audit':audit,'mask':mask,'draws':draws,'replicates':replicates,'excluded_AA':frame.loc[~keep,'aa_variant_id'].tolist()}


def run(inputs,attempt,config):
    ensemble=pd.read_csv(inputs['esm-score/ESM_RAW435.csv'])
    scenarios={}
    for panel in ['RAW435','AUTHOR_REPLACEMENT435']:
        frame=pd.read_csv(inputs['panels/'+panel+'.csv'])
        if set(frame.aa_variant_id)!=set(ensemble.aa_variant_id):raise ValueError('ESM/phenotype AA masks differ')
        frame=frame.merge(ensemble[['aa_variant_id','esm1v_ensemble']],on='aa_variant_id',validate='one_to_one')
        scenarios[panel]=frame
        scenarios[panel+'_NO_HS2']=frame.loc[frame.hotspot!='HS2'].copy()
    all_summary=[];all_zero=[];all_oofpos=[];all_pred=[];all_masks=[];audits=[];draws={};replicates={};influence=[]
    results={}
    for scenario,frame in scenarios.items():
        for endpoint in ENDPOINTS:
            result=analyse(frame,endpoint);results[(scenario,endpoint)]=result
            for key,store in [('summaries',all_summary),('zero_positions',all_zero),('OOF_positions',all_oofpos),('OOF',all_pred),('mask',all_masks)]:
                f=result[key].copy();f['scenario']=scenario;f['endpoint']=endpoint;store.append(f)
            audits.extend([dict(r,scenario=scenario,endpoint=endpoint) for r in result['audit']])
            for k,v in result['draws'].items():draws[scenario+'/'+endpoint+'/'+k]=v
            for k,v in result['replicates'].items():replicates[scenario+'/'+endpoint+'/'+k]=v
        print(json.dumps({'comparison_scenario_complete':scenario}),flush=True)
    for panel in ['RAW435','AUTHOR_REPLACEMENT435']:
        frame=scenarios[panel]
        for position in sorted(frame.aa_position.unique()):
            reduced=frame.loc[frame.aa_position!=position]
            for endpoint in ENDPOINTS:
                deleted_result=analyse(reduced,endpoint,uncertainty=False)
                if any(int(position) in a['train_positions'] or a['heldout_position']==int(position) for a in deleted_result['audit']):
                    raise ValueError('Deleted position leaked into rebuilt LOPO')
                r=deleted_result['summaries']
                selected=r[(r.metric=='median_paired_position_rho_difference') | r.metric.isin(['delta_pooled_Q2','delta_equal_position_Q2'])]
                for row in selected.to_dict('records'):influence.append(dict(row,scenario=panel,deleted_position=int(position),training_rebuilt=True))
        print(json.dumps({'all_position_deletion_complete':panel}),flush=True)
    changes=[];paired_draw_notes=[]
    for suffix in ['', '_NO_HS2']:
        for endpoint in ENDPOINTS:
            a=results[('RAW435'+suffix,endpoint)]['summaries'];b=results[('AUTHOR_REPLACEMENT435'+suffix,endpoint)]['summaries']
            keys=['analysis','contrast_or_model','metric'];joined=a.merge(b,on=keys,suffixes=('_raw','_replacement'),validate='one_to_one')
            joined['point_change']=joined.point_replacement-joined.point_raw
            joined['endpoint']=endpoint;joined['sensitivity']='AUTHOR_REPLACEMENT'+suffix;changes.extend(joined.to_dict('records'))
            for kind,models,contrasts,key in [('zero-shot',list(PREDICTORS),CONTRASTS,'zero_positions'),('supervised-LOPO',list(OOF),[('GEMME_ESM_Evo2','GEMME_ESM')],'OOF_positions')]:
                rawrec=results[('RAW435'+suffix,endpoint)][key];reprec=results[('AUTHOR_REPLACEMENT435'+suffix,endpoint)][key]
                common=sorted(set(paired_positions(rawrec,models))&set(paired_positions(reprec,models)))
                hs=scenarios['RAW435'+suffix].groupby('aa_position').hotspot.first().to_dict()
                shared=draw_positions(common,hs)
                _,ar=rank_summaries(rawrec[rawrec.aa_position.isin(common)],models,contrasts,shared)
                _,br=rank_summaries(reprec[reprec.aa_position.isin(common)],models,contrasts,shared)
                for label in ar:
                    delta=br[label]-ar[label];finite=delta[np.isfinite(delta)]
                    lo,hi=np.quantile(finite,[.025,.975]) if len(finite) else [np.nan,np.nan]
                    paired_draw_notes.append({'endpoint':endpoint,'sensitivity':suffix or 'FULL','analysis':kind,'contrast_or_model':label,'paired_positions':len(common),'conditional_change_low':lo,'conditional_change_high':hi,'interval_type':'conditional paired position-resampling change; fixed OOF records'})
                    replicates['REPLACEMENT_CHANGE'+suffix+'/'+endpoint+'/'+kind+'/'+label]=delta
                draws['PAIRED_REPLACEMENT'+suffix+'/'+endpoint+'/'+kind]=shared
    total_drug=[]
    for scenario,frame in scenarios.items():
        for drug in DRUGS:
            endpoint='fitness_'+drug
            keep=np.isfinite(frame[[*PREDICTORS.values(),endpoint]].to_numpy()).all(axis=1)
            data=frame.loc[keep]
            records=position_records(data,endpoint,PREDICTORS)
            positions=paired_positions(records,list(PREDICTORS))
            for model,col in PREDICTORS.items():
                rho,_=signed_rho(data[col],data[endpoint])
                total_drug.append({'scenario':scenario,'endpoint':endpoint,'model':model,'pooled_Spearman':rho,
                   'median_paired_position_rho':records[(records.model==model)&records.aa_position.isin(positions)].rho.median(),
                   'n_AA':len(data),'n_paired_positions':len(positions),'role':'total-drug-fitness descriptive; not primary D-C comparison'})
    files={
      'COMPARISON_SUMMARIES.csv':pd.concat(all_summary,ignore_index=True),
      'ZERO_SHOT_POSITION_RESULTS.csv':pd.concat(all_zero,ignore_index=True),
      'LOPO_POSITION_RESULTS.csv':pd.concat(all_oofpos,ignore_index=True),
      'OOF_PREDICTIONS.csv':pd.concat(all_pred,ignore_index=True),
      'COMPARISON_MASK_LEDGER.csv':pd.concat(all_masks,ignore_index=True),
      'POSITION_DELETION_INFLUENCE.csv':pd.DataFrame(influence),
      'METRIC_CHANGE_LEDGER.csv':pd.DataFrame(changes),
      'AUTHOR_REPLACEMENT_PAIRED_UNCERTAINTY.csv':pd.DataFrame(paired_draw_notes),
      'TOTAL_DRUG_FITNESS_DESCRIPTIVE.csv':pd.DataFrame(total_drug)}
    for name,frame in files.items():frame.to_csv(attempt/name,index=False)
    atomic_json(attempt/'FOLD_AUDIT.json',audits)
    np.savez_compressed(attempt/'CONDITIONAL_RESAMPLING_DRAWS.npz',**draws)
    np.savez_compressed(attempt/'CONDITIONAL_RESAMPLING_SUMMARIES.npz',**replicates)
    atomic_json(attempt/'COMPARISON_INTERPRETATION.json',{'primary_endpoint':'fitness_none','secondary_endpoints':ENDPOINTS[1:],
       'new_p_values':False,'interval':'Conditional hotspot-stratified position-resampling interval on observed panel/fixed OOF; not generalization CI; no upstream measurement uncertainty.',
       'overlapping_training_folds':True,'all_position_deletions_refit':True,'HS2_removed_from_train_and_test_in_sensitivity':True,
       'no_outcome_dependent_model_or_mask_choice':True,'negative_effects_allowed':True})
    return [*files,'FOLD_AUDIT.json','CONDITIONAL_RESAMPLING_DRAWS.npz','CONDITIONAL_RESAMPLING_SUMMARIES.npz','COMPARISON_INTERPRETATION.json']
