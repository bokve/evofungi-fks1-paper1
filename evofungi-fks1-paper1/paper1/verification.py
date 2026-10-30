"""Clean-process, source-specific reproduction; no historical model refits."""
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from .stages import atomic_json,digest


def run(inputs,attempt,config):
    summary=pd.read_csv(inputs['compare/COMPARISON_SUMMARIES.csv'])
    oof=pd.read_csv(inputs['compare/OOF_PREDICTIONS.csv'])
    checks=[]
    for (scenario,endpoint),data in oof.groupby(['scenario','endpoint']):
        for model in ['GEMME_ESM','GEMME_ESM_Evo2']:
            err=(data[endpoint]-data[model+'_OOF'])**2
            baseline=(data[endpoint]-data[model+'_train_mean_baseline'])**2
            q=1-err.sum()/baseline.sum() if baseline.sum()!=0 else np.nan
            expected=summary[(summary.scenario==scenario)&(summary.endpoint==endpoint)&(summary.contrast_or_model==model)&(summary.metric=='pooled_Q2')].point.iloc[0]
            np.testing.assert_allclose(q,expected,rtol=0,atol=1e-10,equal_nan=True)
            checks.append({'check':'independent_saved_OOF_pooled_Q2','scenario':scenario,'endpoint':endpoint,'model':model,'actual':float(q),'expected':float(expected)})
        differences=[]
        for position,g in data.groupby('aa_position'):
            if len(g)<3 or g[endpoint].nunique()<2 or g.GEMME_ESM_OOF.nunique()<2 or g.GEMME_ESM_Evo2_OOF.nunique()<2:continue
            a=float(spearmanr(g.GEMME_ESM_Evo2_OOF,g[endpoint]).statistic)
            b=float(spearmanr(g.GEMME_ESM_OOF,g[endpoint]).statistic)
            differences.append(a-b)
        point=float(np.median(differences)) if differences else np.nan
        expected=summary[(summary.scenario==scenario)&(summary.endpoint==endpoint)&(summary.analysis=='supervised-LOPO')&(summary.metric=='median_paired_position_rho_difference')].point.iloc[0]
        np.testing.assert_allclose(point,expected,rtol=0,atol=1e-10,equal_nan=True)
        checks.append({'check':'independent_saved_OOF_median_paired_rho','scenario':scenario,'endpoint':endpoint,'actual':point,'expected':float(expected)})
    arrays=np.load(inputs['compare/CONDITIONAL_RESAMPLING_SUMMARIES.npz'],allow_pickle=False)
    interval_checks=0
    for row in summary[np.isfinite(summary.conditional_low)].itertuples():
        key=row.scenario+'/'+row.endpoint+'/'+row.analysis+'/'+row.contrast_or_model
        if key not in arrays:raise ValueError('Missing exact conditional draw summaries: '+key)
        low,high=np.quantile(arrays[key],[.025,.975],method='linear')
        np.testing.assert_allclose([low,high],[row.conditional_low,row.conditional_high],rtol=0,atol=1e-10)
        interval_checks+=1
    raw=pd.read_csv(inputs['panels/RAW435.csv']).set_index('aa_variant_id')
    replacement=pd.read_csv(inputs['panels/AUTHOR_REPLACEMENT435.csv']).set_index('aa_variant_id').loc[raw.index]
    for col in ['fitness_none','fitness_anidulafungin','fitness_caspofungin','fitness_micafungin']:
        changed=set(raw.index[raw[col]!=replacement[col]])
        if changed!={'FKS1:V641W','FKS1:L642K'}:raise ValueError('Replacement change mask differs')
    members=pd.read_csv(inputs['esm-score/ESM_MEMBER_SCORES.csv'])
    ensemble=pd.read_csv(inputs['esm-score/ESM_ENSEMBLE.csv']).set_index('aa_variant_id')
    pivot=members.pivot(index='aa_variant_id',columns='checkpoint',values='score')
    if pivot.shape!=(532,5) or not np.isfinite(pivot.to_numpy()).all():raise ValueError('Incomplete actual ESM ensemble')
    np.testing.assert_allclose(pivot.mean(axis=1),ensemble.loc[pivot.index,'esm1v_ensemble'],rtol=0,atol=1e-10)
    def independent_ranks(frame,endpoint,predictors):
        values=[]
        for position,g in frame.groupby('aa_position'):
            if len(g)<3 or g[endpoint].nunique()<2 or any(g[col].nunique()<2 for col in predictors.values()):continue
            values.append({label:float(spearmanr(g[col],g[endpoint]).statistic) for label,col in predictors.items()})
        return pd.DataFrame(values,columns=list(predictors))
    zero_checks=0
    predictors={'Evo2':'evo2_primary_AA_score','GEMME':'gemme_primary','ESM':'esm1v_ensemble'}
    for scenario in summary.scenario.unique():
        panel='AUTHOR_REPLACEMENT435' if scenario.startswith('AUTHOR') else 'RAW435'
        frame=pd.read_csv(inputs['panels/'+panel+'.csv']).merge(ensemble[['esm1v_ensemble']],left_on='aa_variant_id',right_index=True,validate='one_to_one')
        if scenario.endswith('_NO_HS2'):frame=frame.loc[frame.hotspot!='HS2']
        for endpoint in summary.endpoint.unique():
            frame_used=frame.loc[np.isfinite(frame[[endpoint,*predictors.values()]]).all(axis=1)]
            ranks=independent_ranks(frame_used,endpoint,predictors)
            rows=summary[(summary.scenario==scenario)&(summary.endpoint==endpoint)&(summary.analysis=='zero-shot')]
            for row in rows.itertuples():
                terms=row.contrast_or_model.split(' minus ')
                value=ranks[terms[0]] if len(terms)==1 else ranks[terms[0]]-ranks[terms[1]]
                expected=float(value.median()) if len(value) else np.nan
                np.testing.assert_allclose(expected,row.point,rtol=0,atol=1e-10,equal_nan=True)
                if len(ranks)!=row.n_positions:raise ValueError('Independent zero-shot estimability count differs')
                zero_checks+=1
    support_checks=0
    if 'support/SUPPORT_SCORE_ASSOCIATIONS.csv' in inputs:
        panel=pd.read_csv(inputs['support/EXPLORATORY_SUPPORT_AA_PANEL.csv'])
        if len(panel)!=208 or not (panel.T0_HC_codons<=panel.designed_codons).all():raise ValueError('Categorical support universe/count mismatch')
        np.testing.assert_allclose(panel.support_fraction,panel.T0_HC_codons/panel.designed_codons,rtol=0,atol=1e-12)
        ranks=independent_ranks(panel,'support_fraction',{m:m for m in ['Evo2_theoretical_H_median','GEMME','ESM']})
        for row in pd.read_csv(inputs['support/SUPPORT_SCORE_ASSOCIATIONS.csv']).itertuples():
            terms=row.contrast_or_model.split(' minus ');value=ranks[terms[0]] if len(terms)==1 else ranks[terms[0]]-ranks[terms[1]]
            np.testing.assert_allclose(value.median(),row.point,rtol=0,atol=1e-10,equal_nan=True)
            if len(ranks)!=row.n_positions:raise ValueError('Independent categorical estimability count differs')
            support_checks+=1
    folds=json.loads(inputs['compare/FOLD_AUDIT.json'].read_text())
    for fold in folds:
        if set(fold['train_AA_ids'])&set(fold['test_AA_ids']) or fold['heldout_position'] in fold['train_positions']:
            raise ValueError('Fold leakage found in saved audit')
    influence=pd.read_csv(inputs['compare/POSITION_DELETION_INFLUENCE.csv'])
    for (scenario,endpoint),g in influence.groupby(['scenario','endpoint']):
        if set(g.deleted_position)!=set(raw.aa_position):raise ValueError('Position deletion coverage incomplete')
        if not g.training_rebuilt.all():raise ValueError('Deletion must rebuild training')
    actual_weight_repeat=json.loads(inputs['esm-finalize/ESM_REAL_REPRODUCTION.json'].read_text())
    if actual_weight_repeat['status']!='EXECUTED_VALIDATED':raise ValueError('Actual post-inference checkpoint reproduction required')
    acquisition=json.loads(inputs['esm-score/CHECKPOINT_ACQUISITION_MANIFEST.json'].read_text())
    checkpoint=Path(acquisition[0]['remote_cache_path'])
    if checkpoint.is_file():
        import torch
        if torch.cuda.is_available():
            from .esm_execution import load_checkpoint
            from .esm_contract import read_protein
            from .esm_scoring import score_masked_position
            model,alphabet,info=load_checkpoint(acquisition[0]['checkpoint'])
            if info['sha256']!=acquisition[0]['sha256']:raise ValueError('Cached checkpoint identity changed')
            seq=read_protein(inputs['protein'],config['protein_record_id']);p=int(raw.aa_position.min())
            values,_=score_masked_position(model,alphabet,seq,p)
            selected=members[(members.checkpoint==info['checkpoint'])&(members.aa_position==p)]
            np.testing.assert_allclose([values[r.mut_aa] for r in selected.itertuples()],selected.score.to_numpy(),rtol=1e-5,atol=1e-5)
            actual_weight_repeat={'status':'EXECUTED_VALIDATED','checkpoint':info['checkpoint'],'sha256':info['sha256'],'position':p,'AA_scores_repeated':len(selected),'atol':1e-5,'rtol':1e-5}
            del model;torch.cuda.empty_cache()
    code_root=Path(__file__).resolve().parents[1]
    tests=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(code_root/'tests'),'-v'],cwd=code_root,capture_output=True,text=True)
    (attempt/'tests.stdout.txt').write_text(tests.stdout);(attempt/'tests.stderr.txt').write_text(tests.stderr)
    if tests.returncode:raise RuntimeError('Clean-process tests failed; see preserved stderr')
    visual=json.loads(inputs['visual_review'].read_text())
    figure_index=json.loads(inputs['figures/FIGURE_INDEX.json'].read_text())
    if set(visual['inspected_figures'])!={x['figure'] for x in figure_index} or visual['status']!='INSPECTED':
        raise ValueError('Actual visual inspection of every rendered figure required')
    for fig in figure_index:
        png=inputs['figures/'+fig['figure']+'.png']
        if visual['PNG_sha256'][fig['figure']]!=digest(png):raise ValueError('Reviewed figure changed')
        fig['visual_inspection']='INSPECTED_HASH_VERIFIED'
    atomic_json(attempt/'FIGURE_INDEX_VERIFIED.json',figure_index)
    manuscript=inputs['writeup/MANUSCRIPT.md'].read_text()
    for forbidden in ['TODO_SCIENCE','Evo2 proves clinical resistance','H3 confidence intervals','all analyses were preregistered']:
        if forbidden in manuscript:raise ValueError('Unsupported manuscript assertion: '+forbidden)
    for check in checks:
        for key in ['actual','expected']:
            if key in check and not np.isfinite(check[key]):check[key]=None
    report={'status':'EXECUTED_VALIDATED','scope':'clean-process saved-metric/interval/ensemble/fold/deletion/source checks and selected cached checkpoint repeat; not full historical/model-inference rerun',
       'selected_metric_checks':checks,'zero_shot_point_checks':zero_checks,'categorical_support_point_checks':support_checks,'conditional_interval_checks':interval_checks,'actual_checkpoint_repeat':actual_weight_repeat,
       'clean_tests_returncode':tests.returncode,'source_fixture_changed_input_invalidation':'covered by clean test_stages/test_pipeline_cli; no source inputs mutated',
       'visual_review':visual,'source_raw_replacement_difference_gate':True,'negative_effects_allowed':True,'original_Evo2_API_recalled':False,
       'hosted_Evo2_immutable_SHA':'UNVERIFIED','claim_source_ledger':str(inputs['writeup/CLAIM_LEDGER.csv'])}
    atomic_json(attempt/'REPRODUCIBILITY_REPORT.json',report)
    return ['REPRODUCIBILITY_REPORT.json','FIGURE_INDEX_VERIFIED.json','tests.stdout.txt','tests.stderr.txt']
