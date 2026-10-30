"""CPU replay of already archived statistics. No inference or new predictions."""
from pathlib import Path
import argparse, hashlib, json, platform
import importlib.metadata
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
def rho(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
    return float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
def close(a,b):
    np.testing.assert_allclose(a,b,rtol=0,atol=1e-10,equal_nan=True)
def run(data,out):
    results=data/'results'; out.mkdir(parents=True,exist_ok=True)
    summaries=pd.read_csv(results/'compare/COMPARISON_SUMMARIES.csv')
    oof=pd.read_csv(results/'compare/OOF_PREDICTIONS.csv')
    ensemble=pd.read_csv(results/'esm-score/ESM_ENSEMBLE.csv').set_index('aa_variant_id')
    primary=pd.read_csv(results/'panels/RAW435.csv')
    assert len(primary)==435 and primary.aa_position.nunique()==28
    assert int(primary.n_experimental_codons.sum())==589 and primary.aa_variant_id.is_unique
    replacement=pd.read_csv(results/'panels/AUTHOR_REPLACEMENT435.csv').set_index('aa_variant_id').loc[primary.aa_variant_id]
    for c in ['fitness_none','fitness_anidulafungin','fitness_caspofungin','fitness_micafungin']:
        assert set(primary.aa_variant_id[np.abs(primary[c].to_numpy()-replacement[c].to_numpy())>1e-12])=={'FKS1:V641W','FKS1:L642K'}
    for c in ['evo2_primary_AA_score','gemme_primary']:
        close(primary[c],replacement[c].to_numpy())
    members=pd.read_csv(results/'esm-score/ESM_MEMBER_SCORES.csv')
    pivot=members.pivot(index='aa_variant_id',columns='checkpoint',values='score')
    assert pivot.shape==(532,5) and np.isfinite(pivot.to_numpy()).all()
    close(pivot.mean(axis=1),ensemble.loc[pivot.index,'esm1v_ensemble'])
    predictions=['GEMME_ESM','GEMME_ESM_Evo2']; fixed={'Evo2':'evo2_primary_AA_score','GEMME':'gemme_primary','ESM':'esm1v_ensemble'}
    records=[]; checks=[]; point_count=0; rank_count=0; interval_count=0; max_diff=0.
    with np.load(results/'compare/CONDITIONAL_RESAMPLING_DRAWS.npz',allow_pickle=False) as draws, np.load(results/'compare/CONDITIONAL_RESAMPLING_SUMMARIES.npz',allow_pickle=False) as distributions:
      for (scenario,endpoint),g in oof.groupby(['scenario','endpoint'],sort=True):
        panel='AUTHOR_REPLACEMENT435' if scenario.startswith('AUTHOR') else 'RAW435'
        f=pd.read_csv(results/f'panels/{panel}.csv').merge(ensemble[['esm1v_ensemble']],left_on='aa_variant_id',right_index=True,validate='one_to_one')
        if scenario.endswith('_NO_HS2'):f=f[f.hotspot!='HS2']
        assert set(f.aa_variant_id)==set(g.aa_variant_id)
        for drug in ['anidulafungin','caspofungin','micafungin']:close(f['drug_specific_'+drug],f['fitness_'+drug]-f.fitness_none)
        table={analysis:{} for analysis in ['zero-shot','supervised-LOPO']}
        for pos,h in f.groupby('aa_position',sort=True):
          table['zero-shot'][pos]={m:rho(h[col],h[endpoint]) for m,col in fixed.items()};rank_count+=3
        for pos,h in g.groupby('aa_position',sort=True):
          table['supervised-LOPO'][pos]={m:rho(h[m+'_OOF'],h[endpoint]) for m in predictions}
        selected=summaries[(summaries.scenario==scenario)&(summaries.endpoint==endpoint)]
        for row in selected.itertuples():
          terms=row.contrast_or_model.split(' minus ')
          if row.metric in ['median_signed_position_rho','median_paired_position_rho_difference']:
            t=pd.DataFrame(table[row.analysis]).T.dropna()
            values=t[terms[0]] if len(terms)==1 else t[terms[0]]-t[terms[1]]
            actual=float(values.median());close(actual,row.point);point_count+=1
            assert len(t)==row.n_positions
            if np.isfinite(row.conditional_low):
              key=scenario+'/'+endpoint+'/'+row.analysis
              draw=draws[key];assert set(draw.ravel())<=set(t.index)
              look={p:i for i,p in enumerate(t.index)}
              idx=np.array([[look[p] for p in line] for line in draw],dtype=int)
              dist=np.median(values.to_numpy()[idx],axis=1)
              close(dist,distributions[key+'/'+row.contrast_or_model])
              quant=np.quantile(dist,[.025,.975]);close(quant,[row.conditional_low,row.conditional_high]);interval_count+=1
          else:
            metrics={}
            for m in predictions:
              err=(g[endpoint]-g[m+'_OOF'])**2;den=(g[endpoint]-g[m+'_train_mean_baseline'])**2
              by=pd.DataFrame({'p':g.aa_position,'e':err,'d':den,'a':abs(g[endpoint]-g[m+'_OOF'])}).groupby('p').mean()
              metrics[m]={'pooled_Q2':1-err.sum()/den.sum(),'equal_position_Q2':1-by.e.mean()/by.d.mean(),'MSE':err.mean(),'MAE':abs(g[endpoint]-g[m+'_OOF']).mean(),'equal_position_MSE':by.e.mean(),'pooled_Spearman':rho(g[m+'_OOF'],g[endpoint])}
            metric=row.metric.removeprefix('delta_')
            if metric=='pooled_Spearman' and terms[0] in fixed:
              metrics.update({m:{'pooled_Spearman':rho(f[col],f[endpoint])} for m,col in fixed.items()})
            if metric not in metrics[terms[0]]:raise ValueError('Unrecognized saved metric '+row.metric)
            actual=metrics[terms[0]][metric]
            if len(terms)==2:actual-=metrics[terms[1]][metric]
            close(actual,row.point);point_count+=1
          if np.isfinite(actual) and np.isfinite(row.point):max_diff=max(max_diff,abs(actual-row.point))
          checks.append({'scenario':scenario,'endpoint':endpoint,'model_or_contrast':row.contrast_or_model,'metric':row.metric,'recomputed':actual,'reported':row.point})
        for analysis,t in table.items():
          for p,v in t.items():
            for m,r in v.items():records.append({'scenario':scenario,'endpoint':endpoint,'analysis':analysis,'position':int(p),'model':m,'rho':r})
    fold_audit=json.loads((results/'compare/FOLD_AUDIT.json').read_text())
    for fold in fold_audit:
        assert not set(fold['train_AA_ids'])&set(fold['test_AA_ids'])
        assert fold['heldout_position'] not in fold['train_positions']
    pd.DataFrame(checks).to_csv(out/'RECOMPUTED_SUMMARIES.csv',index=False)
    pd.DataFrame(records).to_csv(out/'RECOMPUTED_POSITION_CORRELATIONS.csv',index=False)
    primary_rows=summaries[(summaries.scenario=='RAW435')&(summaries.endpoint=='fitness_none')]
    primary_rows.to_csv(out/'DRUG_FREE_PRIMARY_RESULT.csv',index=False)
    report={'status':'PASS','scope':'Fixed-score and saved-OOF replay only; no inference, phenotype refitting, R/H3 refitting, or new biological predictions','AA':435,'positions':28,'exact_genotypes':589,'summary_points_checked':point_count,'fixed_score_position_correlations':rank_count,'conditional_intervals_checked':interval_count,'saved_fold_audits_checked':len(fold_audit),'maximum_summary_absolute_difference':max_diff,'python':platform.python_version(),'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy']}}
    (out/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',type=Path,default=ROOT/'data');p.add_argument('--out',type=Path,default=ROOT/'verification_run');a=p.parse_args();run(a.data_root,a.out)
