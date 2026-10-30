"""Locked signed position summaries and conditional position resampling."""
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def signed_rho(x,y):
    x,y=np.asarray(x,dtype=float),np.asarray(y,dtype=float)
    if len(x)<3:return np.nan,'fewer_than_three_AA'
    if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Shared finite mask required')
    if np.ptp(x)==0:return np.nan,'constant_predictor'
    if np.ptp(y)==0:return np.nan,'constant_outcome'
    return float(np.corrcoef(rankdata(x,method='average'),rankdata(y,method='average'))[0,1]),'estimable'


def position_records(data,endpoint,predictors):
    rows=[]
    for pos,g in data.groupby('aa_position',sort=True):
        for label,col in predictors.items():
            rho,reason=signed_rho(g[col],g[endpoint])
            rows.append({'endpoint':endpoint,'aa_position':int(pos),'hotspot':g.hotspot.iloc[0],'model':label,'n_AA':len(g),'rho':rho,'reason':reason})
    return pd.DataFrame(rows)


def paired_positions(records,models):
    pivot=records.pivot(index='aa_position',columns='model',values='rho').reindex(columns=models)
    return pivot.loc[np.isfinite(pivot).all(axis=1)].index.to_numpy(dtype=int)


def draw_positions(positions,hotspots,B=10000,seed=20261018):
    positions=np.asarray(sorted(positions),dtype=int)
    if not len(positions):return np.empty((B,0),dtype=int)
    rng=np.random.default_rng(seed);parts=[]
    for hs in sorted(set(hotspots[p] for p in positions)):
        ids=np.asarray([p for p in positions if hotspots[p]==hs],dtype=int)
        parts.append(rng.choice(ids,size=(B,len(ids)),replace=True))
    return np.concatenate(parts,axis=1)


def rank_summaries(records,models,contrasts,draws):
    positions=paired_positions(records,models)
    pivot=records.pivot(index='aa_position',columns='model',values='rho')
    lookup={int(p):i for i,p in enumerate(pivot.index)}
    output=[];replicates={}
    for label,terms in [(m,(m,)) for m in models]+[(a+' minus '+b,(a,b)) for a,b in contrasts]:
        values=pivot[terms[0]].to_numpy(dtype=float)
        if len(terms)==2:values=values-pivot[terms[1]].to_numpy(dtype=float)
        point=float(np.median([values[lookup[p]] for p in positions])) if len(positions) else np.nan
        if draws.shape[1]:
            if set(draws.ravel())-set(positions):raise ValueError('Draw contains ineligible position')
            idx=np.vectorize(lookup.__getitem__)(draws)
            distribution=np.median(values[idx],axis=1)
            lo,hi=np.quantile(distribution,[.025,.975],method='linear')
        else:distribution=np.full(len(draws),np.nan);lo=hi=np.nan
        replicates[label]=distribution
        output.append({'contrast_or_model':label,'point':point,'conditional_low':lo,'conditional_high':hi,'n_positions':len(positions),
                       'interval_type':'95% conditional hotspot-stratified position-resampling interval; fixed records, not generalization CI','B':len(draws),'seed':20261018})
    return output,replicates


def prediction_metrics(data,endpoint,model):
    y=data[endpoint].to_numpy(dtype=float);prediction=data[model+'_OOF'].to_numpy(dtype=float)
    baseline=data[model+'_train_mean_baseline'].to_numpy(dtype=float)
    if not np.isfinite(prediction).all():raise ValueError('Nonfinite OOF predictions')
    errors=(y-prediction)**2;den=(y-baseline)**2
    table=pd.DataFrame({'position':data.aa_position,'error':errors,'baseline_error':den,'absolute_error':np.abs(y-prediction)})
    g=table.groupby('position').mean(numeric_only=True)
    ratio=lambda a,b:1-float(a)/float(b) if np.isfinite(b) and b!=0 else np.nan
    return {'MSE':float(errors.mean()),'MAE':float(np.abs(y-prediction).mean()),
            'pooled_Q2':ratio(errors.sum(),den.sum()),'equal_position_MSE':float(g.error.mean()),
            'equal_position_Q2':ratio(g.error.mean(),g.baseline_error.mean())}
