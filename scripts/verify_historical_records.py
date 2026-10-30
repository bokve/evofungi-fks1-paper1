"""Verify archived R diagnostics and H3 quantiles; never refit historical models."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy.stats import chi2,kurtosis
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',type=Path,required=True);p.add_argument('--out',type=Path,default=ROOT/'historical_verification');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
r=a.data_root/'results';legacy=a.data_root/'legacy';h=legacy/'05_results/paper1_hierarchical'
qc=pd.read_csv(h/'P1_NEW2C0A_F0_ML_FIT_QC_v1.1.csv').set_index('model')
assert qc.hard_numeric_QC_PASS.all() and qc.n_warnings.eq(0).all()
lr=2*(qc.loc['FULL_ML_F0','logLik']-qc.loc['REDUCED_ML_F0','logLik'])
reported=pd.read_csv(r/'history/NEW2_summary.csv');rr=reported[reported.result_type=='primary_omnibus_LRT'].iloc[0]
np.testing.assert_allclose([lr,chi2.sf(lr,3)],[rr.estimate,rr.formal_p_value],atol=1e-10,rtol=0)
res=pd.read_csv(h/'P1_NEW2C2C1_ROW_LEVEL_DIAGNOSTICS_RAW_3884_v1.0.csv').normalized_residual_level3
q=pd.read_csv(r/'history/H3_Q2.csv');d=pd.read_csv(r/'history/H3_deltas.csv');paired=pd.read_csv(r/'history/h3_paired.csv')
components=['total','codon_local','downstream50'];models=['M0_G0','M1_G0_plus_G1','M2_G0_plus_G2','M3_G0_plus_G1_plus_G2'];contrasts={'M1-M0':(1,0),'M2-M0':(2,0),'M3-M0':(3,0),'M3-M1':(3,1),'M3-M2':(3,2)};checks=[]
def check(kind,key,values,row):
 for name,value in zip(['sensitivity_low_025','sensitivity_median','sensitivity_high_975'],np.quantile(values,[.025,.5,.975])):
  diff=abs(float(value)-float(row[name]));assert diff<1e-12
  checks.append({'kind':kind,'key':key,'summary':name,'recomputed':float(value),'reported':float(row[name]),'absolute_difference':diff})
with np.load(legacy/'05_results/P1_NEW3_codon_attribution/P1_NEW3H3_DIRICHLET_REPLICATES_AND_WEIGHTS_v1.0.npz',allow_pickle=False) as ar:
 for row in q.to_dict('records'):
  arr=ar[row['panel']+'_Q2'];assert arr.shape==(10000,3,4) and np.isfinite(arr).all()
  check('Q2',row['panel']+':'+row['component']+':'+row['model'],arr[:,components.index(row['component']),models.index(row['model'])],row)
 for row in d.to_dict('records'):
  arr=ar[row['panel']+'_Q2'];c=components.index(row['component']);j,k=contrasts[row['contrast']]
  check('delta',row['panel']+':'+row['component']+':'+row['contrast'],arr[:,c,j]-arr[:,c,k],row)
 for row in paired.to_dict('records'):
  x=ar['COMMON124_COMPLETE_Q2'];y=ar['OBSERVED278_Q2'];c=components.index(row['component'])
  if row['metric']=='Q2':
   j=models.index(row['model_or_contrast']);v=x[:,c,j]-y[:,c,j]
  else:
   j,k=contrasts[row['model_or_contrast']];v=x[:,c,j]-x[:,c,k]-y[:,c,j]+y[:,c,k]
  check('paired',row['component']+':'+row['metric']+':'+row['model_or_contrast'],v,row)
 weight_keys=[k for k in ar.files if k.endswith('_weights')]
 for k in weight_keys:assert np.isfinite(ar[k]).all() and (ar[k]>0).all()
 assert np.array_equal(ar['OBSERVED278_weights'],ar['COMMON124_COMPLETE_weights'])
pd.DataFrame(checks).to_csv(a.out/'H3_QUANTILE_CHECKS.csv',index=False)
report={'status':'PASS','scope':'Saved-object diagnostics and quantile arithmetic, not independent R or H3 optimization','hierarchical_LR':float(lr),'df':3,'hierarchical_p':float(chi2.sf(lr,3)),'residual_rows':len(res),'abs_normalized_residual_over_3':int((res.abs()>3).sum()),'excess_kurtosis':float(kurtosis(res,fisher=True,bias=True)),'H3_statistics_checked':len(checks),'H3_max_difference':max(v['absolute_difference'] for v in checks),'historical_fits_reexecuted':False}
(a.out/'HISTORICAL_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
