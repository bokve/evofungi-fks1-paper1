"""Render CBC figure compositions directly from their frozen CSV sources.
Recreated plot layout; numerical values, uncertainty and memberships are unchanged.
"""
from pathlib import Path
import argparse,json,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',type=Path,default=ROOT/'data');p.add_argument('--out',type=Path,default=ROOT/'rebuilt_cbc_figures');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
src=a.data_root/'cbc_figures';index=[]
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
eps=['fitness_none','drug_specific_anidulafungin','drug_specific_caspofungin','drug_specific_micafungin'];epnames=['Control','AF D−C','CF D−C','MF D−C'];colors=['#0072B2','#D55E00','#009E73'];models=['GEMME_ESM','GEMME_ESM_Evo2']
def read(n):return pd.read_csv(src/(n+'_source.csv'))
def save(fig,n):
 fig.tight_layout()
 for ext in ['pdf','svg','png']:fig.savefig(a.out/(n+'.'+ext),dpi=300,bbox_inches='tight')
 plt.close(fig);shutil.copy2(src/(n+'_source.csv'),a.out/(n+'_source.csv'));index.append({'figure':n,'source':n+'_source.csv','method':'frozen source table; recreated CBC layout; no refitting'})
def intervals(ax,rows,point,low,high,labels,color='#0072B2'):
 for i,r in enumerate(rows.to_dict('records')):ax.hlines(i,r[low],r[high],color=color);ax.plot(r[point],i,'o',color=color,ms=4)
 ax.set_yticks(range(len(rows)),labels);ax.invert_yaxis();ax.axvline(0,color='.65',lw=.8)
# Figure 1 retains the complete design at every physical site.
d=read('Fig1');fig,axs=plt.subplots(2,1,figsize=(10,5.8),sharex=True)
for col,label,c in [('designed_codons','NNK design','.6'),('PRE_HC','PRE high confidence','#D55E00'),('T0_HC_FINAL','T0 high confidence','#0072B2')]:axs[0].plot(range(len(d)),d[col],'.-',label=label,color=c)
for col,label,c in [('AA_designed','Designed amino acids','.6'),('AA_T0_HC','Represented at T0','#0072B2'),('AA_T0_multi_codon','Retaining multiple codons','#009E73')]:axs[1].plot(range(len(d)),d[col],'.-',label=label,color=c)
for ax in axs:ax.legend(ncol=3,frameon=False,fontsize=7);ax.set_ylabel('Design or support count')
axs[1].set_xticks(range(len(d)),d.canonical_aa_position,rotation=90,fontsize=7);axs[1].set_xlabel('Physical Fks1 position');save(fig,'Fig1')
# Figure 2 explicitly keeps ranking and error metrics together.
d=read('Fig2');fig=plt.figure(figsize=(10,7.2));gs=fig.add_gridspec(2,2,height_ratios=[1,1.2]);ax=fig.add_subplot(gs[0,:])
for j,m in enumerate(['Evo2','GEMME','ESM']):
 rows=d[(d.analysis=='zero-shot')&(d.metric=='median_signed_position_rho')&(d.contrast_or_model==m)].set_index('endpoint').loc[eps]
 x=np.arange(4)+(j-1)*.17
 ax.vlines(x,rows.conditional_low,rows.conditional_high,color=colors[j]);ax.plot(x,rows.point,'o',color=colors[j],label='ESM-1v' if m=='ESM' else m)
ax.axhline(0,color='.6',lw=.8);ax.set_xticks(range(4),epnames);ax.set_ylabel('Median signed within-position ρ');ax.legend(frameon=False,ncol=3);ax.set_title('A  Matched fixed scores',loc='left')
rows=d[(d.analysis=='supervised-LOPO')&(d.metric=='median_paired_position_rho_difference')].set_index('endpoint').loc[eps]
ax=fig.add_subplot(gs[1,0]);intervals(ax,rows,'point','conditional_low','conditional_high',epnames);ax.set_title('B  Added Evo2 rank',loc='left');ax.set_xlabel('Median paired OOF Δρ')
ax=fig.add_subplot(gs[1,1])
for j,m in enumerate(models):
 rows=d[(d.metric=='pooled_Q2')&(d.contrast_or_model==m)].set_index('endpoint').loc[eps]
 ax.plot(rows.point,np.arange(4)+(j-.5)*.09,'s',color=['.4','#0072B2'][j],label=['GEMME + ESM','Plus Evo2'][j],ms=4)
ax.set_yticks(range(4),epnames);ax.invert_yaxis();ax.axvline(0,color='.6',lw=.8);ax.legend(frameon=False,fontsize=7);ax.set_xlabel('Pooled predictive Q²');ax.set_title('C  Same saved predictions',loc='left');save(fig,'Fig2')
d=read('Fig3');d=d[d.analysis=='C1_primary_REML'];fig,axs=plt.subplots(1,2,figsize=(9,3.4))
for ax,typ,labels in [(axs[0],'condition_specific_slope',['Control','Anidulafungin','Caspofungin','Micafungin']),(axs[1],'planned_contrast',['AF − control','CF − control','MF − control','Mean drug − control'])]:
 rows=d[d.result_type==typ];intervals(ax,rows,'estimate','CI95_low','CI95_high',labels,color='.2');ax.set_xlabel('Standardized effect and model-based Wald interval')
axs[0].set_title('A  Condition-specific slopes',loc='left');axs[1].set_title('B  Planned contrasts',loc='left');save(fig,'Fig3')
def landscape(n,score_rows):
 d=read(n);fig,axs=plt.subplots(len(score_rows),3,figsize=(10,3.0*len(score_rows)),squeeze=False,sharex=True,sharey=True)
 for r,score in enumerate(score_rows):
  for j,ep in enumerate(eps[1:]):
   ax=axs[r,j];kwargs={'s':9,'alpha':.8}
   if score:kwargs.update(c=d[score],cmap='viridis',vmin=d[score].min(),vmax=d[score].max())
   else:kwargs['color']='.55'
   scatter=ax.scatter(d.fitness_none,d[ep],**kwargs);ax.axvline(0,color='.7',lw=.5);ax.axhline(0,color='.7',lw=.5)
   if r==0:ax.set_title(['Anidulafungin','Caspofungin','Micafungin'][j])
   if score:fig.colorbar(scatter,ax=ax,shrink=.7,label=score)
   ax.set_xlabel('Drug-free fitness');ax.set_ylabel('Drug − control')
 save(fig,n)
landscape('Fig4',['evo2_primary_AA_score']);landscape('S2a',[None,'gemme_primary']);landscape('S2b',['esm1v_ensemble']);landscape('S2c',[None,'evo2_primary_AA_score']);landscape('S2d',['gemme_primary','esm1v_ensemble'])
d=read('Fig5');fig=plt.figure(figsize=(10,7));gs=fig.add_gridspec(2,2);ax=fig.add_subplot(gs[0,:]);modelorder=['M0_G0','M1_G0_plus_G1','M2_G0_plus_G2','M3_G0_plus_G1_plus_G2'];panels=['OBSERVED278','COMMON124_COMPLETE','THEORETICAL490','THEORETICAL_HISTORICAL27']
for panel in panels:
 rows=d[(d.source_module=='H2')&(d.component=='total')&(d.panel==panel)].set_index('model').loc[modelorder];ax.plot(range(4),rows.Q2_OOF,'o-',label=panel,ms=3)
ax.set_xticks(range(4),['M0 Geometry','M1 + usage','M2 + tRNA','M3 + both']);ax.set_ylabel('Held-out Q² for Evo2 score');ax.legend(frameon=False,fontsize=7);ax.set_title('A  Within-amino-acid score attribution',loc='left')
ax=fig.add_subplot(gs[1,0]);rows=d[(d.source_module=='H3')&(d.component=='total')&(d.contrast=='M1-M0')].set_index('panel').loc[panels];intervals(ax,rows,'point_delta_Q2','sensitivity_low_025','sensitivity_high_975',['Observed 278','Same 124 complete','Theoretical 490','Historical 27 sites']);ax.set_xlabel('Usage increment ΔQ²');ax.set_title('B  Dirichlet sensitivity',loc='left')
ax=fig.add_subplot(gs[1,1]);rows=d[d.test_id.notna()];intervals(ax,rows,'observed_partial_spearman_hotspot_adjusted','bootstrap_95CI_low','bootstrap_95CI_high',['Glucan / control','Caspofungin / CF D−C'],color='.25');ax.set_xlabel('Hotspot-adjusted signed partial ρ');ax.set_title('C  Structural context',loc='left');save(fig,'Fig5')
d=read('S4');fig,axs=plt.subplots(1,2,figsize=(10,8),gridspec_kw={'width_ratios':[2,1]});scenarios=['RAW435','AUTHOR_REPLACEMENT435','RAW435_NO_HS2','AUTHOR_REPLACEMENT435_NO_HS2'];rows=[];labels=[]
for ep,lab in zip(eps,epnames):
 for sc in scenarios:
  row=d[(d.endpoint==ep)&(d.scenario==sc)&(d.metric=='median_paired_position_rho_difference')].iloc[0];rows.append(row);labels.append(lab+' / '+sc)
rows=pd.DataFrame(rows);intervals(axs[0],rows,'point','conditional_low','conditional_high',labels)
q=[d[(d.endpoint==row.endpoint)&(d.scenario==row.scenario)&(d.metric=='delta_pooled_Q2')].point.iloc[0] for row in rows.itertuples()]
axs[1].plot(q,range(len(q)),'.',color='.25');axs[1].set_yticks(range(len(q)),[]);axs[1].invert_yaxis();axs[1].axvline(0,color='.6',lw=.8);axs[0].set_xlabel('Δρ');axs[1].set_xlabel('ΔQ²');save(fig,'S4')
# Remaining unchanged supplementary compositions are provided by original code.
legacy=a.out.parent/'rebuilt_figures'
if not all((legacy/(n+'.pdf')).is_file() for n in ['FigS1_paired_zero_shot','FigS3_position_influence','FigS5_codon_representation','FigS6_historical_robustness']):
 import subprocess,sys
 subprocess.run([sys.executable,str(ROOT/'scripts/rebuild_figures.py'),'--data-root',str(a.data_root),'--out',str(legacy)],check=True)
for n,old in [('S1','FigS1_paired_zero_shot'),('S3','FigS3_position_influence'),('S5','FigS5_codon_representation'),('S6','FigS6_historical_robustness')]:
 if not (legacy/(old+'.pdf')).is_file():continue
 for ext in ['pdf','svg','png']:shutil.copy2(legacy/(old+'.'+ext),a.out/(n+'.'+ext))
 shutil.copy2(src/(n+'_source.csv'),a.out/(n+'_source.csv'));index.append({'figure':n,'source':n+'_source.csv','method':'original plotting module; name mapping only'})
(a.out/'CBC_FIGURE_REBUILD_INDEX.json').write_text(json.dumps(index,indent=2)+'\n');print('Recreated CBC figures:',len(index))
