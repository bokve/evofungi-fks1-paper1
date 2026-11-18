from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import argparse
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--package-root', type=Path, default=ROOT)
args=parser.parse_args()
OUT=args.package_root
(OUT / 'figures').mkdir(exist_ok=True)
T=OUT/'tables'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
blue='#176A93';orange='#CC7832';green='#32866B';purple='#86579E';gray='#64717A'
colors=[blue,orange,green,purple,gray]
eps=['fitness_none','drug_specific_anidulafungin','drug_specific_caspofungin','drug_specific_micafungin']
labels=['Drug free','AF D−C','CF D−C','MF D−C']
def save(fig,name):
    d=OUT/'figures'
    for ext in ['pdf','svg','png']:
        fig.savefig(d/f'{name}.{ext}',dpi=350,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
def panel(ax,label,title):
    ax.set_title(title,loc='left',pad=12)
    ax.text(-.14,1.06,label,transform=ax.transAxes,fontweight='bold',fontsize=12)
def forest(ax,point,lo,hi,y,color=blue):
    ax.errorbar(point,y,xerr=[np.array(point)-np.array(lo),np.array(hi)-np.array(point)],fmt='o',color=color,ms=4,elinewidth=1.25,capsize=2)
    ax.axvline(0,color='#A4ABB0',lw=.8,zorder=0)

# Figure 1: design coverage and distinct phenotype units.
d=pd.read_csv(T/'DESIGN_FIGURE_SOURCE.csv').sort_values('canonical_aa_position')
x=np.arange(len(d));site=d.canonical_aa_position.astype(int).astype(str)
fig,(a,b)=plt.subplots(2,1,figsize=(7.2,6.3),sharex=True,gridspec_kw={'hspace':.65})
for ax in [a,b]:
    for left,right,c in [(-.5,8.5,blue),(8.5,19.5,green),(19.5,27.5,orange)]:ax.axvspan(left,right,facecolor=c,alpha=.055,zorder=-2)
    ax.set_xlim(-.7,27.7);ax.grid(axis='y',color='#E7EBED',lw=.6,zorder=-1)
a.bar(x-.24,d.designed_codons,.24,label='Designed',color='#CED5DA')
a.bar(x,d.PRE_HC,.24,label='PRE high confidence',color=blue)
a.bar(x+.24,d.T0_HC_FINAL,.24,label='T0 high confidence',color=orange)
a.set_ylabel('Exact DNA genotypes');a.legend(loc='upper center',bbox_to_anchor=(.5,-.08),ncol=3,frameon=False)
panel(a,'A','Exact-codon support at each targeted position')
b.bar(x-.15,d.AA_designed,.3,color='#CED5DA',label='Designed')
b.bar(x+.15,d.AA_T0_HC,.3,color=green,label='Measured primary panel')
b.set_ylabel('Amino acid substitutions');b.set_xticks(x,site,rotation=90);b.set_xlabel('Fks1 residue position');b.legend(frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,-.46))
panel(b,'B','Measured amino acid panel within the targeted design')
for start,end,text in [(0,8,'HS1'),(9,19,'HS3'),(20,27,'HS2')]:a.text((start+end)/2,1.32,text,transform=a.get_xaxis_transform(),ha='center',fontsize=9)
fig.subplots_adjust(top=.84,bottom=.23,left=.12,right=.98)
save(fig,'Figure_1')

# Figure 2: matched signed associations; chronology is explicit in panel titles.
z=pd.read_csv(T/'ZERO_SHOT_FIGURE_SOURCE.csv')
fig,axes=plt.subplots(1,2,figsize=(7.4,4.3),sharex=True,sharey=True)
for ax,models,prefix,title in [(axes[0],['Evo2','GEMME','ESM-1v'],'A','Original fixed comparisons'),(axes[1],['ESMC600M','SaProt650M'],'B','Additional comparisons')]:
    for j,m in enumerate(models):
        q=z[z.model==m].set_index('endpoint').loc[eps]
        offset=(j-(len(models)-1)/2)*.18
        cc={'Evo2':blue,'GEMME':orange,'ESM-1v':green,'ESMC600M':purple,'SaProt650M':gray}[m]
        forest(ax,q.point,q.low,q.high,np.arange(4)+offset,cc)
        ax.plot([],[],color=cc,marker='o',ms=4,label={'Evo2':'Evo 2','ESMC600M':'ESM-C','SaProt650M':'SaProt'}.get(m,m))
    ax.set_xlim(-.43,.44);ax.set_yticks(range(4),labels);ax.set_xlabel('Median within-position Spearman ρ')
    ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.24),ncol=2 if len(models)==2 else 3)
    ax.grid(axis='x',color='#EBEFF1',lw=.6);panel(ax,prefix,title)
axes[0].invert_yaxis();fig.subplots_adjust(left=.13,right=.98,bottom=.29,top=.84,wspace=.25)
save(fig,'Figure_2')

# Figure 3: ranking, prediction errors, and their distinct uncertainty procedures.
q=pd.read_csv(T/'INCREMENTAL_FIGURE_SOURCE.csv')
fig,axes=plt.subplots(1,3,figsize=(8.4,4.1),gridspec_kw={'wspace':.6})
y=np.arange(4)
forest(axes[0],q.rank_delta,q.rank_low,q.rank_high,y)
axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[0].set_xlabel('Paired rank Δρ')
panel(axes[0],'A','Within-position ranking')
for i,r in q.iterrows():
    axes[1].plot([r.q_baseline,r.q_augmented],[i,i],color='#AAB4BB',lw=1.4)
axes[1].scatter(q.q_baseline,y,color=orange,label='GEMME + ESM-1v',s=25)
axes[1].scatter(q.q_augmented,y,color=blue,label='Plus Evo 2',s=25)
axes[1].set_xlabel('Pooled predictive Q²');axes[1].set_yticks(y,[]);axes[1].invert_yaxis();axes[1].axvline(0,color='#A4ABB0',lw=.8);axes[1].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.25),fontsize=7)
panel(axes[1],'B','Numerical prediction')
forest(axes[2],q.q_delta,q.q_refit_low,q.q_refit_high,y)
axes[2].set_yticks(y,[]);axes[2].invert_yaxis();axes[2].set_xlabel('Pooled ΔQ², training refits')
panel(axes[2],'C','Conditional uncertainty')
for ax in axes:ax.grid(axis='x',color='#EBEFF1',lw=.6)
fig.subplots_adjust(left=.10,right=.99,bottom=.30,top=.84)
save(fig,'Figure_3')

# Figure 4: effect estimates remain separate from conditional calibration.
h=pd.read_csv(T/'HIERARCHICAL_FIGURE_SOURCE.csv');u=pd.read_csv(T/'U9_NULL_DRAWS.csv');ut=pd.read_csv(T/'U9_CALIBRATION.csv')
fig=plt.figure(figsize=(7.4,6.5));gs=fig.add_gridspec(2,2,height_ratios=[1,1.1],hspace=.65,wspace=.6)
a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,:])
sl=h[(h.analysis=='C1_primary_REML')&(h.result_type=='condition_specific_slope')]
forest(a,sl.estimate,sl.CI95_low,sl.CI95_high,range(4));a.set_yticks(range(4),['Drug free','AF','CF','MF']);a.invert_yaxis();a.set_xlabel('Standardized slope');panel(a,'A','Condition-specific associations')
co=h[(h.analysis=='C1_primary_REML')&(h.result_type=='planned_contrast')]
forest(b,co.estimate,co.CI95_low,co.CI95_high,range(4),orange);b.set_yticks(range(4),['AF − C','CF − C','MF − C','Mean drug − C']);b.invert_yaxis();b.set_xlabel('Standardized slope difference');panel(b,'B','Planned contrasts')
for arm,col,lab in [('G',blue,'Gaussian'),('T5E',orange,'t5 errors'),('T5ALL',green,'t5 errors and effects')]:
    c.hist(u[u.arm==arm].LR,bins=np.linspace(0,max(ut.observed_LR.max()+2,u.LR.max()+1),45),histtype='step',color=col,lw=1.3,density=True,label=lab)
obs=ut.observed_LR.iloc[0];c.axvline(obs,color='#202A31',ls='--',lw=1.3,label=f'Observed LR = {obs:.3f}')
c.set_xlabel('Likelihood-ratio statistic after refitting');c.set_ylabel('Simulated density');c.legend(frameon=True,facecolor='white',framealpha=1,edgecolor='white',fontsize=8,ncol=2,loc='upper right');panel(c,'C','Conditional null calibration: 999 valid draws per scenario')
fig.subplots_adjust(left=.13,right=.98,top=.94,bottom=.09)
save(fig,'Figure_4')

# Figure 5: attribution predicts computational scores, not biological fitness.
v=pd.read_csv(T/'CODON_STRUCTURE_FIGURE_SOURCE.csv');model_order=['M0_G0','M1_G0_plus_G1','M2_G0_plus_G2','M3_G0_plus_G1_plus_G2']
fig,(a,b)=plt.subplots(1,2,figsize=(7.6,3.9),gridspec_kw={'wspace':.58})
for j,(p,label,col) in enumerate([('OBSERVED278','Observed 278',blue),('THEORETICAL490','Theoretical 490',orange)]):
    vv=v[(v.panel==p)&(v.component=='total')&v.model.isin(model_order)].set_index('model').loc[model_order]
    a.bar(np.arange(4)+(j-.5)*.33,vv.Q2_OOF,.33,color=col,label=label)
a.set_xticks(range(4),['Geometry','+ usage','+ tRNA','+ both'],rotation=30,ha='right');a.set_ylabel('Held-out Q² for Evo 2 score');a.set_ylim(0,.69);a.legend(frameon=False,fontsize=8,loc='upper left');panel(a,'A','Predicting codon-score variation')
vv=v[(v.component=='total')&(v.contrast=='M1-M0')].set_index('panel').loc[['OBSERVED278','COMMON124_COMPLETE','THEORETICAL490','THEORETICAL_HISTORICAL27']]
forest(b,vv.point_delta_Q2,vv.sensitivity_low_025,vv.sensitivity_high_975,range(4),green)
b.set_yticks(range(4),['Observed 278','Same 124 complete','Theoretical 490','Historical 27 sites']);b.invert_yaxis();b.set_xlabel('Usage increment ΔQ²');panel(b,'B','Support-preserving sensitivity')
fig.subplots_adjust(left=.1,right=.98,bottom=.24,top=.86)
save(fig,'Figure_5')
