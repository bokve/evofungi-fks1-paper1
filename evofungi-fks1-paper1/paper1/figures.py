"""Publication exports generated exclusively from saved source tables."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from .panels import DRUGS,ENDPOINTS
from .stages import atomic_json,digest

COLORS={'Evo2':'#0072B2','GEMME':'#D55E00','ESM':'#009E73','GEMME_ESM':'#D55E00','GEMME_ESM_Evo2':'#0072B2'}
LABELS={'fitness_none':'Drug-free','drug_specific_anidulafungin':'Anidulafungin D−C','drug_specific_caspofungin':'Caspofungin D−C','drug_specific_micafungin':'Micafungin D−C'}


def run(inputs,attempt,config):
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
    index=[];captions=[]
    def save(fig,name,source,caption):
        source_path=attempt/(name+'_source.csv');source.to_csv(source_path,index=False)
        fig.tight_layout()
        for suffix in ['pdf','svg','png']:fig.savefig(attempt/(name+'.'+suffix),dpi=300,bbox_inches='tight')
        plt.close(fig)
        index.append({'figure':name,'source_CSV':source_path.name,'source_sha256':digest(source_path),
                      'exports':[name+'.'+s for s in ['pdf','svg','png']],'code':'paper1/figures.py',
                      'code_sha256':digest(Path(__file__)),'caption':caption,'visual_inspection':'PENDING_ACTUAL_RENDER_REVIEW'})
        captions.append('## '+name+'\n\n'+caption+'\n')
    positions=pd.read_csv(inputs['hs2-core/HS2_CORE_SOURCE_TABLE.csv'])
    fig,axes=plt.subplots(2,1,figsize=(10,6.8),sharex=True);ax=axes[0];x=np.arange(len(positions))
    fig.suptitle(f"NNK design: {positions.designed_codons.sum()} codons / {positions.AA_designed.sum()} AA; T0-HC: {positions.T0_HC_FINAL.sum()} codons / {positions.AA_T0_HC.sum()} AA",fontsize=11)
    ax.bar(x,positions.designed_codons,color='#DDDDDD',label='NNK design')
    ax.bar(x,positions.PRE_HC,width=.65,color='#56B4E9',label='PRE high confidence')
    ax.bar(x,positions.T0_HC_FINAL,width=.35,color='#0072B2',label='T0-HC = FINAL')
    ax.set_xticks(x);ax.set_xticklabels(positions.canonical_aa_position,rotation=90)
    ax.set(ylabel='Exact designed codons')
    ax.legend(ncol=3,frameon=False,loc='upper center',bbox_to_anchor=(.5,1.24))
    ax=axes[1]
    ax.bar(x,positions.AA_designed,color='#DDDDDD',label='Designed AA groups')
    ax.bar(x,positions.AA_T0_HC,width=.65,color='#56B4E9',label='T0-HC represented AA')
    ax.bar(x,positions.AA_T0_multi_codon,width=.35,color='#009E73',label='T0-HC AA with ≥2 codons')
    ax.set_xticks(x);ax.set_xticklabels(positions.canonical_aa_position,rotation=90)
    ax.set(xlabel='Canonical Fks1 position (hotspot blocks)',ylabel='AA substitution groups')
    ax.legend(ncol=3,frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,1.24))
    for axis in axes:
        for hs,group in positions.groupby('hotspot',sort=False):
            ids=group.index.to_numpy();axis.text(float(ids.mean()),1.015,hs,transform=axis.get_xaxis_transform(),ha='center',fontsize=8)
            if ids[0]>0:axis.axvline(ids[0]-.5,color='.75',lw=.6,ls=':')
    save(fig,'Fig1_ascertainment',positions,'Exact-genotype and AA-group representation across the targeted design. Top: exhaustive NNK missense design, source PRE high-confidence and T0 high-confidence support; FINAL589 equals T0-HC and is not a later survival stage. Bottom: designed AA groups, represented AA groups and groups retaining multiple exact codons. L1360 remains represented in the theoretical design but has no observed informative multi-codon AA group. Canonical positions derive from sequence. Support loss is not lethality. No uncertainty bars: exhaustive design descriptions; historical B3 reused.')
    summary=pd.read_csv(inputs['compare/COMPARISON_SUMMARIES.csv'])
    primary=summary[(summary.scenario=='RAW435') & (summary.analysis=='zero-shot') & (summary.metric=='median_signed_position_rho')]
    fig,axes=plt.subplots(1,4,figsize=(11,3.2),sharey=True)
    for ax,endpoint in zip(axes,ENDPOINTS):
        d=primary[primary.endpoint==endpoint]
        for i,model in enumerate(['Evo2','GEMME','ESM']):
            row=d[d.contrast_or_model==model].iloc[0]
            ax.vlines(i,row.conditional_low,row.conditional_high,color=COLORS[model],lw=2)
            ax.scatter(i,row.point,color=COLORS[model],s=32,zorder=3)
            ax.text(i,-1.06,str(int(row.n_positions)),ha='center',fontsize=8)
        ax.axhline(0,color='.7',lw=.7);ax.set_xticks(range(3),['Evo2','GEMME','ESM'],rotation=35)
        ax.set_title(LABELS[endpoint]);ax.set_ylim(-1.13,1.03)
    axes[0].set_ylabel('Median signed within-position Spearman ρ')
    save(fig,'Fig2_matched_benchmark',primary,'RAW435 matched zero-shot benchmark. Equal-position median signed Spearman correlations, using positions estimable for all three predictors on the same finite AA panel (counts beneath points). Lines are pointwise 95% conditional hotspot-stratified position-resampling intervals (10,000 draws); they condition on the observed panel and exclude phenotype measurement uncertainty. They are not formal superiority tests or new-landscape generalization CIs. Drug−control endpoints are secondary.')
    contrasts=summary[(summary.scenario=='RAW435') & (summary.analysis=='zero-shot') & (summary.metric=='median_paired_position_rho_difference')]
    fig,ax=plt.subplots(figsize=(8,4));labels=[]
    for i,row in enumerate(contrasts.itertuples()):
        ax.hlines(i,row.conditional_low,row.conditional_high,color='#777777');ax.scatter(row.point,i,color='#0072B2',s=22)
        labels.append(LABELS[row.endpoint]+' | '+row.contrast_or_model)
    ax.axvline(0,color='.5',lw=.7);ax.set_yticks(range(len(labels)),labels,fontsize=8);ax.invert_yaxis();ax.set_xlabel('Median paired position-level difference in signed ρ')
    save(fig,'FigS1_paired_zero_shot',contrasts,'All planned paired zero-shot contrasts. Median of position-level rho differences, not difference of medians. Conditional resampling intervals as in Figure 2; all planned contrasts/endpoints shown without outcome selection.')
    new2=pd.read_csv(inputs['history/NEW2_summary.csv'])
    old=new2[(new2.analysis=='C1_primary_REML') & (new2.result_type=='condition_specific_slope')]
    inc=summary[(summary.scenario=='RAW435') & (summary.analysis=='supervised-LOPO') & (summary.metric=='median_paired_position_rho_difference')]
    fig,axes=plt.subplots(1,2,figsize=(10,3.4))
    for i,row in enumerate(old.itertuples()):
        axes[0].hlines(i,row.CI95_low,row.CI95_high,color='#777777');axes[0].scatter(row.estimate,i,color='#0072B2')
    axes[0].set_yticks(range(len(old)),[s.replace('evo2_slope_','') for s in old.term]);axes[0].set_xlabel('Standardized within-position slope');axes[0].set_title('Historical NEW-2: model-based Wald intervals')
    for i,row in enumerate(inc.itertuples()):
        axes[1].hlines(i,row.conditional_low,row.conditional_high,color='#777777');axes[1].scatter(row.point,i,color='#009E73')
    axes[1].set_yticks(range(len(inc)),[LABELS[s] for s in inc.endpoint]);axes[1].set_xlabel('Median paired OOF ρ increment');axes[1].set_title('New fixed LOPO: conditional intervals')
    for ax in axes:ax.axvline(0,color='.5',lw=.7);ax.invert_yaxis()
    source=pd.concat([old.assign(source_module='historical_NEW2'),inc.assign(source_module='new_ESM_LOPO')],ignore_index=True)
    save(fig,'Fig3_slopes_and_increment',source,'Left: reused NEW-2 condition-specific standardized within-position slopes, with original model-based 95% Wald intervals; repeated-measure hierarchical model, not refitted. Right: new supervised position-held-out OLS comparison, GEMME+ESM+Evo2 minus GEMME+ESM, summarized as median paired OOF rho difference. Right-hand intervals resample fixed OOF position records; overlapping LOPO training folds and upstream phenotype uncertainty prevent interpretation as validated generalization CIs. Panels answer different questions and have different units.')
    for panel,fig_name in [('RAW435','Fig4_NEW4_landscape'),('AUTHOR_REPLACEMENT435','FigS2_author_landscape')]:
        frame=pd.read_csv(inputs['panels/'+panel+'.csv'])
        esm=pd.read_csv(inputs['esm-score/ESM_RAW435.csv'])
        frame=frame.merge(esm[['aa_variant_id','esm1v_ensemble']],on='aa_variant_id',validate='one_to_one')
        fig,axes=plt.subplots(3,4,figsize=(12,8),sharex=True,sharey='row')
        overlays=[None,'evo2_primary_AA_score','gemme_primary','esm1v_ensemble']
        for j,col in enumerate(overlays):
            norm=Normalize(frame[col].min(),frame[col].max()) if col else None
            for i,drug in enumerate(DRUGS):
                ax=axes[i,j]
                scat=ax.scatter(frame.fitness_none,frame['drug_specific_'+drug],c=frame[col] if col else '#777777',cmap='viridis' if col else None,norm=norm,s=12,alpha=.8,edgecolors='none')
                ax.axhline(0,color='.7',lw=.6);ax.axvline(0,color='.7',lw=.6)
                if j==0:ax.set_ylabel(drug+' D−C')
                if i==2:ax.set_xlabel('Drug-free fitness')
                if i==0:ax.set_title(['Phenotypes','Evo2','GEMME','ESM'][j])
                if col:fig.colorbar(scat,ax=ax,shrink=.7,pad=.02)
        save(fig,fig_name,frame,'Original NEW-4 functional-fitness × drug-selective-effect landscape, '+panel+'. One point per AA substitution. Three drug−control outcomes versus drug-free fitness, with plain phenotype and fixed model overlays. All values shown; each score uses a consistent full-range scale across drugs. Zero axes do not define clinical categories. No classifier, optimized threshold, clinical resistance or causal trade-off inference. Author replacements are a separately labeled source sensitivity.')
    h2=pd.read_csv(inputs['history/H2_points.csv']);h3=pd.read_csv(inputs['history/H3_deltas.csv']);structure=pd.read_csv(inputs['history/NEW1_primary.csv'])
    fig,axes=plt.subplots(1,3,figsize=(12,3.6))
    components=['total','codon_local','downstream50'];colors=['#0072B2','#009E73','#D55E00']
    for comp,color in zip(components,colors):
        d=h2[(h2.panel=='THEORETICAL490') & (h2.component==comp)].sort_values('model')
        axes[0].plot(range(len(d)),d.Q2_OOF,'o-',label=comp,color=color)
    axes[0].set_xticks(range(4),['M0','M1','M2','M3']);axes[0].set_ylabel('Within-AA model-score OOF Q²');axes[0].set_title('Historical H2: theoretical490');axes[0].legend(fontsize=7,frameon=False)
    d=h3[(h3.component=='total') & (h3.contrast=='M1-M0')]
    for i,row in enumerate(d.itertuples()):
        axes[1].hlines(i,row.sensitivity_low_025,row.sensitivity_high_975,color='.5');axes[1].scatter(row.point_delta_Q2,i,color='#009E73')
    axes[1].set_yticks(range(len(d)),d.panel,fontsize=7);axes[1].set_xlabel('ΔQ² (M1−M0)');axes[1].set_title('Historical H3: sensitivity intervals');axes[1].axvline(0,color='.5',lw=.7)
    for i,row in enumerate(structure.itertuples()):
        axes[2].hlines(i,row.bootstrap_95CI_low,row.bootstrap_95CI_high,color='.5');axes[2].scatter(row.observed_partial_spearman_hotspot_adjusted,i,color='#0072B2')
    axes[2].set_yticks(range(len(structure)),structure.test_id);axes[2].set_xlabel('Hotspot-adjusted signed partial ρ');axes[2].set_title('Historical NEW-1: bootstrap intervals');axes[2].axvline(0,color='.5',lw=.7)
    source=pd.concat([h2.assign(source_module='H2'),h3.assign(source_module='H3'),structure.assign(source_module='NEW1')],ignore_index=True)
    save(fig,'Fig5_codon_and_structure',source,'Reused historical evidence, not newly refitted. H2 predicts within-AA Evo2 score variation rather than experimental fitness; G0 edit geometry is technical attribution. H3 bars are support-preserving Dirichlet sensitivity intervals, not sampling CIs, posteriors or p-values. NEW-1 shows the original signed structural-context tests, preserving uncertain/negative findings and original bootstrap intervals. Caspofungin structural context does not establish equivalent binding geometry for the other drugs. No new structural features or mutant simulations.')
    influence=pd.read_csv(inputs['compare/POSITION_DELETION_INFLUENCE.csv'])
    selected=influence[(influence.scenario=='RAW435') & (influence.analysis=='supervised-LOPO') & (influence.metric=='median_paired_position_rho_difference')]
    fig,axes=plt.subplots(1,4,figsize=(11,3.2),sharey=True)
    for ax,ep in zip(axes,ENDPOINTS):
        d=selected[selected.endpoint==ep].sort_values('deleted_position');ax.scatter(range(len(d)),d.point,color='#0072B2',s=15)
        ax.axhline(0,color='.5',lw=.7);ax.set_title(LABELS[ep]);ax.set_xticks(range(len(d)));ax.set_xticklabels(d.deleted_position,rotation=90,fontsize=6)
    axes[0].set_ylabel('LOPO ρ increment after deletion');save(fig,'FigS3_position_influence',selected,'Exhaustive single-position deletion stability. Each point removes the named position from training and testing and rebuilds both LOPO models. All positions are shown; these are stability summaries, not independent tests or confidence intervals. Full RAW435 remains primary.')
    sensitivity=summary[(summary.analysis=='supervised-LOPO') & (summary.metric=='median_paired_position_rho_difference')]
    fig,ax=plt.subplots(figsize=(9,4));labels=[]
    for i,row in enumerate(sensitivity.itertuples()):
        ax.hlines(i,row.conditional_low,row.conditional_high,color='.6');ax.scatter(row.point,i,color='#009E73',s=20)
        labels.append(row.scenario+' | '+LABELS[row.endpoint])
    ax.set_yticks(range(len(labels)),labels,fontsize=7);ax.invert_yaxis();ax.axvline(0,color='.5',lw=.7);ax.set_xlabel('Median paired OOF ρ increment')
    save(fig,'FigS4_fixed_sensitivities',sensitivity,'All prespecified new LOPO sensitivity panels/endpoints: primary raw, author replacement, and removal of all HS2 positions from both training and testing. Fixed features/models and shared AA masks within comparisons. Conditional fixed-OOF intervals as in Figure 3; no sensitivity was selected as a new primary.')
    multiplicity=pd.read_csv(inputs['hs2-core/HS2_AA_MULTIPLICITY.csv'])
    multi=multiplicity.loc[multiplicity.n_theoretical_codons>=2].copy()
    multi['T0_state']=np.where(multi.n_T0_HC_codons==0,'No HC codon',np.where(multi.n_T0_HC_codons==1,'One HC codon','Multiple HC codons'))
    counts=pd.crosstab(multi.hotspot,multi.T0_state).reindex(columns=['No HC codon','One HC codon','Multiple HC codons'],fill_value=0)
    fig,ax=plt.subplots(figsize=(6,3.4));bottom=np.zeros(len(counts))
    for label,color in zip(counts.columns,['#999999','#E69F00','#0072B2']):
        ax.bar(counts.index,counts[label],bottom=bottom,label=label,color=color);bottom+=counts[label].to_numpy()
    ax.set_ylabel('Theoretical multi-codon AA groups');ax.legend(frameon=False,fontsize=8)
    save(fig,'FigS5_codon_representation',counts.reset_index(),'Categorical representation of the theoretical multi-codon AA universe, separated by hotspot. Designed codons and observed support are distinct units. These reused B3 descriptions do not quantify absolute survival or continuous depletion; the H208 feature gate is separate.')
    m5=pd.read_csv(inputs['history/M5_robustness.csv']);rel=pd.read_csv(inputs['history/method3_reliability.csv'])
    rel=rel[rel.scope_type=='global_equal_AA']
    fig,axes=plt.subplots(1,2,figsize=(10,3.5))
    for i,row in enumerate(m5.itertuples()):
        axes[0].hlines(i,row.bootstrap_CI_low,row.bootstrap_CI_high,color='.5');axes[0].scatter(row.median_delta_combo_minus_GEMME,i,color='#0072B2')
    axes[0].set_yticks(range(len(m5)),[r.analysis+' (n='+str(r.n_positions)+')' for r in m5.itertuples()],fontsize=7)
    axes[0].set_xlabel('Historical median incremental ρ');axes[0].set_title('M5 robustness: original masks')
    for i,(_,row) in enumerate(rel.iterrows()):
        axes[1].hlines(i,row['bootstrap_ci_2.5'],row['bootstrap_ci_97.5'],color='.5');axes[1].scatter(row.mean_group_tau_b,i,color='#009E73')
    axes[1].set_yticks(range(len(rel)),rel.condition);axes[1].set_xlabel('Historical mean within-AA replicate τ-b');axes[1].set_title('Method 3: measurement reliability')
    for ax in axes:ax.axvline(0,color='.5',lw=.7);ax.invert_yaxis()
    save(fig,'FigS6_historical_robustness',pd.concat([m5.assign(source_module='M5'),rel.assign(source_module='Method3')],ignore_index=True),'Reused M5 and Method 3 summaries with their original bootstrap intervals. The M5 combined/independent/epistatic estimable-position denominators differ and are printed explicitly; a favorable alternative component is not substituted as the primary. Method 3 describes the existing within-AA measurement-resolution check, not an independent codon-fitness ground truth or new clinical validation. No historical refit or uncertainty recalculation.')
    atomic_json(attempt/'FIGURE_INDEX.json',index)
    (attempt/'CAPTIONS.md').write_text('# Figure captions\n\n'+'\n'.join(captions))
    return [p.name for p in attempt.iterdir() if p.name!='started.json']
