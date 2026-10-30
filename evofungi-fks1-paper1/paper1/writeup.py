"""English drafts and claim ledger generated from actual frozen/new outputs."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .panels import ENDPOINTS,DRUGS
from .stages import atomic_json,digest


def markdown_table(frame,columns=None):
    f=frame[columns] if columns else frame
    def value(v):
        if pd.isna(v):return 'NA'
        if isinstance(v,(float,np.floating)):return f'{v:.6g}'
        return str(v).replace('|','/').replace('\n',' ')
    return '| '+' | '.join(f.columns)+' |\n| '+' | '.join(['---']*len(f.columns))+' |\n'+'\n'.join('| '+' | '.join(value(v) for v in row)+' |' for row in f.itertuples(index=False,name=None))+'\n'


def run(inputs,attempt,config):
    summary=pd.read_csv(inputs['compare/COMPARISON_SUMMARIES.csv'])
    primary=summary[summary.scenario=='RAW435']
    zero=primary[(primary.analysis=='zero-shot')&(primary.metric=='median_signed_position_rho')]
    paired=primary[(primary.analysis=='zero-shot')&(primary.metric=='median_paired_position_rho_difference')]
    increments=primary[(primary.analysis=='supervised-LOPO')&(primary.metric=='median_paired_position_rho_difference')]
    q2=primary[(primary.analysis=='supervised-LOPO')&primary.metric.isin(['pooled_Q2','equal_position_Q2','MSE','MAE'])]
    n2=pd.read_csv(inputs['history/NEW2_summary.csv'])
    h2=pd.read_csv(inputs['history/H2_points.csv']);h3=pd.read_csv(inputs['history/H3_deltas.csv'])
    structure=pd.read_csv(inputs['history/NEW1_primary.csv']);m5=pd.read_csv(inputs['history/M5_robustness.csv'])
    reliability=pd.read_csv(inputs['history/method3_reliability.csv'])
    influence=pd.read_csv(inputs['compare/POSITION_DELETION_INFLUENCE.csv'])
    raw=pd.read_csv(inputs['panels/RAW435.csv']);ledger=pd.read_csv(inputs['panels/CORRECTION_LEDGER.csv'])
    coverage=json.loads(inputs['esm-score/ESM_COVERAGE.json'].read_text())
    acquisition=json.loads(inputs['esm-score/CHECKPOINT_ACQUISITION_MANIFEST.json'].read_text())
    protocol=json.loads(inputs['protocol'].read_text())
    def estimate(row):
        if not np.isfinite(row.point):return f'not estimable under the fixed mask ({int(row.n_positions)} eligible positions)'
        return f"{row.point:+.3f} (95% conditional position-resampling interval {row.conditional_low:+.3f} to {row.conditional_high:+.3f}; {int(row.n_positions)} positions)"
    control={m:zero[(zero.endpoint=='fitness_none')&(zero.contrast_or_model==m)].iloc[0] for m in ['Evo2','GEMME','ESM']}
    added=increments[increments.endpoint=='fitness_none'].iloc[0]
    drugtext='; '.join(d+': '+estimate(increments[increments.endpoint=='drug_specific_'+d].iloc[0]) for d in DRUGS)
    controltext='; '.join(m+': '+estimate(control[m]) for m in control)
    oldlrt=n2[n2.result_type=='primary_omnibus_LRT'].iloc[0]
    influences=influence[(influence.scenario=='RAW435')&(influence.endpoint=='fitness_none')&(influence.analysis=='supervised-LOPO')&(influence.metric=='median_paired_position_rho_difference')]
    def value_at(scenario,endpoint,metric,label='GEMME_ESM_Evo2 minus GEMME_ESM'):
        return float(summary[(summary.scenario==scenario)&(summary.endpoint==endpoint)&(summary.metric==metric)&(summary.contrast_or_model==label)].point.iloc[0])
    baseline_control=value_at('RAW435','fitness_none','pooled_Q2','GEMME_ESM')
    augmented_control=value_at('RAW435','fitness_none','pooled_Q2','GEMME_ESM_Evo2')
    negative_q2_increments=int(((primary.analysis=='supervised-LOPO')&(primary.metric=='delta_pooled_Q2')&(primary.point<0)).sum())
    replacement_control=summary[(summary.scenario=='AUTHOR_REPLACEMENT435')&(summary.endpoint=='fitness_none')&(summary.analysis=='supervised-LOPO')&(summary.metric=='median_paired_position_rho_difference')].iloc[0]
    nohs_control=summary[(summary.scenario=='RAW435_NO_HS2')&(summary.endpoint=='fitness_none')&(summary.analysis=='supervised-LOPO')&(summary.metric=='median_paired_position_rho_difference')].iloc[0]
    historical=json.loads(inputs['history/HISTORICAL_SOURCE_LEDGER.json'].read_text())
    f2=next(r for r in historical if r['role']=='history_P1_NEW3F2_BOOTSTRAP_UNCERTAINTY_v1.0')
    f2_primary=f2['record']['primary_same_AA_bootstrap'];f2_m0=f2_primary['diagnostics']['M0']
    lineage=json.loads(inputs['lineage/METHOD3_UNIT_AUDIT.json'].read_text())
    columns=['endpoint','contrast_or_model','point','conditional_low','conditional_high','n_positions','n_AA']
    support_note='The separate categorical feature/identity branch is blocked: '+json.loads(inputs['support_disposition'].read_text()).get('block_reason','see terminal disposition')+'. No score-support association is substituted from the measured-AA survivor panel.'
    if 'support/HS2_SUPPORT_SCORE_GATE.json' in inputs:
        gate=json.loads(inputs['support/HS2_SUPPORT_SCORE_GATE.json'].read_text())
        if gate['status']=='EXECUTED_VALIDATED':
            support_results=pd.read_csv(inputs['support/SUPPORT_SCORE_ASSOCIATIONS.csv'])
            support_note=f"The independent categorical feature gate passed on {gate['candidate_AA']} theoretical multi-codon AA groups ({gate['exact_codons']} exact codons), including {gate['zero_T0_support_AA']} groups with no T0-HC codon. The eligible panel and protocol were hash-locked before these associations. The Evo2 predictor here is the explicitly new theoretical-H codon median, not the primary phenotype-panel AA score. All three exploratory score associations and paired contrasts are shown below. The response is categorical designed-codon representation fraction, not continuous PRE→T0 depletion or biological replicate survival.\n\n"+markdown_table(support_results)
        else:support_note='Categorical score-support branch: '+gate['status']+'. '+gate.get('reason','See feature-universe gate.')
    claims=[];table_index=[]
    # Every numerical output table remains linked to actual row and byte hash.
    for role,path in inputs.items():
        if path.suffix!='.csv' or role=='protocol':continue
        frame=pd.read_csv(path)
        source_hash=digest(path)
        table_index.append({'role':role,'path':str(path),'sha256':source_hash,'rows':len(frame)})
        if len(frame)>1000 or ('aa_variant_id' in frame and len(frame)>50):continue
        for i,row in frame.iterrows():
            numeric={k:(None if pd.isna(v) else v) for k,v in row.items() if isinstance(v,(int,float,np.number))}
            if not numeric:continue
            claims.append({'claim_id':f'C{len(claims)+1:05d}','source_type':'reused' if role.startswith('history/') else 'new',
                'source_path':str(path),'source_sha256':source_hash,'CSV_row_including_header':i+2,
                'row_identity':json.dumps({k:str(row[k]) for k in ['scenario','endpoint','analysis','contrast_or_model','panel','component','model','term','test_id','aa_variant_id'] if k in row},sort_keys=True),
                'quantitative_record':json.dumps(numeric,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),
                'status':'SUPPORTED_AS_REPORTED_ESTIMATE','protocol_sha256':config['protocol_sha256'],
                'claim_text_or_template':'Reported table record: '+role+'; '+', '.join(f'{k}={row[k]}' for k in ['scenario','endpoint','analysis','contrast_or_model','metric','panel','component','model','term','test_id'] if k in row),
                'figure_or_table':'MANUSCRIPT/SI and TABLE_INDEX; figure sources retain the same row keys',
                'uncertainty_rule':'original definition for reused records; conditional fixed-record position resampling for new rank intervals',
                'prohibited_inference':'No clinical resistance/lethality/causal mechanism/equivalence; no superiority test from pointwise intervals.'})
    for role,pointer,text in [
        ('panels/RAW435.csv','count(rows),nunique(aa_position)',f'Primary phenotype panel contains {len(raw)} AA substitutions at {raw.aa_position.nunique()} physical positions; codons/repeated rows are separate units.'),
        ('esm-score/ESM_COVERAGE.json','/','Actual all-five ensemble coverage and 140 position jobs; primary 435 and theoretical 532 AA are separate universes.'),
        ('lineage/METHOD3_UNIT_AUDIT.json','/',f"All {lineage['rows_checked']} exported Method 3 values match original NT coefficients; {lineage['groups_with_distinct_NT_coefficients']} multi-codon AA/condition/replicate groups have distinct coefficients, not copied AA aggregates."),
        ('history/HISTORICAL_SOURCE_LEDGER.json','role=history_P1_NEW3F2_BOOTSTRAP_UNCERTAINTY_v1.0/record/primary_same_AA_bootstrap/diagnostics/M0',f"Historical F2: valid-only Q2 uncertainty with {f2_m0['Q2_valid_replicates']} valid / {f2_m0['Q2_nonestimable_replicates']} nonestimable draws and {f2_m0['LOPO_nonestimable_fit_events']} nonestimable fold-fit events for M0."),
        ('compare/POSITION_DELETION_INFLUENCE.csv','scenario=RAW435;endpoint=fitness_none;analysis=supervised-LOPO;metric=median_paired_position_rho_difference','Control incremental-rho deletion range and positive-estimate count are aggregates of all 28 saved deletion records.')]:
        path=inputs[role]
        claims.append({'claim_id':f'C{len(claims)+1:05d}','source_type':'reused' if role.startswith('history/') else 'new','claim_text_or_template':text,
            'source_path':str(path),'source_sha256':digest(path),'row_identity':pointer,'status':'SUPPORTED_AS_REPORTED_ESTIMATE',
            'figure_or_table':'MANUSCRIPT Methods/Results; SI; linked source JSON/CSV','protocol_sha256':config['protocol_sha256'],
            'uncertainty_rule':'F2 retains original valid-only estimability limits; deletion range is stability, not uncertainty interval',
            'prohibited_inference':'No independent-replication or population-generalization claim.'})
    interpretations=[('General functional/evolutionary constraint explains the observed associations','UNCERTAIN','Associations alone do not identify a causal biological component.'),
      ('Frozen Evo2 adds no information beyond GEMME and ESM','NOT_ESTABLISHED','Effect estimates/conditional intervals cannot establish equivalence or absence of all information.'),
      ('The present models predict clinical MIC or patient resistance','NOT_TESTED','No clinical prediction task or external clinical validation.'),
      ('Loss of T0 high-confidence support proves death','NOT_TESTED','Support is assay ascertainment, not absolute viability.'),
      ('NEW-3I synthetic panels were executed','CONTRADICTED','The locked dated protocol explicitly defers this scope.'),
      ('Structural proximity causally explains model behavior','UNCERTAIN','Original signed tests and intervals do not establish mechanism.')]
    for text,status,reason in interpretations:
        claims.append({'claim_id':f'C{len(claims)+1:05d}','source_type':'interpretation','claim_text':text,'status':status,
            'source_path':str(inputs['protocol']),'source_sha256':digest(inputs['protocol']),'row_identity':'protocol claim limits and linked full result tables',
            'prohibited_inference':reason,'protocol_sha256':config['protocol_sha256']})
    pd.DataFrame(claims).to_csv(attempt/'CLAIM_LEDGER.csv',index=False)
    atomic_json(attempt/'STUDY_COUNTS.json',{'AA_substitutions':len(raw),'physical_positions':int(raw.aa_position.nunique()),
      'raw_source_path':str(inputs['panels/RAW435.csv']),'raw_source_sha256':digest(inputs['panels/RAW435.csv']),
      'actual_ESM_coverage':coverage,'coverage_source_path':str(inputs['esm-score/ESM_COVERAGE.json']),'coverage_source_sha256':digest(inputs['esm-score/ESM_COVERAGE.json'])})
    pd.DataFrame(table_index).to_csv(attempt/'TABLE_INDEX.csv',index=False)
    admin='Authors, order, affiliations, corresponding author, funding, conflicts, ethics statement and CRediT roles: [AUTHOR INPUT REQUIRED]. No administrative information has been inferred.'
    manuscript=f'''# Genomic sequence likelihood across functional and echinocandin-selective Fks1 landscapes

Computational manuscript draft — 18 October 2026. Scientific results are separated into reused historical modules and the newly executed ESM extension. {admin}

## Abstract

Sequence likelihood can reflect evolutionary constraints without directly measuring a drug-response phenotype. We investigated this distinction in the published *Saccharomyces cerevisiae* Fks1 hotspot deep-mutational scan [1], retaining drug-free fitness, total drug fitness and drug-minus-control selective effects as separate quantities. The primary panel contained {len(raw)} amino-acid substitutions at {raw.aa_position.nunique()} physical positions. Frozen genomic Evo2 scores and GEMME combined scores were compared with a {coverage['checkpoints']}-checkpoint ESM-1v masked-marginal ensemble. Drug-free median signed within-position Spearman estimates were {controltext}. Adding Evo2 to fixed GEMME+ESM ordinary least squares under position-held-out validation gave a median paired out-of-fold rho increment of {estimate(added)}. Drug-free pooled Q² nevertheless decreased from {baseline_control:.4f} to {augmented_control:.4f}, showing metric-dependent incremental behavior. Three drug-minus-control increments are reported separately in the Results. Intervals describe conditional position-composition sensitivity; they are not validated generalization confidence intervals. Source replacements and deletion sensitivities were retained regardless of direction. Historical codon attribution and structural tests provide limited, distinct context. These analyses quantify endpoint-specific alignment within one selected Fks1 landscape; they do not establish clinical prediction, lethality or a causal constraint mechanism.

## Introduction

An antifungal-target substitution can affect growth without drug, performance during drug exposure, and the difference between those conditions. These phenotypes need not coincide with a model's sequence likelihood. The published Fks1 scan [1] provides matched drug-free and echinocandin conditions, making it possible to compare these components within a bounded target and assay rather than treating a resistance annotation as universal ground truth.

Genomic language models assign probabilities to nucleotide sequence context [2]. A local mutant-minus-reference score can therefore combine amino-acid consequences, nucleotide editing geometry and other sequence preferences. Protein language models provide a comparator at amino-acid resolution [4], while the frozen GEMME analysis supplies an existing evolutionary baseline [7]. We ask how these fixed scores align with each phenotype component and whether frozen Evo2 contributes to a deliberately limited supervised comparison.

Fks1 is the catalytic component of β-1,3-glucan synthase, an echinocandin target [1,6]. This study integrates previously executed project modules and a dated ESM extension. The historical primary benchmark, score composition and orientation analyses, codon reliability analysis, hierarchical condition-slope model, structural context tests and codon-attribution analyses are reused with their original definitions. They are not independent replications of one another. The new comparison protocol was locked before new ESM–phenotype associations, after historical results were known; the entire project is not presented as preregistered confirmation.

## Methods

### Source identity, units and phenotype panels

Original DMS sources were pinned to Landrylab/Durand_et_al_2026 commit d08f8da6d4925ef0847e23140de4eb0a204c6092. Primary observations are unique amino-acid substitutions, not codon realizations or repeated rows. RAW435 preserves the original project MASTER and phenotype-component tables. The historical source audit reproduced high-confidence exact-codon medians within AA/replicate and medians across replicates. The source publication uses both averaging language and median figure-caption language; the computational lineage is stated from pinned code and actual tables, not inferred from prose alone. The gyōza source version is 1.2.1, resolved commit 2b690e2fe78872b98dbdd0700b5bfc5942a9e726 [3]; this does not reconstruct its historical installed environment.

Drug-specific endpoints are the frozen raw differences fitness_drug minus fitness_none for anidulafungin, caspofungin and micafungin. Drug-free fitness is the new extension's primary endpoint; these three differences are secondary. Total drug fitness is descriptive. No performance-dependent scaling, sign choice, imputation, classifier or clinical breakpoint was introduced.

AUTHOR_REPLACEMENT435 is a separate AA-aggregate source sensitivity. Only V641W and L642K use the author validation-derived estimates in each condition, including control; all differences are recomputed. The original rows and historical repeated-measure fits remain unchanged. The exact correction ledger retains raw/replacement values, sequence-plus-condition joins and source hashes. Classified-source annotations retain their author replacement provenance and are not independent validation of raw scores. Operational classified-table identities include WT references and stop variants; the article's wording of 465 single mutants is retained with that table/version caveat rather than interpreted as 465 missense observations.

### Frozen genomic and evolutionary scores

Evo2 is the existing hosted NVIDIA 7B, coding-orientation ALT-minus-REF local53 score, with a 2053-nt context and the correct autoregressive shift. The primary AA score is the frozen median over experimentally represented exact mutant codons. Context, direction and checkpoint were not selected from new performance. The immutable hosted checkpoint SHA remains UNVERIFIED; the existing overlap/sentinel bridge is numerical consistency evidence, not immutable-weight identity. GEMME uses the original raw normalized combined score without sign reversal. Finite AA coverage is distinguished from within-position correlation estimability; constant alternative GEMME components must not be described as missing variants.

### Reused hierarchical and resolution modules

The frozen NEW-2 specification defines 3884 repeated AA/condition/replicate measurements, not 3884 independent mutants. It includes position×condition fixed terms, sample-cell nuisance contrasts, within-position Evo2 predictors, condition-specific position random slopes, nested AA and AA×condition random intercepts, and hotspot×condition residual variance groups. The original maximum-likelihood full/reduced comparison supplies the single three-degree-of-freedom omnibus test; REML supplies effect estimates and model-based Wald intervals. Original rank-scale and 1740-row source-aggregate sensitivities are reused rather than refitted. They share the experiment and are not independent validation. Gaussian residual behavior and limited position/hotspot support remain limitations.

Methods 1 and 2 retain their frozen score-component identity and physical-window orientation definitions; components' explanatory metrics are not additive and rank concordance is not strand invariance. A targeted source-unit audit compares every Method 3 exported exact-nucleotide/condition/replicate coefficient with the original NT source and checks its within-AA residual arithmetic, without refitting its correlations or uncertainty. This distinguishes nucleotide-level source coefficients from copied AA aggregates while retaining shared-assay and measurement-reliability limitations. The audit and original summaries are separate release artifacts.

### Actual ESM-1v scoring

The verified genomic BY4741-derived Fks1 protein was reused after byte/hash checks and comparison with the authoritative query; no alternative strain sequence was substituted. All target WT residues and prepared jobs were checked. Fixed checkpoints esm1v_t33_650M_UR90S_1 through _5 were loaded sequentially with fair-esm 2.0.0 and real CUDA execution. For a canonical one-based target p1, p0=p1−1, w=min(1022,L), start0=max(0,min(p0−511,L−w)), end0=start0+w, and target token=p0−start0+1 for BOS. Windows preserve the WT context and mask exactly one target. Scores are log P(ALT)−log P(WT) at that mask [4], averaged arithmetically across all fixed checkpoints; checkpoint spread is not biological uncertainty.

Real-weight smoke checks preceded full inference: model.eval and float32 inference mode, actual token capacity and logits, finite canonical-AA scores, exact WT-self zero, window/BOS indexing, repeated numerical agreement and comparison with the official masked-marginal formula. Actual downloaded file SHA256, acquisition URL, package versions, device and precision are recorded. Each position/member was checkpointed under input/code/parameter/software/output identity. The same forwards provide theoretical AA scores without creating experimental phenotypes for unmeasured variants.

### Matched signed zero-shot comparisons

Canonical AA keys are joined one-to-one. Each endpoint uses the shared finite mask across Evo2, GEMME, ESM and phenotype. A within-position Spearman correlation requires at least three matched AA observations and nonconstant predictor and outcome vectors; ties use average ranks. Undefined correlations remain NA. Main summaries use positions estimable for every compared predictor and report their counts. The summary is median signed within-position rho; paired contrasts are medians of position-level rho differences, not differences between medians. All planned ESM−GEMME, ESM−Evo2 and Evo2−GEMME contrasts are reported. Pooled correlations and hotspot descriptions are descriptive.

### Fixed supervised comparison and stability

Unweighted ordinary least squares compares intercept+GEMME+ESM with the same model plus frozen Evo2. Every AA at the held-out physical position is excluded from training. Predictor means and population SDs are calculated in the training fold only; a constant feature uses scale one and is flagged. Deterministic least squares retains rank-deficient solutions with singular-value/conditioning diagnostics. No tuning, interaction search, held-out preprocessing or codon phenotype duplication is performed. The primary incremental effect is the median paired position-level OOF rho difference. Secondary MSE, MAE and pooled/equal-position Q² use the fold training-mean baseline. Equal-position Q² is a ratio of equal-position mean squared errors, not a mean of fold Q² values. Exactly zero or nonfinite denominators yield NA; negative Q² is retained.

Uncertainty uses {protocol['uncertainty']['B']:,} position draws with replacement within hotspot and seed {protocol['uncertainty']['seed']}, retaining the eligible hotspot composition and paired records. Percentile limits are 95% conditional position-resampling intervals. Default LOPO intervals resample fixed OOF records without refitting training; the intervals exclude upstream phenotype measurement error. Overlapping training folds preclude treating them as independent errors [5]. These intervals are not validated frequentist CIs for prediction on a new landscape, and pointwise intervals do not provide family-wise coverage. No new formal superiority p-values or post-result equivalence margins are used. Separately, every physical position is deleted from the entire analysis and both LOPO models are rebuilt. All HS2 positions are also removed from training and testing in the prescribed sensitivity. Raw/replacement comparisons use explicit paired eligible-position draws and fixed features/models.

### HS2, codon attribution and structural context

The mandatory HS2 core reuses canonical sequence-derived PRE-HC, T0-HC and FINAL membership, designed-codon multiplicity and the previously closed B3 support descriptions. FINAL589 equals T0-HC, not a later survival stage. L1360 is distinguished between theoretical H and observed same-AA informative panels. The bounded continuous PRE→T0 gate remains BLOCKED_INPUT because comparable biological material/replicate relationships and raw count/denominator provenance are not established by the existing exports. No depletion ratio, extra pseudocount or lethality inference is produced. A separate theoretical-H208 score-support feature gate is evaluated independently; any eligible association is exploratory and does not treat codons as independent biological trials or condition only on measured AA survivors.

Historical NEW-3 M0–M3 estimate within-AA Evo2 codon-score attribution, not experimental resistance. All H panels and components are retained. F2 structural non-estimability and valid-only resampling limits remain disclosed. F2B/H3 support-preserving Dirichlet intervals are sensitivity intervals, not sampling CIs, posterior intervals or p-values. H4B supersedes H4A's invalid outcome-dependent success assertion and closes H only. NEW-3I full-standard-code and separate WT-synonymous synthetic panels remain DEFERRED_SCOPE. Historical H3 NPZ arrays and hosted requests were not freshly reconstructed. NEW-1's existing signed structural tests and NEW-2's original hierarchical fits are reused without performance-driven refitting. Experimental PDB context is verified against RCSB; caspofungin context is not substituted for measured atomic binding geometry of the other drugs [6].

## Results

### Verified source panels and model coverage

The primary panel contains {len(raw)} AA substitutions across {raw.aa_position.nunique()} physical positions. Separate author-source replacements changed only the intended two AA substitutions across four aggregate conditions; no repeated measurement was manufactured. ESM scoring completed {coverage['position_jobs']} position/checkpoint jobs, {coverage['checkpoints']} complete members and {coverage['primary_AA']} primary AA matches. Theoretical score coverage is {coverage['theoretical_AA_no_phenotype']} AA substitutions and is reported separately from measured phenotypes. Individual member scores and the checkpoint acquisition manifest are released. Software execution is not independent biological validation.

### Endpoint-matched zero-shot alignment

Drug-free estimates were {controltext}. All secondary D−C estimates and planned paired contrasts are shown below; pooled summaries and total-drug descriptions remain separate source tables. A negative signed relationship is retained as such, not converted into absence of all information or reversed to improve apparent prediction.

{markdown_table(zero,columns)}
{markdown_table(paired,columns)}

### Limited incremental held-position comparison

The control increment was {estimate(added)}. Secondary increments were {drugtext}. These are phenotype-trained, held-position estimates, not zero-shot scores or independent validation across proteins. A negative Q² means greater squared prediction error than the fold training-mean reference; these values are retained.

Rank and squared-error metrics gave different messages. For drug-free fitness, pooled Q² changed from {baseline_control:.4f} to {augmented_control:.4f} after adding Evo2 (ΔQ²={augmented_control-baseline_control:+.4f}), despite the positive median rank increment. Of the four full-panel endpoints, {negative_q2_increments} pooled-Q² increments were negative. This prevents a claim of uniform improvement and is not resolved by preferring the favorable metric.

{markdown_table(increments,columns)}
{markdown_table(q2,['endpoint','contrast_or_model','metric','point','n_positions','n_AA'])}

### Source and position sensitivities

The full author-replacement and HS2-excluded comparison table is supplied in SI and Figure S4. These sensitivities retain the primary RAW435 designation. Exhaustive control-position deletions gave an incremental-rho range of {influences.point.min():+.3f} to {influences.point.max():+.3f}; {int((influences.point>0).sum())} deletion estimates were positive out of {len(influences)}. These are training-rebuilt stability summaries, not independent tests. Full endpoint-specific deletion and metric-change ledgers preserve every deletion and source-panel difference.

Author replacements yielded a control rank increment of {estimate(replacement_control)} and pooled-Q² increment {value_at('AUTHOR_REPLACEMENT435','fitness_none','delta_pooled_Q2'):+.4f}. Excluding HS2 from both training and testing yielded {estimate(nohs_control)} and pooled-Q² increment {value_at('RAW435_NO_HS2','fitness_none','delta_pooled_Q2'):+.4f}. Thus the squared-error increment changes direction under the prespecified hotspot deletion. The source correction and ascertainment sensitivities answer different questions; neither redefines the full-panel primary analysis. The micafungin rank increment changes from {value_at('RAW435','drug_specific_micafungin','median_paired_position_rho_difference'):+.4f} to {value_at('AUTHOR_REPLACEMENT435','drug_specific_micafungin','median_paired_position_rho_difference'):+.4f} under author replacement, with full paired sensitivity intervals in the source table.

### Reused hierarchical, codon and structural evidence

The historical NEW-2 omnibus condition-slope likelihood-ratio statistic was {oldlrt.estimate:.6g}, with its original df=3 p-value {oldlrt.formal_p_value:.6g}. This tests standardized within-position slope heterogeneity under the original repeated-measure model; it does not establish a drug-specific predictive mechanism or remove every control-subtraction coupling concern. Original rank-scale and aggregate sensitivities remain in SI, alongside heavy-tail and assay/position limitations. M5 alternative-component sensitivities retain their original, differing estimable-position counts and all effect directions.

H2 and H3 tables retain all panels/components/models/contrasts. Technical edit geometry explains a portion of the model's own within-AA score variation; Q² components are not additive and do not partition causal biological mechanisms. G1 and tAI remain external reference proxies rather than measured translation rates. Structural-context tests remain uncertain under their original intervals and multiplicity definitions; no new cutoff, absolute-rho reinterpretation or mechanistic simulation was introduced. Method 3 reliability remains a resolution limitation rather than proof of independent codon-fitness truth or absence of codon effects.

The historical F2 M0 ordinary-bootstrap Q² summary had {f2_m0['Q2_valid_replicates']} valid and {f2_m0['Q2_nonestimable_replicates']} nonestimable draws out of {f2_primary['B']}, with {f2_m0['LOPO_nonestimable_fit_events']} nonestimable fold-fit events. These are separate draw/fit-event units; valid-only intervals condition on estimability. Support-preserving H3 sensitivity addresses a different perturbation and does not convert those intervals into sampling confidence intervals. The new source-unit audit verified {lineage['rows_checked']} Method 3 rows and found distinct original NT coefficients in all {lineage['groups_with_distinct_NT_coefficients']} multi-codon AA/condition/replicate groups. It closes the copied-AA-aggregate ambiguity while preserving the historical reliability limitation.

### Separate categorical representation extension

{support_note}

## Discussion

The main contribution is an endpoint-resolved, source-audited comparison within a matched Fks1 landscape. The tables quantify where genomic, evolutionary and protein likelihood scores align with drug-free and selective phenotypes; the actual directions and uncertainty remain visible. A score correlated with drug-free competition can be consistent with general sequence constraint while failing to distinguish specific drug-mediated effects, but correlation alone cannot identify that mechanism. Conversely, negative D−C association can carry structured information about costs or shared-control effects. It is not permission for an outcome-selected sign flip.

Adding Evo2 to a fixed GEMME+ESM model asks a restricted incremental question. Its effect estimate, conditional intervals and exhaustive deletion range must be considered together. A weak or inconsistent increment does not prove equivalence or that Evo2 contains no additional information. A positive estimate or a pointwise interval excluding zero is likewise not a confirmatory superiority test. We do not promote a favorable sensitivity, component or hotspot into a new primary model.

In the present results, the modest positive drug-free Evo2 and GEMME signed-rank estimates coexist with an ESM estimate near zero. The D−C pattern varies across drugs, with negative Evo2/GEMME point estimates for caspofungin and micafungin. This is compatible with functional constraint and shared-control coupling contributing to the landscape, but no direct endpoint-superiority test or causal decomposition was performed. A protein comparator's weak alignment also cautions against describing every sequence score as a common, validated functional-constraint measurement. The new LOPO control rank increment is positive across all individual position deletions, yet its conditional resampling interval crosses zero, full-panel squared error worsens, and the Q² increment changes direction after HS2 removal. These observations address different stability/metric questions and jointly support a limited, metric-dependent incremental interpretation.

Historical codon attribution separates nucleotide-level preferences from protein-level interpretation. The geometry baseline is close to the autoregressive scoring procedure, so explanatory power for an Evo2 score is not a measured translational or resistance mechanism. Expanded theoretical panels characterize scoring and representation, not additional DMS phenotypes. The ESM ensemble adds one comparator; it is not a model-search campaign or five biological replicates.

Ascertainment particularly limits HS2 interpretation. High-confidence support reflects construction, abundance, sampling, thresholds, platform and biological effects jointly. Neither z-scoring nor excluding a hotspot repairs missing-not-at-random selection. The core presentation is therefore retained even when a continuous analysis is blocked, and any categorical score-support extension has a separate eligible universe and exploratory interpretation.

## Limitations

This is one gene in one model-organism landscape, with hotspot and assay/library history partly confounded. LOPO does not establish transfer to another gene, fungal species, patient, or unobserved hotspot class. Phenotype uncertainty is not fully propagated into new rank summaries. Historical modules reuse the same experiment and are not independent replications. LOPO training sets overlap; the reported resampling intervals condition on fixed predictions and observed position/hotspot composition. Source replacements reuse validation-derived data and are not independent validation of those estimates. The hosted Evo2 immutable checkpoint SHA remains unverified. Historical H3 arrays/API calls are not freshly reproduced. Continuous PRE→T0 remains input-blocked. Negative and contradictory findings are retained, and no equivalence, clinical MIC, lethality or causal mechanism conclusion is claimed.

Zero-shot here means no task-phenotype training; it does not mean that the species or homologous sequence was absent from pretraining. Pretraining overlap was not independently quantified in this study.

## Data and code availability

Original source: pinned Landrylab/Durand_et_al_2026 commit. The release contains code/config/protocol lock, source locators and actual hashes, new ESM member/ensemble scores, source/correction/mask ledgers, OOF predictions/fold audits, conditional draws, all-position sensitivities, figure source tables/exports, manuscripts and verification records. Large raw inputs and all persisted scientific outputs remain in authenticated Drive. Model weights were downloaded only to the dedicated Colab GPU runtime; its cache is ephemeral and is not assumed available after that runtime was released. Checkpoint acquisition URLs and actual file hashes are provided for fresh remote acquisition when needed. The source-complete context snapshot is separately hash-addressed. No credentials or authorization state are released. Release verification scope and one actual resume command appear in RUN_INDEX and REPRODUCIBILITY_REPORT. No public deposition or manuscript submission has been performed.

## Administrative fields

{admin}

## References

1. Durand R et al. Mutational landscape and molecular bases of echinocandin resistance in *Saccharomyces cerevisiae*. GENETICS 233(1), iyag055 (2026). https://doi.org/10.1093/genetics/iyag055 . Bibliography/method excerpts checked in the provided original PDF; pinned author code governs computational lineage.
2. Brixi G, Durrant MG et al. Genome modelling and design across all domains of life with Evo 2. Nature (2026). https://doi.org/10.1038/s41586-026-10176-5 . Provided original PDF; hosted immutable-weight identity is a separate unresolved issue.
3. Durand R, Pageau A, Landry CR. gyōza: a Snakemake workflow for modular analysis of deep-mutational scanning data. GENETICS 232(1), iyaf199 (2026; advance access 2025). https://doi.org/10.1093/genetics/iyaf199 . Provided PDF and resolved source commit.
4. Meier J et al. Language models enable zero-shot prediction of the effects of mutations on protein function. NeurIPS 34 (2021). https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html . Official scoring implementation: https://github.com/facebookresearch/esm/blob/main/examples/variant-prediction/predict.py . These sources support the comparator/scoring family, not validation of the present long-protein window.
5. Bengio Y, Grandvalet Y. No Unbiased Estimator of the Variance of K-Fold Cross-Validation. JMLR 5:1089–1105 (2004). https://www.jmlr.org/papers/v5/grandvalet04a.html . Supports caution about naive independent-fold variance, not validity of the present conditional intervals as generalization CIs.
6. Ren Z et al. Structural basis of fungal β-1,3-glucan synthase inhibition by caspofungin. Nature 654:547–555 (2026). https://doi.org/10.1038/s41586-026-10409-7 . Entry metadata checked at https://www.rcsb.org/structure/9PE1 , https://www.rcsb.org/structure/9PE3 and https://www.rcsb.org/structure/9ZTC ; not a claim that every page of the full text was reread.
7. Laine E, Karami Y, Carbone A. GEMME: A Simple and Fast Global Epistatic Model Predicting Mutational Effects. Molecular Biology and Evolution 36(11):2604–2619 (2019). https://doi.org/10.1093/molbev/msz179 . Primary publication metadata checked; the project uses its existing frozen webserver matrices, not a newly fitted GEMME model.
'''
    (attempt/'MANUSCRIPT.md').write_text(manuscript)
    si='# Supplementary information\n\n'+admin+'\n\n'
    si+='## Locked extension methods\n\n'+json.dumps(protocol,indent=2)+'\n\n'
    for heading,frame in [('All new comparison summaries',summary),('Historical NEW-2 primary and sensitivities',n2),('Historical H2 all points',h2),('Historical H3 all contrasts',h3),('Historical M5 robustness (different position denominators retained)',m5),('Historical structural tests',structure),('Historical Method 3 reliability',reliability),('Source correction ledger',ledger)]:
        si+='## '+heading+'\n\n'+markdown_table(frame)+'\n'
    si+='## Additional source tables\n\n'+markdown_table(pd.DataFrame(table_index))+'\n\nFull per-AA OOF predictions, every fold and position deletion, all draws and every figure source are released as machine-readable files rather than represented as independent replicate observations.\n'
    si+='## HS2 and NEW-3 disposition\n\nContinuous PRE→T0: BLOCKED_INPUT. Categorical H208 score-support: see its independent gate/output, not the continuous count gate. NEW-3I: DEFERRED_SCOPE. H3 intervals remain support-preserving sensitivity intervals. No all-historical-refit or full-study-preregistration claim.\n'
    (attempt/'SI.md').write_text(si)
    (attempt/'INTEGRATED_RESULTS.md').write_text('# Executed results\n\n'+controltext+'\n\nControl LOPO increment: '+estimate(added)+'\n\nSecondary increments: '+drugtext+'\n\n'+markdown_table(increments,columns)+'\nNew versus reused and limitations are specified in MANUSCRIPT.md, CLAIM_LEDGER.csv and TABLE_INDEX.csv.\n')
    (attempt/'BIBLIOGRAPHY.md').write_text(manuscript.split('## References\n\n')[1])
    docx_status={}
    try:
        from docx import Document
        for name in ['MANUSCRIPT','SI']:
            document=Document()
            for line in (attempt/(name+'.md')).read_text().splitlines():
                if line.startswith('#'):document.add_heading(line.lstrip('# ').strip(),level=min(len(line)-len(line.lstrip('#')),3))
                elif line.strip():document.add_paragraph(line)
            document.save(attempt/(name+'.docx'))
        docx_status={'status':'EXECUTED_VALIDATED','format':'editable draft; Markdown tables remain plain text in DOCX'}
    except ImportError:
        docx_status={'status':'NOT_APPLICABLE','reason':'python-docx absent; complete editable Markdown supplied'}
    atomic_json(attempt/'WRITEUP_STATUS.json',{'status':'MANUSCRIPT_DRAFT_COMPLETE','docx':docx_status,'author_admin_fields':'AUTHOR_INPUT_REQUIRED','new_ESM_coverage':coverage,'new_results_generated_from_actual_tables':True,'manuscript_numbers_source_ledger':'CLAIM_LEDGER.csv'})
    return [p.name for p in attempt.iterdir() if p.name!='started.json']
