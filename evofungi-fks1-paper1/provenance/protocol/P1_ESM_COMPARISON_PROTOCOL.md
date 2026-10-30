# Adopted extension protocol — 18 October 2026

# New ESM comparison protocol and claim rules

This specification supplies bounded defaults for the previously open new-analysis choices. Adopt it in a dated live protocol BEFORE computing new ESM–phenotype associations or performance. These are extension decisions, not edits to historical freezes. Identity and schema validation are allowed before lock. Existing historical effect estimates, CIs, p-values and multiplicity families remain as originally defined.

A mathematical impossibility or documented source incompatibility may require a pre-result technical amendment. Weak results are not a reason to select another metric or model. If new associations were already viewed in the live project, record that chronology rather than claiming outcome blindness.

## 1. Questions and endpoints

- Q1: On the observed 435-AA panel, what signed within-position association does each fixed Evo2, GEMME and ESM score have with drug-free fitness? What are the paired performance differences? This benchmarks one landscape, not a clinical task.
- Q2: Repeat fixed comparisons for each of the three historical drug-specific D−C endpoints. These are secondary. Total drug fitness is descriptive, not a replacement primary.
- Q3: Does adding frozen Evo2 improve position-held-out predictions over GEMME+ESM? Drug-free fitness is primary; three D−C endpoints are secondary. This stage fits phenotype and is not zero-shot.

RAW435 is primary. AUTHOR_REPLACEMENT435 is source sensitivity. Use the same feature matrices/masks where outcomes are finite. Do not select a panel by performance. Replacements can change both fitness and D−C; their comparison is not independent validation of the fitted model.

## 2. Predictor and ESM contract

Evo2 primary AA score and GEMME combined primary keep their existing definitions and directions. Confirm sign and codon-to-AA aggregation from original derivations, not observed correlations. Independent/epistatic GEMME are existing context, not candidates for choosing a better baseline.

For substitution WpA, ESM score is logP(A | WT context with position p masked) minus logP(W | that same masked context). Average these scores across the five fixed ESM-1v checkpoints. This is a mean log ratio, not a ratio formed from averaged probabilities. Checkpoint spread is not biological measurement uncertainty.

The source-verified protein length is currently 1876 AA. Define p0=p1−1 for canonical one-based position p1. With w=min(1022,L), use:

    start0 = max(0, min(p0 - 511, L - w))
    end0 = start0 + w
    target_local0 = p0 - start0
    target_token = target_local0 + 1  # verified BOS

The end coordinate is Python-exclusive. Compare this with existing 140 job records and the live code to preserve the 17 October convention. Do not substitute p1 for p0 and shift the window by one. If an actual earlier frozen job explicitly differs, resolve that source identity before inference and without performance-based selection. No residue padding or truncation; verify BOS/EOS and actual checkpoint capacity.

Use eval/inference mode and record precision/device. No fine-tuning, embedding probes or mutant structure prediction. One masked WT forward per position/checkpoint supplies all 20 canonical AA values. WT-self zero is an arithmetic check, not biological validation.

## 3. Keys, common masks and estimability

Use canonical AA ID as row key, with validated position/hotspot mapping. Join one-to-one, never by row order. For each endpoint, use the shared finite-AA mask across all three primary predictors and the outcome. Both nested supervised models use the same mask, including Evo2 availability even for the baseline that omits Evo2. Log exclusions and mask hash. No imputation.

Expected coverage is 435 AA and 28 positions, but actual ESM outputs must prove it. Incomplete checkpoint coverage is a technical failure. Any justified biological reduction requires a recorded amendment; a reduced analysis is not the full 435 benchmark.

Within-position Spearman requires at least three matched AA observations and nonconstant predictor/outcome vectors. Use average ranks for ties. Undefined correlations are NA, not zero. Main paired zero-shot summaries use positions estimable for all three models on the shared AA panel. Also publish each model's available-position results and reasons for exclusions. For supervised contrasts, both OOF prediction vectors and the outcome must be estimable at the position. Publish these rules before association and report counts.

The known 28/24/23 GEMME-component difference concerns correlation estimability, not finite AA coverage.

## 4. Zero-shot summaries

For each endpoint and position, calculate signed Spearman rho for each model. Main summary is the median within-position rho on the paired position panel. The paired effect for A versus B is:

    median over positions of [rho(A,p) - rho(B,p)]

This is not generally median(rho A) minus median(rho B). Release component correlations and paired differences.

Planned contrasts are ESM−GEMME, ESM−Evo2 and Evo2−GEMME. Report all three for control and all three for each D−C endpoint; select none after results. Pooled Spearman on the shared AA mask is descriptive and may include between-position variation. Publish hotspot counts and stratified descriptions. Do not combine all conditions into one headline benchmark.

Retain signed correlations. No absolute-rho ranking or performance-based sign flip. Only a source-documented direction convention fixed before association may transform a comparator, preserving raw and analysis columns. Historical Evo2 sign remains unchanged.

## 5. Fixed limited supervised comparison

Baseline is intercept + GEMME combined + ESM ensemble. Augmented is the same plus Evo2 primary AA score. Use unweighted OLS on AA rows, as prescribed in 17 October READINESS_DECISIONS. No hyperparameter tuning, spline, new interaction, hotspot covariate or feature search.

Use physical-position LOPO. Exclude every AA at the held-out position from training. No held-out outcome may enter preprocessing, transformations, missingness treatment or model selection. Do not duplicate AA phenotype into codon/replicate observations.

Compute predictor means/SDs from training only. For a constant training feature, use scale=1 and flag it. Use deterministic least squares and record rank, singular values and conditioning. Preserve rank-deficient folds with diagnostics; do not silently drop features or switch model. Numerical overflow is a technical failure, not a reason for outcome-driven winsorization. The fold baseline is the training-outcome mean.

Primary new incremental effect on control is:

    median over eligible positions of
    [rho(augmented OOF,p) - rho(baseline OOF,p)]

Three corresponding D−C differences are secondary. Report both models' absolute OOF performance. This equal-position rank summary avoids letting positions with more measured AA dominate.

Secondary summaries: pooled OOF MSE, MAE and Q² relative to the fold training-mean baseline; also equal-position MSE and corresponding Q². Let e_i=(y_i−prediction_i)^2 and b_i=(y_i−training_mean_i)^2:

    pooled Q² = 1 - sum(e_i)/sum(b_i)
    equal-position Q² = 1 - mean_p(mean_i_in_p(e_i)) /
                              mean_p(mean_i_in_p(b_i))

The latter is not mean per-position Q². If the denominator is zero or numerically undefined under a documented fixed tolerance, return NA. Do not choose an epsilon from observed performance. Retain negative Q². Both models share the same denominator; delta Q² is augmented minus baseline. These are secondary metrics, not replacements when rank gain is weak.

Save per-AA OOF predictions, fold baselines, outcomes, source panel and fold ID. Save training/test AA IDs, preprocessing and rank audits. Training fit must never be presented as held-out accuracy.

## 6. Uncertainty and position influence

Keep historical uncertainty definitions unchanged. For this bounded extension, use **effect estimation with explicitly conditional position-resampling intervals and influence analysis**, not new confirmatory superiority tests.

Use B=10,000 draws and seed 20261018. Sample physical positions with replacement within each hotspot, retaining the original number of eligible positions per hotspot. Preserve paired models/endpoints, all rows of each position and raw/replacement pairing. Save the draw matrix or exact reproducible draw definition and per-draw summaries. If paired masks differ by endpoint, record separate eligible-position sets and how draws are generated; do not reuse invalid IDs.

For zero-shot scores, recompute medians/differences from fixed score/outcome position records. For the default LOPO intervals, resample fixed OOF position records; do not claim that training was refitted. The 2.5th/97.5th percentiles form the 95% conditional resampling interval.

State that this conditions on the observed panel, computed OOF predictions and hotspot composition. It excludes upstream phenotype measurement uncertainty. Do not independently bootstrap 435 AA or compute independent-fold SE/sqrt(28).

**CV limitation:** LOPO training sets overlap. The default interval describes sensitivity to observed position composition conditional on existing OOF predictions. It is NOT a validated frequentist CI for an algorithm's generalization error on a new landscape. Do not use zero exclusion as a formal test or claim independent replication. Three hotspots do not support a reliable broad population inference from a three-cluster bootstrap. Put this limitation in Methods and captions.

Assess training/position influence by exhaustive deletion of each physical position from the entire new analysis. For every deletion, rebuild LOPO on the remaining positions, refit both models and recompute the primary effect. Include all positions, not just L1360 or influential favorable cases. These are stability summaries, not independent tests or CIs. Also rerun the fixed comparison with all HS2 positions excluded from both training and testing. Report ranges and sign stability without optimizing exclusions.

Formal CV generalization inference would require a separately justified group-aware method and validation; it is not silently added here. This package explicitly chooses a restrained effect-estimation route for the limited Paper 1 extension. Bootstrap is not being claimed to solve CV uncertainty generally.

## 7. Multiplicity and permissible claims

No new formal p-value is required for the ESM extension. Publish all planned contrasts/endpoints with point estimates and correctly labeled conditional intervals. Multiple pointwise intervals do not have family-wise 95% coverage. Do not announce superiority from an unadjusted interval excluding zero or highlight only favorable contrasts. Historical families retain their original adjustments. If a valid later live protocol already defines formal tests, reconcile its authority before changing it; do not mix two approaches silently.

Five checkpoints are not five biological replicates. Do not infer statistically established control-versus-drug superiority from separate endpoint significances. Introduce no equivalence margin after results. 'No clear consistent incremental gain detected' can be considered against effect, uncertainty and stability; 'Evo2 adds no information' is too strong.

## 8. Required sensitivities and optional data

AUTHOR_REPLACEMENT: rerun new zero-shot and LOPO summaries using the two-AA replacements across all four conditions, with fixed masks, models, folds and paired uncertainty draws. Do not refit the historical NEW-2 replicate model with copied aggregates. Quantify metric changes rather than assuming two variants are negligible.

HS2 exclusion: remove all HS2 positions and rebuild training/testing on the remainder. Full panel remains primary regardless of sensitivity direction. Hotspot-specific descriptions have small-group and ascertainment limitations. L1360 membership in the codon panels is separate from its influence in the new AA comparison.

Do not use confidence labels as inverse-variance weights, apply automatic disattenuation, or choose alternative GEMME components from performance. Optional FKS2-repressed Fks1 and individual-validation records are nonblocking checks only when the existing intake provides ready locators. They require documented units and overlap limitations. Do not count the validation data used for replacements as independent validation of those replacements.

## 9. Primary methodological references

Verify metadata and claims when writing; these references do not validate our implementation automatically.

- Meier et al. (2021), Language models enable zero-shot prediction of the effects of mutations on protein function. Official implementation: https://github.com/facebookresearch/esm/tree/main/examples/variant-prediction . Supports the model/scoring family, not clinical interpretation or validation of our long-protein window.
- Bengio and Grandvalet (2004), No Unbiased Estimator of the Variance of K-Fold Cross-Validation. JMLR 5:1089–1105. https://www.jmlr.org/papers/v5/grandvalet04a.html . Supports avoiding naive independent-fold variance, not a guarantee for our bootstrap.
- Bates, Hastie and Tibshirani, Cross-validation: what does it estimate and how well does it do it? DOI:10.1080/01621459.2023.2197686. https://doi.org/10.1080/01621459.2023.2197686 . Distinguishes CV targets and uncertainty; its methods are not automatically implemented here.
- Roberts et al. (2017), Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. DOI:10.1111/ecog.02881. Supports structured holdout reasoning, not transfer proof.
- Durand et al. (2026), original DMS, DOI:10.1093/genetics/iyag055. Bundled PDF and pinned author code define phenotype provenance.


## Implementation lock

The accompanying MACHINE_READABLE_PROTOCOL.json fixes actual column names, float32 GPU arithmetic, 1e-5 numerical repeat/formula tolerances, exact-zero/nonfinite Q² denominator handling, seed and source direction conventions before new associations. SHA256: 168bc66ad0646f338d00d7edd58832ecb1e912a38c31a4513fbdb0d75c328ae7. Existing historical results were known. This is not whole-study preregistration.
