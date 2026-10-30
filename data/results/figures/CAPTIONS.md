# Figure captions

## Fig1_ascertainment

Exact-genotype and AA-group representation across the targeted design. Top: exhaustive NNK missense design, source PRE high-confidence and T0 high-confidence support; FINAL589 equals T0-HC and is not a later survival stage. Bottom: designed AA groups, represented AA groups and groups retaining multiple exact codons. L1360 remains represented in the theoretical design but has no observed informative multi-codon AA group. Canonical positions derive from sequence. Support loss is not lethality. No uncertainty bars: exhaustive design descriptions; historical B3 reused.

## Fig2_matched_benchmark

RAW435 matched zero-shot benchmark. Equal-position median signed Spearman correlations, using positions estimable for all three predictors on the same finite AA panel (counts beneath points). Lines are pointwise 95% conditional hotspot-stratified position-resampling intervals (10,000 draws); they condition on the observed panel and exclude phenotype measurement uncertainty. They are not formal superiority tests or new-landscape generalization CIs. Drug−control endpoints are secondary.

## FigS1_paired_zero_shot

All planned paired zero-shot contrasts. Median of position-level rho differences, not difference of medians. Conditional resampling intervals as in Figure 2; all planned contrasts/endpoints shown without outcome selection.

## Fig3_slopes_and_increment

Left: reused NEW-2 condition-specific standardized within-position slopes, with original model-based 95% Wald intervals; repeated-measure hierarchical model, not refitted. Right: new supervised position-held-out OLS comparison, GEMME+ESM+Evo2 minus GEMME+ESM, summarized as median paired OOF rho difference. Right-hand intervals resample fixed OOF position records; overlapping LOPO training folds and upstream phenotype uncertainty prevent interpretation as validated generalization CIs. Panels answer different questions and have different units.

## Fig4_NEW4_landscape

Original NEW-4 functional-fitness × drug-selective-effect landscape, RAW435. One point per AA substitution. Three drug−control outcomes versus drug-free fitness, with plain phenotype and fixed model overlays. All values shown; each score uses a consistent full-range scale across drugs. Zero axes do not define clinical categories. No classifier, optimized threshold, clinical resistance or causal trade-off inference. Author replacements are a separately labeled source sensitivity.

## FigS2_author_landscape

Original NEW-4 functional-fitness × drug-selective-effect landscape, AUTHOR_REPLACEMENT435. One point per AA substitution. Three drug−control outcomes versus drug-free fitness, with plain phenotype and fixed model overlays. All values shown; each score uses a consistent full-range scale across drugs. Zero axes do not define clinical categories. No classifier, optimized threshold, clinical resistance or causal trade-off inference. Author replacements are a separately labeled source sensitivity.

## Fig5_codon_and_structure

Reused historical evidence, not newly refitted. H2 predicts within-AA Evo2 score variation rather than experimental fitness; G0 edit geometry is technical attribution. H3 bars are support-preserving Dirichlet sensitivity intervals, not sampling CIs, posteriors or p-values. NEW-1 shows the original signed structural-context tests, preserving uncertain/negative findings and original bootstrap intervals. Caspofungin structural context does not establish equivalent binding geometry for the other drugs. No new structural features or mutant simulations.

## FigS3_position_influence

Exhaustive single-position deletion stability. Each point removes the named position from training and testing and rebuilds both LOPO models. All positions are shown; these are stability summaries, not independent tests or confidence intervals. Full RAW435 remains primary.

## FigS4_fixed_sensitivities

All prespecified new LOPO sensitivity panels/endpoints: primary raw, author replacement, and removal of all HS2 positions from both training and testing. Fixed features/models and shared AA masks within comparisons. Conditional fixed-OOF intervals as in Figure 3; no sensitivity was selected as a new primary.

## FigS5_codon_representation

Categorical representation of the theoretical multi-codon AA universe, separated by hotspot. Designed codons and observed support are distinct units. These reused B3 descriptions do not quantify absolute survival or continuous depletion; the H208 feature gate is separate.

## FigS6_historical_robustness

Reused M5 and Method 3 summaries with their original bootstrap intervals. The M5 combined/independent/epistatic estimable-position denominators differ and are printed explicitly; a favorable alternative component is not substituted as the primary. Method 3 describes the existing within-AA measurement-resolution check, not an independent codon-fitness ground truth or new clinical validation. No historical refit or uncertainty recalculation.
