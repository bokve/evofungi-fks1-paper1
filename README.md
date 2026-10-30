# Fixed sequence scores in a yeast Fks1 fitness landscape

Code and frozen numerical records for the computational reanalysis by Muhammed Ali Boz, prepared for Computational Biology and Chemistry. Release candidate v1.0.0, 29 October 2026. This is not a clinical resistance classifier or a new experimental study.

The primary panel contains 435 amino acid substitutions, represented by 589 exact genotypes, at 28 selected positions in one Saccharomyces cerevisiae protein. Fixed Evo2, GEMME and five-member ESM-1v scores are evaluated against drug-free fitness and separately against drug-minus-control outcomes. Zero-shot scores and supervised score combinations are distinct analyses.

In the drug-free panel, adding Evo2 gives a numerical median paired within-position rank increment of 0.108, with a conditional interval of -0.017 to 0.133, while pooled predictive Q2 decreases from 0.0289 to 0.0139. These findings do not establish consistent incremental benefit or transfer to another protein/species. The theoretical panels are not additional fitness observations. Continuous PRE-to-T0 inference remains unperformed.

## Offline replay

Use Python 3.12 in an isolated environment. The dependency file records the environment tested for this release replay, not the historical GPU/CPU generating environment; the latter is retained in provenance/SOFTWARE_LOCK.json. No account, API key, GPU or model download is needed for these commands.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-replay.txt
python -m unittest discover -s tests -v
python scripts/replay_frozen_results.py
python scripts/rebuild_cbc_figures.py
```

The replay checks the panel, eight source replacement cells, shared predictor identities, the five-member ensemble, 384 saved summary points, 1,152 fixed-score position correlations, 144 conditional intervals and 768 saved fold-audit records. It recomputes summaries from fixed scores and saved OOF predictions; it does not refit models or generate new biological predictions. All 30 existing unit tests pass after restoring three original helper modules needed by test_readiness.py. The tests include artificial fixtures for leakage, input identity and source replacement. Tests and saved-record replay do not validate a biological mechanism.

For historical R diagnostics and H3 array summaries, extract Paper1_Scientific_Archive_v1.0.0.zip to a directory called scientific_archive, then run:

```bash
python scripts/verify_historical_records.py --data-root scientific_archive
```

The archive root must directly contain results/ and legacy/. This checks the recorded likelihood-ratio arithmetic, saved residual summaries and 405 H3 quantile statistics. It does not rerun the historical R or H3 fitting procedure.

## Package map

| Path | Purpose |
| --- | --- |
| paper1/ | Original analysis modules, including historical plotting code |
| tests/ | Original 30-test suite |
| tools/ | Three original helpers restored to satisfy the test imports |
| reference_audit/evidence/ | Small original test fixtures and source records |
| data/results/ | Fixed panels, scores, OOF predictions, interval draws and tables |
| scripts/ | Offline replay and CBC figure rebuilding entry points |
| provenance/ | Protocol, checkpoint/software records and source license |
| FIGURE_MAP.csv | Current CBC figure names mapped to original plot outputs |
| FILE_MANIFEST.json | SHA256 and byte count for release files |

The accompanying scientific archive additionally contains legacy R scripts, model-ready tables, fitted objects and H3 saved arrays; selected frozen evidence; and the full-precision CBC source tables. Its file inventory records exact bytes. The original mixed working ZIP and source-context notebook are not published wholesale; their original identities are retained in the provenance record.

## Reproduction boundaries

The original model-fitting modules use historical configurations and paths. The original remote launcher is not the offline replay entry point. Do not run historical notebooks with Run All. Generating fresh Evo2/ESM/GEMME scores or refitting R/H3 models is outside the tested replay command. The hosted Evo2 immutable checkpoint identity remains unverified; stored scores can be replayed, but exact historical inference is not claimed. GEMME frozen values are included; original web-server output and MSA regeneration are not independently certified. Conditional position-resampling intervals hold OOF predictions fixed and do not estimate full training/generalization uncertainty.

Original experimental source: https://github.com/Landrylab/Durand_et_al_2026/tree/d08f8da6d4925ef0847e23140de4eb0a204c6092

Read RIGHTS.md for license scope and third-party attribution. The public source repository is https://github.com/bokve/evofungi-fks1-paper1. The companion scientific-data record has a reserved version DOI: https://doi.org/10.5281/zenodo.23814961. This DOI is registered and becomes publicly resolvable only when the Zenodo dataset is published. GitHub code availability alone does not establish that the software release has been archived by Zenodo.

The CBC plot renderer recreates the current figure compositions from their exact frozen source CSVs. Original statistical plotting modules are retained; font, layout and file metadata can differ between renders. The renderer does not change any scientific values.
