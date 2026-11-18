# Fks1 fitness landscape summary data

Summary tables and figure sources for *Genomic and protein sequence scores across drug-free and echinocandin-selected Fks1 fitness*, by Muhammed Ali Boz. Manuscript snapshot: 18 November 2026.

The tables describe a retrospective analysis of 435 amino acid substitutions at 28 positions in laboratory *Saccharomyces cerevisiae*. C denotes drug-free growth, D denotes total growth under drug, and D−C is the difference between those readouts.

| Files | Contents |
| --- | --- |
| `tables/Table_1.csv` | Median within-position Spearman correlations for five fixed scores across seven outcomes |
| `tables/Table_2.csv` | Paired differences in signed and absolute Spearman correlations, with conditional intervals |
| `tables/Table_3.csv` | Leave-one-position-out regression comparisons and conditional uncertainty summaries |
| `tables/*_FIGURE_SOURCE.csv` | Aggregate sources for the five main figures |
| `tables/U9_CALIBRATION.csv`, `tables/U9_NULL_DRAWS.csv` | Conditional calibration summaries and simulated likelihood-ratio statistics |
| `render_figures.py` | Renderer for main Figures 1–5 using the frozen summary inputs |

From this directory:

```sh
python -m pip install -r requirements.txt
python render_figures.py
```

The renderer reads the frozen summaries and writes PDF, PNG and SVG figures to `figures/`. Supporting Data S1 and S2 accompanying the manuscript document the scientific calculations behind these summaries. The repository root also documents the October replay and its dependencies.

The original experimental study is Durand et al. (2026), [Genetics, iyag055](https://doi.org/10.1093/genetics/iyag055). The historical analysis is pinned at [commit d215732](https://github.com/bokve/evofungi-fks1-paper1/tree/d21573281a20fd5844a7b1dd9e5776a3b4ec11ad). The renderer is under MIT. The summary tables and this documentation are under CC BY 4.0.
