# Analysis scripts

The scripts are numbered in the intended pipeline order. Each module is executable as a command-line program and writes machine-readable outputs for the next analytical stage or for final manuscript assembly.

| Script | Purpose | Main outputs |
|---|---|---|
| 01_data_harmonization.py | Harmonize de-identified analysis-ready tables; preserve missingness; create protocol z-scores | `outputs/harmonized/*_analysis_ready.csv`, sample-flow CSVs |
| 02_family_graph_reconstruction.py | Build strict social graph, Mendelian-compatible genetic layer, relation subgraphs and parent-defined lines | `outputs/family_graph/*.csv` |
| 03_age_genotype_models.py | HTR1A/HTR2A genotype × stage and genotype × adult-age models | `outputs/age_genotype_models.csv` |
| 04_parent_offspring_permutations.py | Directed mother/father phenotype permutations with covariate adjustment | `outputs/parent_offspring_permutations.csv` |
| 05_family_graph_permutations.py | Component-equal family-clustering tests; max-sample and age-adjusted models are emitted together | `outputs/family_graph_permutations.csv` |
| 06_kinship_scale_decomposition.py | Locus × phenotype tests across family scales | `outputs/kinship_scale_decomposition.csv` |
| 07_multilocus_and_dispersion.py | Multilocus block, interactions, dispersion and genetic-distance QC | `outputs/multilocus_dispersion/*.csv` |
| 08_lasso_analysis.py | Nested-CV LASSO, λmin/λ1se, stability selection | `outputs/lasso/*.csv` |
| 09_grouped_random_forest.py | Leave-one-family-component-out RF and family-level sign-flip tests; max-sample and age-known cohorts are emitted together | `outputs/random_forest/*.csv` |
| 10_generate_tables_figures.py | Assemble manuscript tables, Supplementary registries and selected figures | `outputs/manuscript/*` |

## Reproducibility safeguards

- Family IDs/components are never used as predictive features.
- Parent→offspring direction is taken from explicit de-identified fields; it is not guessed from unordered legacy pairs.
- Stochastic procedures use fixed default seeds.
- LASSO preprocessing is fit only inside training folds.
- RF genotype imputation is fit only on training families.
- Public scripts do not require participant names or raw field identifiers.

## Final release check

Before tagging `v1.0.0`, run the full pipeline with the frozen Supplementary Data package and compare generated manuscript-level values against the frozen tables in the submitted manuscript. Any discrepancy should be resolved before the Zenodo archive is published.
