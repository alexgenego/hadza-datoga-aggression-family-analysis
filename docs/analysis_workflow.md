# Analysis workflow

The public pipeline is organized as ten sequential analysis modules.

1. **01_data_harmonization.py** — harmonizes de-identified child/adult phenotype, demographic and genotype tables while preserving original rows and missingness; creates protocol-standardized phenotype variables and sample-flow summaries.
2. **02_family_graph_reconstruction.py** — builds the strict social/field graph and the nested Mendelian-compatible genetic layer; writes connected-component labels, relation-specific edge tables, and parent-defined offspring lines when explicit parent/offspring direction is available.
3. **03_age_genotype_models.py** — evaluates HTR1A/HTR2A genotype × developmental/protocol-stage effects and adult genotype × continuous-age effects.
4. **04_parent_offspring_permutations.py** — tests covariate-adjusted parent–offspring phenotype similarity with directed family-preserving permutations; optional sex-sensitive permutations preserve offspring sex.
5. **05_family_graph_permutations.py** — tests component-level family clustering using equal-component mean pairwise squared residual distance; emits max-sample and age-adjusted phenotype models in one run.
6. **06_kinship_scale_decomposition.py** — compares locus-specific associations across sibling components, mother/father-defined offspring lines, direct parent-genotype tests, mixed-generation components, and the full graph; q9 is calculated within phenotype × scale screens.
7. **07_multilocus_and_dispersion.py** — evaluates additive multilocus-block ΔR², block × age/population follow-up tests, Brown–Forsythe genotype-group dispersion, and complete-profile genetic-distance QC.
8. **08_lasso_analysis.py** — nested cross-validated LASSO with fold-specific covariate handling, λmin/λ1se selection, 500 stratified 80% stability-selection subsamples, and out-of-fold ΔR² relative to the covariate-only baseline.
9. **09_grouped_random_forest.py** — leave-one-family-component-out Random Forest for both max-sample and age-known cohorts, with equal total training weight per family, training-fold-only genotype imputation, family-balanced error metrics, weighted R², and paired Monte Carlo sign-flip comparison.
10. **10_generate_tables_figures.py** — downstream-only assembly of manuscript tables, Supplementary registries, provenance manifests, and selected figures from outputs of modules 01–09.

## Family-scale principle

The analysis treats the definition of the family unit as part of the biological question. Effects observed among adult offspring from the same parent-defined line are not assumed to be equivalent to effects obtained after mixing parents and offspring or averaging over the full connected pedigree.

## Key null models

- Full-graph phenotype clustering: phenotype residuals are shuffled within population while graph/component membership remains fixed.
- Relation-specific parent–offspring analysis: directed parent structure and offspring counts are preserved; admissible offspring identities are permuted within population, with optional preservation of offspring sex.
- Kinship-scale locus tests: genotype dosage is shuffled within population while the relevant family unit is fixed.
- Genetic-distance QC: complete multilocus profiles are shuffled within population while observed relative pairs stay fixed.
- Grouped RF: the null comparison is empirical out-of-family performance against an identical baseline fold structure; paired family MSE differences are assessed by sign-flip testing.

## Model-specific missingness

Rows are preserved in public analysis-ready inputs. Downstream scripts construct each model-specific sample from the variables required for that model. Age-adjusted adult analyses use behavior-wave age and do not substitute a separate genotype-file age for missing behavioral-wave age.

## Privacy

Only de-identified public inputs are supported. Reconstruction from raw names, free-text field identifiers, or private internal linkage tables is outside the public pipeline.
