# Supplementary package map

The manuscript uses three public Online Resources.

## Online Resource 1 — Supplementary Materials and Figures

File:
- `Online_Resource_1_Supplementary_Materials_and_Figures_Hadza_Datoga.pdf`

Purpose:
- overview of the public reproducibility package;
- technical workflow figure;
- strict adult family-graph component summary;
- kinship-scale decomposition framework;
- focal kinship-scale effect summaries.

Main code dependencies:
- `scripts/02_family_graph_reconstruction.py`
- `scripts/05_family_graph_permutations.py`
- `scripts/06_kinship_scale_decomposition.py`
- `scripts/10_generate_tables_figures.py`

## Online Resource 2 — Supplementary Tables

File:
- `Online_Resource_2_Supplementary_Tables_Hadza_Datoga.xlsx`

Contents and code map:
- S1 / sample structure → `01_data_harmonization.py`, `03_age_genotype_models.py`
- S2 / age-stage and cross-age genotype models → `03_age_genotype_models.py`
- S3 / developmental and vertical parent–offspring analyses → `04_parent_offspring_permutations.py`
- S4 / full family graph and genetic QC → `02_family_graph_reconstruction.py`, `05_family_graph_permutations.py`, `07_multilocus_and_dispersion.py`
- S5 / kinship-scale decomposition → `06_kinship_scale_decomposition.py`
- S6 / multilocus and phenotype dispersion → `07_multilocus_and_dispersion.py`
- S7 / LASSO → `08_lasso_analysis.py`
- S8 / grouped Random Forest → `09_grouped_random_forest.py`
- S9 / reproducibility QC → all scripts 01–10 and `docs/reproducibility_check.md`

## Online Resource 3 — Supplementary Data

File:
- `Online_Resource_3_Supplementary_Data_Hadza_Datoga.xlsx`

Contents:
- S1_Demography — aggregate sample composition, age completeness, and age-group counts;
- S2_AgeBehaviorSummary — aggregate age-group means/SD for physical and verbal aggression, with suppression of behavioral summaries for cells with N<10;
- S3_AnalysisSampleN — model-specific sample sizes used in the manuscript;
- S4_FamilySummary — aggregate family-graph counts only;
- Coding_QC — genotype coding and panel-level QC definitions without participant genotypes.

Participant-level behavioral, demographic, genotype, node/edge, and parent-line records are not publicly distributed. The analysis scripts remain public, but full reruns require the restricted internal analysis-ready inputs governed by the original ethics/data-governance conditions.

## Final assembly

`scripts/10_generate_tables_figures.py` reads outputs of modules 01–09 and assembles manuscript-level tables, supplementary registries, provenance summaries, and selected figures.

## Privacy

The public Supplementary Data package is aggregate-only. It contains no participant-level behavioral values, individual genotypes, node/edge identifiers, parent–offspring identifiers, original field identifiers, or free-text kin labels.
