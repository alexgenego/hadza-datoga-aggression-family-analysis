# Hadza–Datoga aggression family analysis

Reproducible analysis code for the manuscript on age-, family-, population-, and genotype-related components of aggressive behavior in Hadza and Datoga samples.

## Scope

This repository contains the final analysis pipeline used to reproduce the manuscript-level results. The pipeline covers:

1. data harmonization and construction of analysis-ready tables;
2. reconstruction of the family graph from de-identified relational tables;
3. age × genotype analyses for markers shared between child and adult datasets;
4. directed parent–offspring permutation tests;
5. full-family-graph permutation tests;
6. kinship-scale decomposition of locus-specific family associations;
7. multilocus, dispersion, and genetic-distance QC analyses;
8. LASSO with nested cross-validation and stability selection;
9. grouped Random Forest with leave-one-family-component-out validation;
10. generation of manuscript tables, figures, and Supplementary registries.

## Central analysis principle

The main design principle is that genetic associations may emerge at different scales of family organization. The analysis therefore distinguishes:

- horizontal within-generation sibling structure;
- mother- and father-defined offspring lines;
- direct parent-genotype → offspring-phenotype tests;
- mixed-generation family components;
- the full connected pedigree graph.

Changes in effect size across these levels are treated as characteristics of the scale of family organization rather than as a simple sequence of replication versus non-replication.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── data/
│   └── README.md
├── docs/
│   └── analysis_workflow.md
├── outputs/
│   └── README.md
└── scripts/
    ├── 01_data_harmonization.py
    ├── 02_family_graph_reconstruction.py
    ├── 03_age_genotype_models.py
    ├── 04_parent_offspring_permutations.py
    ├── 05_family_graph_permutations.py
    ├── 06_kinship_scale_decomposition.py
    ├── 07_multilocus_and_dispersion.py
    ├── 08_lasso_analysis.py
    ├── 09_grouped_random_forest.py
    ├── 10_generate_tables_figures.py
    └── README.md
```

## Reproducibility and privacy

Raw field archives and participant-level analytical tables are **not** included in this public repository. This includes individual behavioral values, individual genotypes, node/edge registries, parent–offspring identifiers, names, free-text identifiers, and private linkage tables.

The analysis scripts remain public and document the full analytical pipeline. Model-specific samples were generated programmatically from restricted internal analysis-ready inputs; participant-level inputs are not publicly distributed because the field datasets remain active research resources and because the combination of small geographically defined populations, family structure, and phenotype/genotype information creates a nontrivial deductive re-identification risk.

The family pipeline preserves two analytical layers:

- a strict social/field graph used for phenotype analyses;
- a Mendelian-compatible genetic layer in which incompatible parent-offspring links are excluded only from biological-parentage analyses while the social relation remains available for phenotype analyses.

## Environment

The release candidate was checked under Python 3.13.5 with the package versions frozen in `requirements.txt`.

Install dependencies with:

```bash
python -m pip install -r requirements.txt
```

## Execution order

The intended execution order is 01 → 10. Exact input filenames will be fixed when the final Supplementary Data package is assembled.

Typical commands are:

```bash
python scripts/01_data_harmonization.py --input <input.xlsx> --dataset-label adults
python scripts/02_family_graph_reconstruction.py --edges <deidentified_edges.xlsx>
python scripts/03_age_genotype_models.py --data <harmonized_child_adult.csv>
python scripts/04_parent_offspring_permutations.py --links <parent_offspring_lines.csv> --nodes <adult_nodes.csv>
python scripts/05_family_graph_permutations.py --nodes <adult_analysis_with_graph.csv>
python scripts/06_kinship_scale_decomposition.py --nodes <adult_nodes.csv> --sibling-edges <sibling_edges.csv> --mother-lines <mother_defined_offspring_lines.csv> --father-lines <father_defined_offspring_lines.csv>
python scripts/07_multilocus_and_dispersion.py --data <adult_genotype_phenotype.csv> --pairs <family_relation_pairs.csv>
python scripts/08_lasso_analysis.py --data <adult_genotype_phenotype.csv> --panel core9
python scripts/08_lasso_analysis.py --data <adult_genotype_phenotype.csv> --panel sensitivity7
python scripts/08_lasso_analysis.py --data <adult_genotype_phenotype.csv> --panel hadza11
python scripts/09_grouped_random_forest.py --data <adult_analysis_with_graph.csv>
python scripts/10_generate_tables_figures.py --outputs-root outputs
```

See `scripts/README.md`, `docs/analysis_workflow.md`, and `docs/reproducibility_notes.md` for the module-to-output map, methodological notes, and the frozen-vs-reconstructed validation record.

## Important implementation details

- Full-graph phenotype permutations keep graph topology fixed and shuffle adjusted residuals within population.
- Parent–offspring permutations preserve parent role, directed parent→offspring structure, offspring counts, and optionally offspring sex.
- Family-component dispersion uses equal component weight.
- LASSO performs covariate handling inside training folds; the 500 stability-selection subsamples preserve population × sex strata.
- Grouped Random Forest uses leave-one-family-component-out validation; relatives never appear simultaneously in training and test folds.
- RF genotype imputation uses training-family medians only and does not add missingness indicators.
- Genetic-distance QC permutes complete genotype profiles within population while family pairs remain fixed.

## Outputs

The analysis modules write machine-readable CSV outputs under `outputs/`. The final script assembles:

- manuscript Table 1 (sample flow / model-specific N);
- manuscript Table 4 (LASSO summary);
- manuscript Table 5 (full-family-graph + grouped RF summary);
- manuscript Table 6 (kinship-scale profiles);
- Supplementary kinship and RF registries;
- a provenance manifest;
- Figure 6 and technical family-component summaries where source outputs are available.

## Data availability

An aggregate-only Supplementary Data workbook is included in `supplementary/Online_Resource_3_Supplementary_Data_Hadza_Datoga.xlsx`, together with Online Resources 1 and 2. It contains demographic summaries, age-group phenotype summaries, model-specific sample sizes, aggregate family-graph structure, and genotype-coding/QC definitions. Participant-level source and analysis-ready data are not publicly released.

## Release status

**Privacy-safe release candidate.** This repository was rebuilt from a clean snapshot after the public data boundary was tightened. OR3 is aggregate-only; participant-level behavioral, genotype, family-node/edge, and parent-line data are not present in this repository or its public history. Scripts 01–10 and the reproducibility documentation remain public.

## Authors

A.M. Kulikov; L.A. Revyakina; P.A. Proshakov; M.L. Butovskaya; O.E. Lazebny.

Corresponding author in the associated Hadza–Datoga materials: O.E. Lazebny.

## License

The analysis code in this repository is released under the MIT License.

This software license does **not** grant rights to restricted participant-level data, identifiable source records, or third-party datasets.


## Archival citation

Privacy-safe GitHub release: `v1.0.2`.

Zenodo DOI for the privacy-safe `v1.0.2` archive: **10.5281/zenodo.24407283**. The earlier `v1.0.1` record remains a restricted historical version. Concept DOI for all versions: **10.5281/zenodo.24247325**.
