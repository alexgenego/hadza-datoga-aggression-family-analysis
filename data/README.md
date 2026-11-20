# Data

This directory documents the expected structure of the restricted analysis-ready inputs used by the public scripts.

Participant-level behavioral, demographic, genotype, family-node/edge, and parent–offspring data are **not** publicly released. The underlying field datasets remain active research resources and are subject to ethics/data-governance restrictions; in addition, the combination of small geographically defined populations, family structure, and phenotype/genotype information creates a nontrivial risk of deductive re-identification.

The restricted internal inputs used in the analyses include the following logical tables:

1. child phenotype/demography analysis-ready table;
2. adult phenotype/demography analysis-ready table, retaining rows with missing age where appropriate;
3. family-graph node and edge registries, including relation type and directed parent→offspring structure;
4. locus-specific genotype matrices for shared child/adult markers and the expanded adult panel.

Model-specific samples are generated programmatically by the scripts and are not created by manually deleting rows from the source matrices.

The public `supplementary/Online_Resource_3_Supplementary_Data_Hadza_Datoga.xlsx` contains aggregate demographic/age-group summaries, model-specific sample sizes, aggregate family-graph counts, and coding/QC definitions only.
