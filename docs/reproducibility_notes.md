# Reproducibility notes

This repository was validated against the frozen Hadza–Datoga working analysis before the v1.0.0/Zenodo release.

## Exact deterministic reconstruction

The internally validated analysis-ready package reproduces the principal structural and deterministic quantities, including:

- child phenotype sample counts (Hadza/Datoga);
- core9, sensitivity7 and Hadza-only 11-locus complete-case panel sizes;
- the strict adult family graph (109 nodes, 178 relations, 30 connected components);
- all 16 observed vertical mother/father phenotype Spearman correlations;
- the strict HTR2A × child/adult-stage physical-aggression coefficient and partial F statistic;
- the key kinship-scale coefficients for DAT1, HTR1A, HTR2A and SERTPr.

One archived D07 directed parent–offspring relation is retained in the restricted internal analysis-ready data as an explicit directional exception because this is required to reproduce the frozen vertical parent–offspring analysis. No participant-level edge registry is included in the public OR3.

## Stochastic reconstruction

Historical random-number-generator states and exact CV split objects were not preserved for all analyses. The public code therefore uses fixed release seeds when supplied with the restricted analysis-ready inputs.

Consequently, Monte Carlo p-values, LASSO selection frequencies and some Random Forest quantities may differ slightly from the historical frozen run while preserving the same sample definitions, deterministic estimates, effect direction and substantive conclusion.

Examples include:
- the strict HTR2A × stage permutation p-value;
- family-component permutation p-values;
- LASSO stability frequencies / out-of-fold estimates;
- RF paired sign-flip p-values.

The frozen manuscript/Supplementary tables remain the archival record of the historical stochastic run.

## Multilocus block

The validated internal reconstruction reproduces the fixed core9/sensitivity7 sample sizes and very closely reproduces the joint-block partial ΔR² values. Small numerical differences remain because the exact historical implementation details were not fully preserved. These differences do not change the null omnibus conclusion.

## Genetic-distance QC

The frozen family-graph workbook describes the validation metric as an allele-mismatch distance averaged across nine core loci using exact valid diploid calls. The precise historical parser/allele-matching implementation for the multiallelic/VNTR markers was not preserved.

The public code documents a reconstructed QC implementation. With the restricted analysis-ready inputs it reproduces the exact informative-pair count (78; 60 Hadza, 18 Datoga), while the historical distance ratio (0.84929) is treated as an archival frozen value rather than tuned by code.

This limitation concerns a secondary graph-validation QC statistic and does not affect the manuscript's primary phenotype, kinship-scale, LASSO or grouped-RF conclusions.

## Release rule

No deterministic or sample-structure mismatch is accepted for v1.0.0. Stochastic differences are accepted only when the historical RNG/split state was not preserved and the public fixed-seed pipeline retains the same effect direction, scale and inferential conclusion.


## Public data boundary

The current public Online Resource 3 is aggregate-only. It contains demographic and age-group summaries, model-specific sample sizes, aggregate family-graph structure, and coding/QC definitions, but no participant-level phenotype, genotype, family-node/edge, or parent-line records. Full reruns of participant-level analyses therefore require restricted internal inputs and are subject to the original ethics and data-governance conditions.
