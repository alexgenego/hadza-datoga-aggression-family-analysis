# Reproducibility check before v1.0.0

Status: **PASS with documented stochastic/reconstruction caveats**

The participant-level internal analysis-ready data were checked against the frozen manuscript/workbook targets before archiving. The public OR3 has subsequently been replaced by an aggregate-only workbook for ethics/data-governance reasons; the checks below document the internal validation record and are not a claim that participant-level inputs are publicly available.

## Exact structural/sample checks

- Child Hadza/Datoga rows: 490 (223 Hadza, 267 Datoga)
- Physical aggression nonmissing: Hadza 214; Datoga 252
- Verbal aggression nonmissing: Hadza 215; Datoga 256
- Core9 complete-case panel: N=300 (293 Hadza, 7 Datoga)
- Sensitivity7 panel: N=325 (307 Hadza, 18 Datoga)
- Hadza-only 11-locus panel: N=290
- Strict family graph: 109 nodes; 178 relations; 30 connected components
- Vertical parent-offspring observed correlations: all 16 mother/father × trait × model Spearman rho values reproduced to numerical precision

## Key deterministic/near-deterministic effects

- HTR2A × developmental/protocol stage, physical aggression: beta=-0.5313146; F=6.02238; N=381 (66 children + 315 adults), exact
- DAT1 verbal, mother-defined offspring lines: rho=-0.8731696; 10 lines / 19 offspring, exact
- HTR1A physical aggression, full graph: rho=+0.4999136, exact
- HTR2A hostility, mother-defined lines: rho=-0.5765154, exact
- SERTPr hostility, mother/father lines: +0.6785758 / -0.6390097, exact
- Core9 physical partial delta-R2: reconstructed ≈0.03389 vs frozen 0.0342517
- Core9 verbal partial delta-R2: reconstructed ≈0.03478 vs frozen 0.0348180

## Stochastic checks

Historical RNG/split objects were not preserved for all analyses. The public code therefore uses fixed release seeds when run with the restricted analysis-ready inputs. Small Monte Carlo/CV differences are acceptable only when direction, scale, sample structure, and inference are preserved.

- HTR2A × stage permutation p: frozen 0.0124; public fixed-seed estimate ≈0.0135
- Pooled max-sample Hostility family clustering: frozen ratio=0.583809, p=0.00997; reconstructed ratio≈0.585, p≈0.010
- Hadza age-adjusted Hostility clustering: frozen ratio=0.484696, p=0.00760; reconstructed ratio≈0.485, p≈0.007–0.008
- LASSO: fixed cohort sizes reproduce exactly; qualitative conclusion remains no stable predictor >=0.70, lambda1se sparse/null, OOF gain near zero
- Grouped RF max-sample Hostility delta-R2: ≈-0.08211
- Grouped RF age-known physical aggression delta-R2: ≈+0.14779; direction persists for min_samples_leaf 2/4/6

## Documented reconstruction limit

The archived Core9 genetic-distance validation reports 78 pairs, ratio=0.84929, p=0.0033. The 78-pair structure is reproduced exactly, but the precise historical multiallelic/VNTR allele-mismatch parser was not preserved. The public code is **not** tuned to force the frozen 0.849 ratio.

## Privacy

The current public Supplementary Data workbook is aggregate-only and contains no participant-level behavioral values, individual genotypes, node/edge identifiers, or directed parent–offspring identifiers. Earlier de-identification checks remain part of the internal validation history, but de-identification alone is not treated as sufficient for public release because deductive re-identification risk remains in these small, geographically defined family samples.

## Release rule

Deterministic or sample-structure mismatches were resolved before the archived analysis was frozen. Stochastic differences are acceptable only under the conditions above. The public repository retains the code and validation record, whereas participant-level inputs remain restricted.
