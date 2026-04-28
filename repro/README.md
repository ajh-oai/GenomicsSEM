# Grotzinger Reproduction Harnesses

This directory contains runnable old-vs-new comparisons for published GenomicSEM workflows.

## Grotzinger et al. 2019, Nature Human Behaviour

`grotzinger_2019_nhb.R` uses the published psychiatric common-factor example from the GenomicSEM wiki / CNS Genomics practical:

- exact published LDSC `S`, `V`, `I`, `N`, and `m` values for ALCH, PTSD, MDD, and ANX;
- the published `usermodel()` common-factor model and expected fit/loadings;
- a paper-shaped four-trait SNP-effect table for `userGWAS(Q_SNP=TRUE)` and `commonfactorGWAS()`, including the first five chromosome-4 SNP IDs shown in the practical.

The original full GWAS input files from the CNS practical are referenced as cluster-local paths and are not web-downloadable from the public HTML page. For that reason, the script separates:

- **exact reproduction**: the no-SNP common-factor model from the published LDSC matrices;
- **performance equivalence**: same-input old R/lavaan loops vs Rust-backed R bindings on the SNP-effect workflow shape used by the paper/practical;
- **prep-engine timing**: generated four-trait raw summary-stat files with the same columns used by `munge()` and `sumstats()`.

Run:

```sh
Rscript repro/grotzinger_2019_nhb.R --model-snps 1000 --prep-snps 50000 --cores 1,4
```

Useful smaller smoke test:

```sh
Rscript repro/grotzinger_2019_nhb.R --model-snps 50 --prep-snps 1000 --cores 1
```

Outputs are written to `repro/results/grotzinger_2019_nhb_<timestamp>.csv`.

## Grotzinger et al. 2025, Nature

`grotzinger_2025_nature.R` targets the 14-disorder workflow from "Mapping the genetic landscape across 14 psychiatric disorders":

- the 14 LDSC genetic covariance estimates, standard errors, and intercepts reported in Supplementary Table 1;
- the five-factor model described in the paper and Supplementary Table 3;
- the paper's `munge() -> ldsc() -> sumstats() -> userGWAS(Q_SNP=TRUE)` execution shape;
- a 14-trait paper-shaped SNP-effect benchmark for the old R/lavaan path versus the Rust-backed R bindings.

The Nature supplement reports rounded LDSC point estimates and standard errors, but not the full GenomicSEM sampling covariance matrix `V`. The script therefore uses the public rounded `S`, a diagonal `V` from the reported LDSC standard errors, and `nearPD` smoothing when needed. This is sufficient for same-input old-vs-new performance and numerical equivalence checks, but it is not an exact full-paper re-run of the original LDSC/model fit.

Run:

```sh
Rscript repro/grotzinger_2025_nature.R --model-snps 100 --prep-snps 50000 --cores 1,4
```

Useful smaller smoke test:

```sh
Rscript repro/grotzinger_2025_nature.R --model-snps 5 --prep-snps 1000 --cores 1
```

Outputs are written to `repro/results/grotzinger_2025_nature_<timestamp>.csv`.
