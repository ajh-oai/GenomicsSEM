# Grotzinger 2019 Nature Human Behaviour Reproduction Harness

This directory contains a runnable old-vs-new comparison for the GenomicSEM workflow introduced by Grotzinger et al. 2019, *Nature Human Behaviour*.

The script uses the published psychiatric common-factor example from the GenomicSEM wiki / CNS Genomics practical:

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
