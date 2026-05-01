# Reproduction and Public Workflow Harnesses

This directory contains runnable old-vs-new comparisons for published GenomicSEM workflows and one
large public-data practical workflow.

## Status at a Glance

| Harness | What it establishes | Main limitation |
|---|---|---|
| `grotzinger_2019_nhb.R` | Exact no-SNP reproduction plus old-vs-new benchmarking of the published four-trait workflow shape | The original full GWAS inputs used by the practical are not publicly downloadable |
| `grotzinger_2025_nature.R` | Same-input old-vs-new benchmarking of the reported 14-disorder workflow shape | Public inputs are rounded/incomplete, so this is not an exact full-paper rerun |
| `pfactor_practical_1m.R` | Real public-data 1M-SNP practical workflow with bounded old-vs-new equivalence checks and a full Rust-backed scan | This is a public tutorial/practical workflow, not a paper reproduction |
| `pfactor_old_lavaan_scaling.R` | Bounded legacy lavaan scaling study used to project the old 1M-SNP runtime | Projection rather than a literal full old 1M run |

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

Recorded result from the current benchmark log:

- exact no-SNP `usermodel()` reproduction matched the practical with max absolute difference
  `2.611839e-06`;
- on the paper-shaped SNP workflow, `userGWAS(Q_SNP=TRUE, sub="F1~SNP")` improved
  `24.259s -> 1.013s` at one core, and `commonfactorGWAS()` improved `53.927s -> 0.056s`.

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

Recorded result from the current benchmark log:

- the public rounded five-factor `usermodel()` check is intentionally not marked equivalent to the
  paper fit because the full unrounded `V` matrix is not public;
- on the paper-shaped 14-disorder workflow, `userGWAS(Q_SNP=TRUE)` improved
  `14.628s -> 1.968s` at one core.

## Public p-Factor Practical, O(1M) SNPs

`pfactor_practical_1m.R` reruns the public GenomicSEM p-factor practical shape on real public
summary statistics for SCZ, BIP, and MDD:

- downloads the public PGC summary-stat files referenced by the practical workflow;
- downloads `GenomicSEMPractical.RData` from the public UT Box practical materials and uses
  `PSYCH_COV`;
- builds a one-million-SNP aligned subset from the real public files using BIP allele frequencies as
  the `sumstats()` reference;
- benchmarks old R `sumstats()` vs the Rust prep engine;
- benchmarks the constrained practical `userGWAS()` model old vs new, and runs the Rust path on the
  full one-million-SNP table.

The practical model is:

```r
F1 =~ SCZ + BIP + MDD
F1 ~ SNP
SCZ ~~ a*SCZ
BIP ~~ b*BIP
MDD ~~ c*MDD
a > .001
b > .001
c > .001
```

`Q_SNP` defaults to `FALSE` in this harness so the benchmark measures the main constrained factor
effect scan used for the 1M-SNP run. Pass `--q-snp` if you explicitly want the Q-enabled variant.

`--old-gwas-snps` controls the bounded old-vs-new comparison window. Keep it small for the normal
1M benchmark, then use `pfactor_old_lavaan_scaling.R` for a defensible legacy full-runtime estimate.
Setting `--old-gwas-snps` equal to `--target-snps` requests the literal full slow old baseline.

Run a small smoke test:

```sh
Rscript repro/pfactor_practical_1m.R --target-snps 1000 --old-gwas-snps 50 --cores 1
```

Run the intended 1M benchmark with a bounded old-path equivalence window:

```sh
Rscript repro/pfactor_practical_1m.R --target-snps 1000000 --old-gwas-snps 100 --cores 1,4,16 --threads 16
```

Outputs are written to `repro/results/pfactor_practical_1m_<timestamp>.csv`.

Recorded remote result from the current benchmark log:

- `sumstats()` on 1,291,369 aligned public SNPs improved `40.821s -> 9.816s`;
- the bounded 100-SNP old-vs-new `userGWAS()` check matched with max absolute numeric difference
  `1.8e-06`;
- the full Rust-backed 1M-SNP scan completed in `12.380s` at 16 threads.

To estimate the full old-lavaan 1M runtime without committing to the full slow run,
`pfactor_old_lavaan_scaling.R` runs only the old `userGWAS()` lavaan path across bounded SNP
counts and fits a linear elapsed-time slope:

```sh
Rscript repro/pfactor_old_lavaan_scaling.R \
  --sizes 1000,5000,10000,25000,50000 \
  --cores 16 \
  --new-full-sec 12.380
```

This script expects `GenomicSEMPractical.RData` and `sumstats_new_1000000.rds` under
`repro/data/pfactor_practical/`. Outputs are written to
`repro/results/pfactor_old_lavaan_scaling_<timestamp>.csv` and
`repro/results/pfactor_old_lavaan_scaling_projection_<timestamp>.csv`.

Recorded remote result from the current benchmark log:

- old lavaan stabilized near `12.2 ms/SNP` by 25k-50k SNPs on the 16-CPU panda pod;
- the fitted old-path projection was `3.363 h` for 1M SNPs (`R^2 = 0.999996`);
- compared with the measured Rust-backed 1M scan at `12.380 s`, the projected model-fitting speedup
  was `978.1x`.
