# Benchmark Log

Append new entries chronologically. Include the timestamp, code change, command, hardware/context, and the result table or summary.

## 2026-04-23 19:20 PDT - Rust kernel and R binding baseline

Change set:

- Added Rust kernels for `.get_V_SNP`, `.get_V_full`, `.get_S_Full`, and `.get_Z_pre`.
- Added `.Call` bindings and R wrappers preserving old R helpers as `*_r`.
- Added `.get_V_SNP_batch` with Rayon.

Local command:

```sh
Rscript benches/compare_backends.R 200000 8 standard 1,2,4,8,16
```

Local result summary:

| backend | threads | elapsed_sec | note |
|---|---:|---:|---|
| old_r_loop | 1 | 11.233 | original R helper in SNP loop |
| r_binding_loop | 1 | 1.255 | R loop calling Rust per SNP |
| r_binding_batch | 1 | 0.079 | one R call into Rust batch |
| r_binding_batch | 16 | 0.007 | Rayon batch |
| rust_loop | 1 | 0.021794 | pure Rust CLI, separate deterministic data |
| rust_batch | 16 | 0.003497 | pure Rust CLI, separate deterministic data |

Remote context:

- brix pod: 16 CPU, panda cluster, flex quota.
- R installed on pod via apt.

Remote command:

```sh
Rscript benches/compare_backends.R 200000 8 standard 1,2,4,8,16
```

Remote result summary:

| backend | threads | elapsed_sec |
|---|---:|---:|
| old_r_loop | 1 | 28.030 |
| r_binding_loop | 1 | 1.847 |
| r_binding_batch | 1 | 0.151 |
| r_binding_batch | 8 | 0.033 |
| r_binding_batch | 16 | 0.043 |
| rust_loop | 1 | 0.071217 |
| rust_batch | 8 | 0.028058 |
| rust_batch | 16 | 0.022602 |

Interpretation:

- The Rust kernel is much faster than the old R matrix fill.
- The full workflow still needs workflow-level work because lavaan/model overhead dominates once the kernel is embedded per SNP.

## 2026-04-23 20:04 PDT - Full workflow synthetic benchmark added

Change set:

- Added `benches/synthetic_inputs.R`.
- Added `benches/bench_usergwas_synthetic.R`.
- The fixture follows the documented workflow shape: LDSC-like `covstruc`, sumstats-like SNP table, and lavaan model syntax.

Command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 4 1,2,4 userGWAS 1
Rscript benches/bench_usergwas_synthetic.R 100 4 1,2,4 commonfactorGWAS 1
```

Results:

| workflow | backend | cores | elapsed_sec | checksum |
|---|---|---:|---:|---:|
| userGWAS | old_r_workflow | 1 | 2.582 | 224714.6 |
| userGWAS | old_r_workflow | 4 | 0.882 | 224714.6 |
| userGWAS | rust_binding_workflow | 1 | 2.541 | 224714.6 |
| userGWAS | rust_binding_workflow | 4 | 0.914 | 224714.6 |
| commonfactorGWAS | old_r_workflow | 1 | 6.118 | 10732.62 |
| commonfactorGWAS | old_r_workflow | 4 | 1.743 | 10732.62 |
| commonfactorGWAS | rust_binding_workflow | 1 | 5.739 | 10732.62 |
| commonfactorGWAS | rust_binding_workflow | 4 | 1.712 | 10732.62 |

Interpretation:

- Checksums match for old R and Rust-backed paths.
- The Rust helper is not the dominant cost in full workflows at this scale.

## 2026-04-23 20:33 PDT - Diagonal inverse optimization

Change set:

- Replaced repeated `solve(diag(values))` with `.diag_inverse_from_values()`.
- The helper falls back to `solve()` for non-finite or ill-conditioned diagonals.
- Added `options(GenomicSEM.fast_diag_inverse)` for A/B checks; default is `TRUE`.

Command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 FALSE,TRUE
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 commonfactorGWAS 1 FALSE,TRUE
```

Results:

| workflow | backend | cores | fast_diag_inverse | elapsed_sec | checksum |
|---|---|---:|---|---:|---:|
| userGWAS | old_r_workflow | 1 | FALSE | 10.735 | 823014.4 |
| userGWAS | old_r_workflow | 1 | TRUE | 11.472 | 823014.4 |
| userGWAS | rust_binding_workflow | 1 | FALSE | 10.656 | 823014.4 |
| userGWAS | rust_binding_workflow | 1 | TRUE | 10.774 | 823014.4 |
| commonfactorGWAS | old_r_workflow | 1 | FALSE | 30.039 | 11616.85 |
| commonfactorGWAS | old_r_workflow | 1 | TRUE | 28.263 | 11616.85 |
| commonfactorGWAS | rust_binding_workflow | 1 | FALSE | 30.022 | 11616.85 |
| commonfactorGWAS | rust_binding_workflow | 1 | TRUE | 27.917 | 11616.85 |
| commonfactorGWAS | old_r_workflow | 4 | FALSE | 7.333 | 11616.85 |
| commonfactorGWAS | old_r_workflow | 4 | TRUE | 7.164 | 11616.85 |

Interpretation:

- Useful for `commonfactorGWAS` at 12 traits: about 6-7% faster single-core in this synthetic run.
- Not a reliable `userGWAS` win; lavaan and output assembly dominate.

## 2026-04-23 20:45 PDT - DWLS diagonal WLS row-scaling experiment

Change set:

- Added opt-in `options(GenomicSEM.fast_diagonal_wls)` around replacing `W %*% delta` with row scaling by `diag(W)`.
- Default is `FALSE` because the workflow benchmark did not show a stable gain.

Command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 TRUE FALSE,TRUE
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 commonfactorGWAS 1 TRUE FALSE,TRUE
```

Results:

| workflow | backend | cores | fast_diagonal_wls | elapsed_sec | checksum |
|---|---|---:|---|---:|---:|
| userGWAS | old_r_workflow | 1 | FALSE | 11.474 | 823014.4 |
| userGWAS | old_r_workflow | 1 | TRUE | 12.075 | 823014.4 |
| userGWAS | rust_binding_workflow | 1 | FALSE | 11.456 | 823014.4 |
| userGWAS | rust_binding_workflow | 1 | TRUE | 11.277 | 823014.4 |
| commonfactorGWAS | old_r_workflow | 1 | FALSE | 32.310 | 11616.85 |
| commonfactorGWAS | old_r_workflow | 1 | TRUE | 32.379 | 11616.85 |
| commonfactorGWAS | rust_binding_workflow | 1 | FALSE | 32.512 | 11616.85 |
| commonfactorGWAS | rust_binding_workflow | 1 | TRUE | 32.244 | 11616.85 |

Interpretation:

- Checksums match, but runtime impact is noise-level.
- Keep as opt-in until larger remote benchmarks justify enabling it.

## 2026-04-23 20:55 PDT - `userGWAS(printwarn=FALSE)` output-size fix

Change set:

- Fixed `.userGWAS_main()` result naming when `printwarn=FALSE`.
- Added `result_cols` and `result_size_bytes` to `bench_usergwas_synthetic.R`.

Command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 TRUE FALSE FALSE,TRUE
```

Results:

| backend | cores | printwarn | elapsed_sec | checksum | result_rows | result_cols | result_size_bytes |
|---|---:|---|---:|---:|---:|---:|---:|
| old_r_workflow | 1 | FALSE | 10.997 | 823014.4 | 3900 | 19 | 1036048 |
| old_r_workflow | 1 | TRUE | 11.333 | 823014.4 | 3900 | 21 | 1140048 |
| rust_binding_workflow | 1 | FALSE | 11.163 | 823014.4 | 3900 | 19 | 1036048 |
| rust_binding_workflow | 1 | TRUE | 10.843 | 823014.4 | 3900 | 21 | 1140048 |

Interpretation:

- Fixes a broken documented memory/output-size control.
- On this small synthetic case, object size drops about 9%; large real runs with distinct warning strings may benefit more.
- Runtime is noise-level; the primary value is output memory pressure.

## 2026-04-23 22:39 PDT - Workflow profiling and high-impact optimization search

Change set:

- Added `benches/profile_workflows.R` for small Rprof profiles.
- Added `benches/profile_commonfactor_phases.R` for manual phase profiling when Rprof is too invasive or unstable.
- Cached the native-symbol availability check in `.genomicssem_use_rust()`; this removes repeated `is.loaded()` checks from hot R loops but is not a headline win.

Commands:

```sh
Rscript benches/profile_workflows.R 100 3 userGWAS 1
Rscript benches/profile_commonfactor_phases.R 10 12 1
```

`userGWAS()` Rprof result, 100 SNPs, 3 traits:

| stack/function | total_pct | self_pct |
|---|---:|---:|
| `.userGWAS_main` | 92.47 | 0.22 |
| `lavaan` | 80.21 | 0.03 |
| `lav_lavaan_step11_estoptim` | 35.65 | 0.05 |
| `nlminb` | 31.79 | 0.05 |
| `solve.default` | 11.45 | 7.26 |

`commonfactorGWAS()` manual phase split, 10 SNPs, 12 traits:

| phase | elapsed_sec | pct |
|---|---:|---:|
| main_lavaan | 5.878 | 57.97 |
| q_lavaan | 3.982 | 39.27 |
| main_sandwich_and_q_setup | 0.131 | 1.29 |
| q_sandwich | 0.127 | 1.25 |
| build_v | 0.017 | 0.17 |
| build_s | 0.005 | 0.05 |

Interpretation:

- Workflow time is dominated by lavaan optimization, not the Rust-backed matrix kernels.
- For `commonfactorGWAS()`, the follow-up Q model is about 39% of the profiled runtime. That is the only obviously high-impact removable cost without replacing lavaan's optimizer or changing the SEM algorithm.
- Rprof was unstable for local `commonfactorGWAS()` profiles at larger trait counts on this macOS/R build, so the phase profiler is the safer k=12 tool.

Rejected experiments:

- `se="none"` inside lavaan matched numeric outputs in paired tests, but it was slower in the end-to-end benchmark and was backed out.
- A userGWAS-style analytic Q shortcut for `commonfactorGWAS()` was close but not exact on an 8-SNP/6-trait check; observed Q differences ranged roughly from `-7.7e-4` to `3.2e-4`, so it was rejected for exact-behavior preservation.

## 2026-04-23 22:46 PDT - `commonfactorGWAS(Q_SNP=FALSE)`

Change set:

- Added `commonfactorGWAS(..., Q_SNP=TRUE)`, defaulting to the original behavior.
- `Q_SNP=FALSE` skips the second per-SNP lavaan Q model and returns missing Q statistics.
- Added `tests/commonfactor-qsnp.R`.

100-SNP/12-trait command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 commonfactorGWAS 1 TRUE FALSE TRUE TRUE,FALSE
```

100-SNP/12-trait results:

| backend | cores | Q_SNP | elapsed_sec | checksum |
|---|---:|---|---:|---:|
| old_r_workflow | 1 | TRUE | 30.967 | 11616.85 |
| old_r_workflow | 4 | TRUE | 7.480 | 11616.85 |
| rust_binding_workflow | 1 | TRUE | 31.697 | 11616.85 |
| rust_binding_workflow | 4 | TRUE | 7.474 | 11616.85 |
| old_r_workflow | 1 | FALSE | 19.242 | 11409.60 |
| old_r_workflow | 4 | FALSE | 4.660 | 11409.60 |
| rust_binding_workflow | 1 | FALSE | 19.523 | 11409.60 |
| rust_binding_workflow | 4 | FALSE | 4.684 | 11409.60 |

Current 20-SNP/12-trait sanity command:

```sh
Rscript benches/bench_usergwas_synthetic.R 20 12 1 commonfactorGWAS 1 TRUE FALSE TRUE TRUE,FALSE
```

Current 20-SNP/12-trait sanity results:

| backend | Q_SNP | elapsed_sec | checksum |
|---|---|---:|---:|
| old_r_workflow | TRUE | 10.532 | 721.1656 |
| rust_binding_workflow | TRUE | 11.653 | 721.1656 |
| old_r_workflow | FALSE | 6.891 | 683.0168 |
| rust_binding_workflow | FALSE | 6.706 | 683.0168 |

Interpretation:

- This is a high-impact optimization only for users who do not need the SNP heterogeneity Q statistic.
- The speedup matches the profile attribution: about 37-38% on the 100-SNP/12-trait benchmark and about 35-42% on the smaller sanity run.
- The checksum changes when `Q_SNP=FALSE` because Q and Q p-value fields are intentionally missing; `tests/commonfactor-qsnp.R` verifies that SNP effect estimates and corrected SEs are unchanged.

## 2026-04-24 00:10 PDT - Experimental Rust common-factor model fit

Change set:

- Added a specialized Rust Gauss-Newton DWLS solver for the generated one-factor `commonfactorGWAS()` model.
- Added `.commonfactor_fit_fast()` and `options(GenomicSEM.fast_commonfactor_fit=TRUE)`.
- The fast path currently activates only for `commonfactorGWAS(Q_SNP=FALSE)`; `Q_SNP=TRUE` still falls back to lavaan because the Q model is a second SEM fit with different free-parameter structure.
- Added `tests/commonfactor-fast-fit.R` comparing the Rust fit against lavaan on a random non-compound covariance fixture.

Correctness check:

```sh
Rscript tests/commonfactor-fast-fit.R
Rscript - <<'RS'
library(GenomicSEM); source('benches/synthetic_inputs.R')
inputs <- make_synthetic_genomicsem_inputs(n_snp=100,k=12,seed=1)
run <- function(fast) {
  options(GenomicSEM.use_rust=TRUE, GenomicSEM.fast_diag_inverse=TRUE, GenomicSEM.fast_commonfactor_fit=fast)
  suppressWarnings(suppressMessages(capture.output({
    out <- commonfactorGWAS(covstruc=inputs$covstruc,SNPs=inputs$SNPs,parallel=FALSE,GC='standard',Q_SNP=FALSE)
  })))
  out
}
a <- run(FALSE); b <- run(TRUE)
cat(max(abs(a$est-b$est)), max(abs(a$se_c-b$se_c)), max(abs(a$Z_Estimate-b$Z_Estimate)), max(abs(a$Pval_Estimate-b$Pval_Estimate)), "\n")
RS
```

100-SNP/12-trait parity result:

| field | max_abs_diff |
|---|---:|
| est | 8.420874e-08 |
| se_c | 3.423040e-10 |
| Z_Estimate | 5.125942e-06 |
| Pval_Estimate | 4.089636e-06 |

Benchmark command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 commonfactorGWAS 1 TRUE FALSE TRUE FALSE FALSE,TRUE
```

Results:

| backend | cores | Q_SNP | fast_commonfactor_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | FALSE | FALSE | 16.349 | 11409.6 |
| old_r_workflow | 4 | FALSE | FALSE | 4.612 | 11409.6 |
| rust_binding_workflow | 1 | FALSE | FALSE | 15.837 | 11409.6 |
| rust_binding_workflow | 4 | FALSE | FALSE | 4.524 | 11409.6 |
| old_r_workflow | 1 | FALSE | TRUE | 15.991 | 11409.6 |
| old_r_workflow | 4 | FALSE | TRUE | 4.490 | 11409.6 |
| rust_binding_workflow | 1 | FALSE | TRUE | 0.690 | 11409.6 |
| rust_binding_workflow | 4 | FALSE | TRUE | 0.387 | 11409.6 |

Interpretation:

- This is the first model-fitting rewrite that materially changes the full workflow: `15.837s -> 0.690s` single-core for the Rust-backed path on this benchmark.
- It is not byte-identical to lavaan; it is a numerical solver matching lavaan within small tolerances on tested fixtures. Keep it opt-in while broadening parity tests.
- The next high-impact Rust work is a corresponding solver for the Q model, or a generalized parameter-table compiler for the constrained `userGWAS()` lavaan syntax.

## 2026-04-24 00:59 PDT - Rust common-factor Q_SNP fit

Change set:

- Added a specialized Rust Gauss-Newton solver for the `commonfactorGWAS(Q_SNP=TRUE)` direct-effect Q model.
- The fast path now uses Rust for both the main common-factor fit and the Q model when `options(GenomicSEM.fast_commonfactor_fit=TRUE)`.
- Added `.commonfactor_q_fit_fast()` and a native `genomicssem_fit_commonfactor_q_call` binding.
- The Q fit fixes the main-model loadings, factor regression, factor residual variance, and SNP variance, then frees SNP direct effects plus trait residual variances. Its sandwich covariance returns the direct-effect covariance block used by the existing Q statistic.
- The Q solver defaults to 500 iterations via `GenomicSEM.fast_commonfactor_q_max_iter`; if it cannot produce a finite converged Q, the code falls back to the lavaan path for that SNP.

Correctness check:

```sh
Rscript tests/commonfactor-fast-fit.R
Rscript - <<'RS'
library(GenomicSEM); source('benches/synthetic_inputs.R')
inputs <- make_synthetic_genomicsem_inputs(n_snp=100L,k=12L,seed=1L)
run <- function(fast) {
  options(GenomicSEM.use_rust=TRUE, GenomicSEM.fast_diag_inverse=TRUE, GenomicSEM.fast_commonfactor_fit=fast)
  suppressWarnings(commonfactorGWAS(covstruc=inputs$covstruc,SNPs=inputs$SNPs,estimation='DWLS',parallel=FALSE,GC='standard',Q_SNP=TRUE))
}
a <- run(FALSE); b <- run(TRUE)
cat(max(abs(a$est-b$est)), max(abs(a$se_c-b$se_c)), max(abs(a$Z_Estimate-b$Z_Estimate)), max(abs(a$Pval_Estimate-b$Pval_Estimate)), max(abs(a$Q-b$Q)), max(abs(a$Q_pval-b$Q_pval)), "\n")
cat("fast warnings nonzero:", sum(b$warning != 0), "\n")
RS
```

100-SNP/12-trait `Q_SNP=TRUE` parity result:

| field | max_abs_diff |
|---|---:|
| est | 8.420874e-08 |
| se_c | 3.423040e-10 |
| Z_Estimate | 5.125942e-06 |
| Pval_Estimate | 4.089636e-06 |
| Q | 8.043244e-07 |
| Q_pval | 1.057687e-08 |

Fast-path fallback count in this check: 0 nonzero warning rows.

Benchmark command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 commonfactorGWAS 1 TRUE FALSE TRUE TRUE FALSE,TRUE
```

Results:

| backend | cores | Q_SNP | fast_commonfactor_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | TRUE | FALSE | 26.370 | 11616.85 |
| old_r_workflow | 4 | TRUE | FALSE | 7.217 | 11616.85 |
| rust_binding_workflow | 1 | TRUE | FALSE | 27.700 | 11616.85 |
| rust_binding_workflow | 4 | TRUE | FALSE | 7.407 | 11616.85 |
| old_r_workflow | 1 | TRUE | TRUE | 29.354 | 11616.85 |
| old_r_workflow | 4 | TRUE | TRUE | 7.363 | 11616.85 |
| rust_binding_workflow | 1 | TRUE | TRUE | 1.170 | 11616.85 |
| rust_binding_workflow | 4 | TRUE | TRUE | 0.545 | 11616.85 |

Interpretation:

- The Q-enabled workflow now gets the same kind of high-impact improvement as `Q_SNP=FALSE`: `27.700s -> 1.170s` single-core and `7.407s -> 0.545s` on 4 cores for the Rust-backed path.
- This remains an opt-in numerical replacement, not byte-identical lavaan output. The tested synthetic workflow matches the reference to sub-micro Q differences.

## 2026-04-24 01:20 PDT - Remove commonfactorGWAS fast-path startup lavaan

Change set:

- The `commonfactorGWAS()` Rust fast path no longer runs lavaan before the SNP loop.
- Vech ordering for the generated SNP-plus-one-factor model is identity, so the fast path now uses `seq_len((k + 1) * (k + 2) / 2)`.
- Added `.commonfactor_fast_start_from_cov()`, which derives fixed startup values from the first eigenvector of the LDSC genetic covariance matrix.
- Slow/fallback paths still keep the original lavaan setup and lavaan fallback behavior.

Validation:

```sh
R CMD INSTALL --install-tests .
Rscript tests/commonfactor-fast-fit.R
Rscript tests/commonfactor-qsnp.R
```

Benchmark command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 commonfactorGWAS 1 TRUE FALSE TRUE TRUE TRUE
```

Results:

| backend | cores | Q_SNP | fast_commonfactor_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | TRUE | TRUE | 29.070 | 11616.85 |
| old_r_workflow | 4 | TRUE | TRUE | 7.647 | 11616.85 |
| rust_binding_workflow | 1 | TRUE | TRUE | 0.916 | 11616.85 |
| rust_binding_workflow | 4 | TRUE | TRUE | 0.318 | 11616.85 |

Interpretation:

- Removing startup lavaan reduced the Rust-backed `Q_SNP=TRUE` benchmark from `1.170s -> 0.916s` single-core and `0.545s -> 0.318s` on 4 cores.
- This mostly removes fixed overhead, so its relative impact is largest for smaller SNP batches and interactive checks.

## 2026-04-24 01:40 PDT - Experimental Rust userGWAS per-SNP SEM fit

Change set:

- Added a generic Rust RAM-matrix SEM solver using numeric deltas and the same sandwich covariance calculation used by the R/lavaan path.
- Added `.sem_fast_compile()` to compile lavaan's setup parameter table into fixed/free RAM matrices.
- Added `.sem_fit_fast()` and `options(GenomicSEM.fast_usergwas_fit=TRUE)`.
- `userGWAS()` now uses the Rust solver inside the SNP loop for supported DWLS parameter-table models. Unsupported syntax still falls back to lavaan.
- Added `tests/usergwas-fast-fit.R`, covering the supported fixed-measurement one-factor fixture with both `Q_SNP=FALSE` and `Q_SNP=TRUE`.

Validation:

```sh
R CMD INSTALL --install-tests .
Rscript tests/usergwas-fast-fit.R
```

100-SNP/12-trait parity result, `Q_SNP=FALSE`:

| field | max_abs_diff |
|---|---:|
| est | 1.369266e-07 |
| SE | 1.210305e-09 |
| Z_Estimate | 2.737585e-04 |
| Pval_Estimate | 6.608180e-06 |
| chisq | 5.307908e-07 |
| chisq_pval | 1.110223e-16 |
| AIC | 5.307908e-07 |

100-SNP/12-trait parity result, `Q_SNP=TRUE`:

| field | max_abs_diff |
|---|---:|
| est | 1.369266e-07 |
| SE | 1.210305e-09 |
| chisq | 5.307908e-07 |
| Q_SNP | 5.307970e-07 |
| Q_SNP_pval | 1.685171e-09 |

Fast-path fallback count in both checks: 0 nonzero warning rows.

Benchmark command, `Q_SNP=FALSE`:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 TRUE FALSE TRUE FALSE FALSE FALSE,TRUE
```

Results:

| backend | cores | Q_SNP | fast_usergwas_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | FALSE | FALSE | 10.539 | 823014.4 |
| old_r_workflow | 4 | FALSE | FALSE | 2.930 | 823014.4 |
| rust_binding_workflow | 1 | FALSE | FALSE | 10.406 | 823014.4 |
| rust_binding_workflow | 4 | FALSE | FALSE | 2.882 | 823014.4 |
| old_r_workflow | 1 | FALSE | TRUE | 10.537 | 823014.4 |
| old_r_workflow | 4 | FALSE | TRUE | 2.966 | 823014.4 |
| rust_binding_workflow | 1 | FALSE | TRUE | 0.758 | 823014.4 |
| rust_binding_workflow | 4 | FALSE | TRUE | 0.535 | 823014.4 |

Benchmark command, `Q_SNP=TRUE`:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 TRUE FALSE TRUE TRUE FALSE FALSE,TRUE
```

Results:

| backend | cores | Q_SNP | fast_usergwas_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | TRUE | FALSE | 10.582 | 824321.7 |
| old_r_workflow | 4 | TRUE | FALSE | 2.915 | 824321.7 |
| rust_binding_workflow | 1 | TRUE | FALSE | 9.894 | 824321.7 |
| rust_binding_workflow | 4 | TRUE | FALSE | 2.877 | 824321.7 |
| old_r_workflow | 1 | TRUE | TRUE | 11.182 | 824321.7 |
| old_r_workflow | 4 | TRUE | TRUE | 3.091 | 824321.7 |
| rust_binding_workflow | 1 | TRUE | TRUE | 0.789 | 824321.7 |
| rust_binding_workflow | 4 | TRUE | TRUE | 0.503 | 824321.7 |

Interpretation:

- This removes the dominant per-SNP lavaan cost for supported `userGWAS()` DWLS models: `10.406s -> 0.758s` single-core with `Q_SNP=FALSE`, and `9.894s -> 0.789s` with `Q_SNP=TRUE` on this benchmark.
- This checkpoint still used lavaan once up front to build the optimized first-SNP basemodel; the next checkpoint removes that pre-loop fit from the fast path.

## 2026-04-24 01:45 PDT - Remove userGWAS fast-path basemodel lavaan fit

Change set:

- The `userGWAS()` Rust fast path no longer runs the extra first-SNP lavaan fit to obtain a basemodel.
- It compiles the setup `ReorderModel` parameter table directly into the Rust RAM solver.
- `.sem_fit_fast()` now reorders `S_Fullrun` to the compiled observed-variable order before entering Rust.
- Slow/fallback paths still keep the original lavaan behavior.

Validation:

```sh
R CMD INSTALL --install-tests .
Rscript tests/usergwas-fast-fit.R
```

Benchmark command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 TRUE FALSE TRUE TRUE FALSE FALSE,TRUE
```

Results:

| backend | cores | Q_SNP | fast_usergwas_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | TRUE | FALSE | 10.173 | 824321.7 |
| old_r_workflow | 4 | TRUE | FALSE | 2.980 | 824321.7 |
| rust_binding_workflow | 1 | TRUE | FALSE | 9.898 | 824321.7 |
| rust_binding_workflow | 4 | TRUE | FALSE | 2.844 | 824321.7 |
| old_r_workflow | 1 | TRUE | TRUE | 11.317 | 824321.7 |
| old_r_workflow | 4 | TRUE | TRUE | 3.141 | 824321.7 |
| rust_binding_workflow | 1 | TRUE | TRUE | 0.762 | 824321.7 |
| rust_binding_workflow | 4 | TRUE | TRUE | 0.366 | 824321.7 |

Interpretation:

- Removing the pre-loop basemodel lavaan fit mostly helps the parallel benchmark: `0.503s -> 0.366s` on 4 cores for `Q_SNP=TRUE`.
- The setup path still calls lavaan to create the initial parameter table and variable ordering. Removing that final dependency would require a lavaan-syntax parser/parameter-table generator rather than a solver replacement.

## 2026-04-24 02:07 PDT - Remove added commonfactorGWAS Q_SNP flag

Change set:

- Removed the public `commonfactorGWAS(..., Q_SNP=...)` argument that was added during the Rust fast-path work.
- Restored upstream API behavior for `commonfactorGWAS()`: the SNP heterogeneity Q model is always computed.
- Left the existing upstream `userGWAS(..., Q_SNP=...)` flag intact.
- Updated common-factor tests and synthetic benchmark reporting so `Q_SNP` is treated as not applicable for `commonfactorGWAS()`.

Validation:

```sh
R CMD INSTALL --install-tests .
Rscript tests/commonfactor-fast-fit.R
Rscript tests/commonfactor-qsnp.R
Rscript tests/usergwas-fast-fit.R
cargo test --workspace
Rscript tests/rust-kernel-parity.R
Rscript tests/usergwas-printwarn.R
R CMD build .
R CMD check --no-manual GenomicSEM_0.0.5.tar.gz
```

`R CMD check` status: 8 WARNINGs, 4 NOTEs, matching the existing package-level noise class; all package tests passed.

Interpretation:

- This is an API cleanup rather than a new performance checkpoint. The latest relevant common-factor performance numbers remain the Rust fast main+Q path from the previous entries.

## 2026-04-24 02:40 PDT - Batch commonfactorGWAS fast path in Rust

Change set:

- Added a native `genomicssem_fit_commonfactor_batch_call` entry point.
- The supported fast `commonfactorGWAS()` DWLS path now builds SNP-specific `S`/`V`, runs the main common-factor fit, runs the direct-effect Q fit, and computes Q inside one Rust call.
- `parallel=TRUE` now maps to Rust/Rayon threads for this fast path instead of R `foreach`; unsupported cases still fall back to the existing R/lavaan loop.

Validation:

```sh
cargo test --workspace
R CMD INSTALL --install-tests .
Rscript tests/commonfactor-fast-fit.R
Rscript tests/commonfactor-qsnp.R
Rscript tests/rust-kernel-parity.R
```

Benchmark command:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 commonfactorGWAS 1 TRUE FALSE TRUE TRUE FALSE,TRUE FALSE
```

Results:

| backend | cores | fast_commonfactor_fit | elapsed_sec | checksum |
|---|---:|---|---:|---:|
| old_r_workflow | 1 | FALSE | 28.949 | 11616.85 |
| old_r_workflow | 4 | FALSE | 7.467 | 11616.85 |
| rust_binding_workflow | 1 | FALSE | 30.000 | 11616.85 |
| rust_binding_workflow | 4 | FALSE | 7.244 | 11616.85 |
| old_r_workflow | 1 | TRUE | 30.525 | 11616.85 |
| old_r_workflow | 4 | TRUE | 7.509 | 11616.85 |
| rust_binding_workflow | 1 | TRUE | 0.774 | 11616.85 |
| rust_binding_workflow | 4 | TRUE | 0.194 | 11616.85 |

Interpretation:

- Moving the common-factor SNP loop into one native call improves the previous Rust fast main+Q result from `0.916s -> 0.774s` single-core and `0.318s -> 0.194s` on 4 cores for this 100-SNP/12-trait benchmark.
- The old R/lavaan baseline remains around `29s` single-core and `7.5s` on 4 cores, so the batched Rust path is about `37x` faster single-core and `38x` faster on 4 cores here.

## 2026-04-24 02:54 PDT - Batch userGWAS fast model fits in Rust

Change set:

- Added a native `genomicssem_fit_generic_sem_batch_call` entry point for supported RAM-matrix DWLS models.
- The supported fast `userGWAS()` path now builds SNP-specific `S`/`V`, runs the generic SEM fit, computes corrected SEs, computes model chi-square with a block-diagonal `V_full` solve, and computes `Q_SNP` values inside one Rust call.
- R still assembles the lavaan-shaped per-SNP result list, while unsupported cases still fall back to the existing `.userGWAS_main()` path.

Validation:

```sh
cargo test --workspace
R CMD INSTALL --install-tests .
Rscript tests/usergwas-fast-fit.R
Rscript tests/usergwas-printwarn.R
Rscript tests/commonfactor-fast-fit.R
Rscript tests/commonfactor-qsnp.R
R CMD build .
R CMD check --no-manual GenomicSEM_0.0.5.tar.gz
```

`R CMD check` status: 8 WARNINGs, 4 NOTEs, matching the existing package-level noise class; all package tests passed.

Benchmark command, `Q_SNP=FALSE`:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 TRUE FALSE TRUE FALSE FALSE FALSE,TRUE
```

Results:

| backend | cores | Q_SNP | fast_usergwas_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | FALSE | FALSE | 10.494 | 823014.4 |
| old_r_workflow | 4 | FALSE | FALSE | 2.974 | 823014.4 |
| rust_binding_workflow | 1 | FALSE | FALSE | 10.360 | 823014.4 |
| rust_binding_workflow | 4 | FALSE | FALSE | 2.879 | 823014.4 |
| old_r_workflow | 1 | FALSE | TRUE | 12.032 | 823014.4 |
| old_r_workflow | 4 | FALSE | TRUE | 3.119 | 823014.4 |
| rust_binding_workflow | 1 | FALSE | TRUE | 0.479 | 823014.4 |
| rust_binding_workflow | 4 | FALSE | TRUE | 0.242 | 823014.4 |

Benchmark command, `Q_SNP=TRUE`:

```sh
Rscript benches/bench_usergwas_synthetic.R 100 12 1,4 userGWAS 1 TRUE FALSE TRUE TRUE FALSE FALSE,TRUE
```

Results:

| backend | cores | Q_SNP | fast_usergwas_fit | elapsed_sec | checksum |
|---|---:|---|---|---:|---:|
| old_r_workflow | 1 | TRUE | FALSE | 10.562 | 824321.7 |
| old_r_workflow | 4 | TRUE | FALSE | 2.970 | 824321.7 |
| rust_binding_workflow | 1 | TRUE | FALSE | 9.787 | 824321.7 |
| rust_binding_workflow | 4 | TRUE | FALSE | 2.930 | 824321.7 |
| old_r_workflow | 1 | TRUE | TRUE | 11.392 | 824321.7 |
| old_r_workflow | 4 | TRUE | TRUE | 3.100 | 824321.7 |
| rust_binding_workflow | 1 | TRUE | TRUE | 0.446 | 824321.7 |
| rust_binding_workflow | 4 | TRUE | TRUE | 0.256 | 824321.7 |

Interpretation:

- Batching the supported `userGWAS()` Rust fast path improves `Q_SNP=FALSE` from the previous `0.758s -> 0.479s` single-core and `0.535s -> 0.242s` on 4 cores.
- For `Q_SNP=TRUE`, the checkpoint improves from `0.762s -> 0.446s` single-core and `0.366s -> 0.256s` on 4 cores.
- On a traced 100-SNP/12-trait `Q_SNP=TRUE` run, the native batch call took about `0.243s`, R result assembly took about `0.063s`, and total elapsed was about `0.473s`; the rest is mostly one-time lavaan setup used to build the parameter table/order.

## 2026-04-24 03:18 PDT - Release-readiness hardening

Change set:

- Added fast-path attributes on returned objects: `GenomicSEM.fast_path`, `GenomicSEM.fast_threads`, and optional `GenomicSEM.fast_fallback_reason`.
- Added `options(GenomicSEM.fast_diagnostics=TRUE)` for explicit fast-path use/fallback messages.
- Added `options(GenomicSEM.fast_strict=TRUE)` to error instead of silently falling back when a requested batched Rust fast path is unavailable.
- Added `tests/fast-path-release.R` to assert that threaded batch paths are used and strict fallback is enforced.
- Documented the Rust fast-path options, source-build Cargo requirement, and result attributes in README/Rd/PATCHNOTES.
- Fixed the touched `commonfactorGWAS.Rd` and `userGWAS.Rd` usage signatures so they match the current function signatures.

Validation:

```sh
R CMD INSTALL --install-tests .
Rscript tests/fast-path-release.R
Rscript tests/commonfactor-fast-fit.R
Rscript tests/commonfactor-qsnp.R
Rscript tests/usergwas-fast-fit.R
Rscript tests/usergwas-printwarn.R
Rscript tests/rust-kernel-parity.R
cargo test --workspace
R CMD build .
R CMD check --no-manual GenomicSEM_0.0.5.tar.gz
```

`R CMD check` status: 8 WARNINGs, 4 NOTEs. The warnings/notes are the remaining pre-existing package-level items plus the Rust static-library `_abort` note; all package tests passed. The touched `commonfactorGWAS` and `userGWAS` Rd mismatch warnings are cleared.

## 2026-04-27 15:24 PDT - Prep-path parser and LDSC block benchmark

Change set:

- Added `benches/bench_prep_fast_paths.R` for `munge()`, `sumstats()`, `ldsc()` reader, and LDSC block-product timing.
- Added `options(GenomicSEM.fast_table_read=TRUE)` default path using `data.table::fread(check.names=TRUE)` with `read.table()` fallback.
- Added `options(GenomicSEM.fast_ldsc_read=TRUE)` default path using `data.table::fread()` for `ldsc()` chromosome, weight, M, and trait summary-stat ingestion.
- Reused the same reader helpers in `s_ldsc()` to replace `plyr::ldply()` during LD-score, M, and weight file ingestion.
- Serial `munge()` and `sumstats()` now stream trait files one at a time instead of preloading all files.
- Tried a custom SNP join helper behind `options(GenomicSEM.fast_snp_join=TRUE)`, but left it disabled by default because it was neutral/slightly slower locally.
- Centralized LDSC block-product code; the attempted cumulative-product variant was backed out after benchmarking because the original BLAS `crossprod()` loop is already faster and much more memory-stable for `s_ldsc()` with many annotations.

Validation:

```sh
R CMD INSTALL --install-tests .
Rscript tests/prep-fast-path.R
```

Local benchmark command:

```sh
Rscript benches/bench_prep_fast_paths.R 100000 2 1 200 1
```

Local result summary:

| workflow | fast_table_read | fast_snp_join | fast_ldsc_read | n_snp | n_traits | elapsed_sec | checksum |
|---|---|---|---|---:|---:|---:|---:|
| sumstats | FALSE | FALSE | NA | 100000 | 2 | 1.496 | 31028.76 |
| munge | FALSE | FALSE | NA | 100000 | 2 | 1.931 | 2000699000 |
| sumstats | TRUE | FALSE | NA | 100000 | 2 | 0.657 | 31028.76 |
| munge | TRUE | FALSE | NA | 100000 | 2 | 1.188 | 2000699000 |
| sumstats | FALSE | TRUE | NA | 100000 | 2 | 1.351 | 31028.76 |
| munge | FALSE | TRUE | NA | 100000 | 2 | 1.992 | 2000699000 |
| sumstats | TRUE | TRUE | NA | 100000 | 2 | 0.731 | 31028.76 |
| munge | TRUE | TRUE | NA | 100000 | 2 | 1.218 | 2000699000 |
| ldsc_block_old_loop | NA | NA | NA | 100000 | NA | 0.003 | 199734.4 |
| ldsc_block_fast | NA | NA | NA | 100000 | NA | 0.001 | 199734.4 |
| ldsc_read | NA | NA | FALSE | 100000 | 1 | 0.213 | 6001087000 |
| s_ldsc_read | NA | NA | FALSE | 100000 | NA | 0.102 | 5001479000 |
| ldsc_read | NA | NA | TRUE | 100000 | 1 | 0.113 | 6001087000 |
| s_ldsc_read | NA | NA | TRUE | 100000 | NA | 0.084 | 5001479000 |

Interpretation:

- The high-impact prep-path change is file ingestion and streaming: on this local 100k-SNP/2-trait synthetic run, `sumstats()` improves `1.496s -> 0.657s` and `munge()` improves `1.931s -> 1.188s`.
- The new `ldsc()` reader improves the synthetic chromosome/trait ingestion benchmark `0.213s -> 0.113s`; the analogous `s_ldsc()` file-list reader improves `0.102s -> 0.084s` on the same local synthetic files.
- The custom SNP join did not justify enabling; it remains opt-in for further experiments only.
- LDSC jackknife block crossproducts are not the current bottleneck for ordinary `ldsc()`; further LDSC work should focus on file ingestion/merging and stratified annotation data movement rather than replacing the existing `crossprod()` loop in R.

## 2026-04-27 15:39 PDT - Avoid full-length OR-detection ifelse

Change set:

- Replaced the `ifelse(rep(condition, nrow(file)), log(effect), effect)` pattern in `munge_main` and `sumstats_main` with a scalar branch.
- This avoids one full-length temporary vector for each trait file; it is mainly a memory/cleanup improvement rather than a clear wall-time shift.

Validation:

```sh
R CMD INSTALL --install-tests .
Rscript tests/prep-fast-path.R
Rscript benches/bench_prep_fast_paths.R 100000 2 1 200 1
```

Local result summary for the default fast reader path:

| workflow | elapsed_sec before-ish | elapsed_sec after | note |
|---|---:|---:|---|
| sumstats | 0.657 | 0.657 | unchanged within local noise |
| munge | 1.188 | 1.233 | unchanged/slightly noisy; gzip and table IO dominate |
| ldsc_read | 0.113 | 0.109 | unrelated reader path noise |
| s_ldsc_read | 0.084 | 0.076 | unrelated reader path noise |

## 2026-04-27 15:58 PDT - Rust row-QC kernels for munge and sumstats

Change set:

- Added native Rust row-QC kernels for supported `munge()` and `sumstats()` post-merge prep work.
- `munge()` now tries Rust for missing-value filtering, OR detection, allele flipping/matching, INFO/MAF filtering, and Z-score generation.
- `sumstats()` now tries Rust for missing-value filtering, MAF/varSNP handling, OR detection, allele flipping/matching, INFO filtering, Z-score generation, and the final beta/SE transforms for OLS, linear-probability, and logistic modes.
- Added `.Call` wrappers and options `GenomicSEM.fast_munge_qc` / `GenomicSEM.fast_sumstats_qc`, both enabled by default when the native library is available.
- Prep parity tests now compare legacy R QC against the Rust QC path. Z values are checked with tolerance because the Rust path uses an independent inverse-normal approximation for `sqrt(qchisq(P, 1, lower.tail = FALSE))`.

Validation:

```sh
cargo test --workspace
R CMD INSTALL --install-tests .
Rscript tests/prep-fast-path.R
Rscript benches/bench_prep_fast_paths.R 100000 2 1 200 1
```

Local benchmark summary:

| workflow | fast_table_read | fast_prep_qc | fast_snp_join | elapsed_sec | checksum |
|---|---|---|---|---:|---:|
| sumstats | FALSE | FALSE | FALSE | 1.376 | 31028.76 |
| munge | FALSE | FALSE | FALSE | 1.838 | 2000699000 |
| sumstats | TRUE | FALSE | FALSE | 0.591 | 31028.76 |
| munge | TRUE | FALSE | FALSE | 1.263 | 2000699000 |
| sumstats | FALSE | TRUE | FALSE | 0.923 | 31028.76 |
| munge | FALSE | TRUE | FALSE | 1.620 | 2000699000 |
| sumstats | TRUE | TRUE | FALSE | 0.185 | 31028.76 |
| munge | TRUE | TRUE | FALSE | 0.987 | 2000699000 |
| sumstats | TRUE | TRUE | TRUE | 0.240 | 31028.76 |
| munge | TRUE | TRUE | TRUE | 0.923 | 2000699000 |

Interpretation:

- The combined fread + Rust QC path improves the 100k-SNP/2-trait synthetic `sumstats()` run `1.376s -> 0.185s`, about `7.4x`.
- The same combined path improves `munge()` `1.838s -> 0.987s`, about `1.9x`; output gzip/write time is now a larger share of the remaining runtime.
- The experimental SNP join is still mixed: it helps `munge()` in this run but hurts `sumstats()`, so it remains disabled by default.
