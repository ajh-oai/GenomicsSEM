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
