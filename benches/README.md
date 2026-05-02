# Benchmarks

This directory has two benchmark layers.

`BENCHMARK_LOG.md` records timestamped benchmark runs, code changes, commands, and interpretation over time.

`compare_backends.R` isolates the rewritten `.get_V_SNP` kernel. It compares:

- `old_r_loop`: original R helper in a SNP loop.
- `r_binding_loop`: R loop calling the Rust-backed `.Call` helper once per SNP.
- `r_binding_batch`: one R call into the Rayon batch kernel.
- `rust_loop` and `rust_batch`: the pure Rust CLI benchmark.

This benchmark is intentionally narrow. It is useful for validating kernel throughput, R-to-native call overhead, and Rayon scaling, but it does not represent a full GenomicSEM workflow.

`bench_usergwas_synthetic.R` exercises the package workflows that upstream performance notes emphasize: `userGWAS()` and `commonfactorGWAS()`. The synthetic fixture produced by `synthetic_inputs.R` matches the documented call shape: an LDSC-like `covstruc` list, a `sumstats`-like SNP table, and a lavaan model with SNP effects. It runs the same workflow twice by toggling `options(GenomicSEM.use_rust)`, so the model fitting, lavaan setup, output assembly, and existing `foreach` parallelism remain in the measurement.

`profile_workflows.R` runs an Rprof profile for a small end-to-end workflow. `profile_commonfactor_phases.R` is a lower-overhead manual phase split for `commonfactorGWAS()`; use it when Rprof is too invasive or unstable at larger trait counts.

`bench_munge_batch.R` and `bench_sumstats_batch.R` compare the old serial/PSOCK prep orchestration
against the batched native prep engines on synthetic multi-file workloads. The `sumstats()` harness
also introduces per-trait missing P-values so the native cross-trait listwise merge is exercised,
not just the per-file transforms.

The upstream repo does not include a runnable benchmark suite or bundled benchmark data. Its existing performance references are in `README.md` and `PATCHNOTES.md`, including 100K-SNP and 1.8M-SNP `userGWAS` timings. This synthetic workflow benchmark is meant to be a reproducible stand-in for those workloads, not a biological validation dataset.

## Current Results

The current clean-clone benchmark medians show the largest prep-engine gains once work is batched
across traits and moved out of the R orchestration layer:

| Workflow | Old path | New path | Speedup |
|---|---:|---:|---:|
| `sumstats()`, 4 files x 100k SNPs | `2.449s` legacy serial | `0.076s` native batch, 4 threads | `32.2x` |
| `munge()`, 4 files x 100k SNPs | `3.628s` legacy serial | `0.375s` native batch, 4 threads | `9.7x` |

The synthetic SEM workflow harnesses also show that the Rust-backed fitting paths remove most of
the old lavaan-loop cost while preserving the same workflow shape:

| Workflow | Old path | New path | Speedup |
|---|---:|---:|---:|
| `userGWAS(Q_SNP=TRUE)`, 1 core/thread | `36.321s` | `0.669s` | `54.3x` |
| `userGWAS(Q_SNP=TRUE)`, 4 cores/threads | `10.415s` | `0.522s` | `20.0x` |
| `commonfactorGWAS()`, 1 core/thread | `86.470s` | `0.228s` | `379x` |
| `commonfactorGWAS()`, 4 cores/threads | `24.045s` | `0.059s` | `408x` |

These are synthetic, reproducible workflow benchmarks. The public-data and paper-shaped workflow
comparisons are summarized separately in [`repro/README.md`](../repro/README.md), and the full
timestamped history lives in [`BENCHMARK_LOG.md`](BENCHMARK_LOG.md).

Example commands:

```sh
Rscript benches/compare_backends.R 200000 8 standard 1,2,4,8,16
Rscript benches/bench_usergwas_synthetic.R 1000 6 1,2,4,8 userGWAS 1
Rscript benches/bench_usergwas_synthetic.R 1000 6 1,2,4,8 commonfactorGWAS 1
Rscript benches/profile_commonfactor_phases.R 50 12 1
Rscript benches/profile_workflows.R 100 3 userGWAS 1
Rscript benches/bench_munge_batch.R 100000 4 4 13
Rscript benches/bench_sumstats_batch.R 100000 4 4 13
```

`bench_usergwas_synthetic.R` accepts optional optimization toggles after the seed:

```sh
Rscript benches/bench_usergwas_synthetic.R <n_snp> <k> <cores> <workflow> <seed> <fast_diag> <fast_wls> <printwarn> <q_snp> <fast_commonfactor_fit> <fast_usergwas_fit>
```

Each toggle can be a comma-separated list such as `FALSE,TRUE`. `printwarn`, `q_snp`, and `fast_usergwas_fit` only apply to `userGWAS()`. `fast_commonfactor_fit` only applies to `commonfactorGWAS()`. The fast common-factor fit is an experimental batched Rust Gauss-Newton replacement for the main lavaan fit and the direct-effect Q model. The fast userGWAS fit is an experimental batched Rust RAM-matrix solver for supported DWLS parameter-table models.
