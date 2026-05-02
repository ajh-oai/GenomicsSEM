# GenomicSEM Rust Fast-Path Update

This branch adds native Rust acceleration to the highest-cost GenomicSEM execution paths while
preserving the existing R-facing API and falling back to the original R/lavaan implementations for
unsupported inputs.

## What Changed

### Model fitting

- Added batched Rust SEM solvers for supported `commonfactorGWAS()` and `userGWAS()` DWLS models.
- Added a native `userGWAS()` fast-path compiler and initializer so supported Rust-backed models no
  longer require a lavaan setup pass before the batch fit.
- Added support for the practical lavaan syntax used by the public p-factor workflow, including
  inequality constraints such as `a > .001`.
- Kept the existing R result shape, `Q_SNP` behavior, `sub=` extraction, and lavaan fallback path.
- Replaced expensive per-SNP R output assembly in the Rust-backed `userGWAS(sub=...)` path with a
  vectorized assembler.

### Summary-stat preparation and LDSC setup

- Added native Rust prep engines for supported `munge()` and `sumstats()` inputs.
- Added batched native `munge()` orchestration for multiple files, sharing one reference index and
  parallelizing compressed output when multiple Rust worker threads are requested.
- Added batched native `sumstats()` orchestration for multiple files, sharing one reference index
  and performing the final cross-trait listwise SNP merge before returning to R.
- Added faster table ingestion and native QC/transformation kernels for supported prep paths.
- Added faster `ldsc()` / `s_ldsc()` readers plus a Rust LDSC block-crossproduct kernel.

### Packaging and validation surface

- Added Cargo workspaces, native R bindings, registration code, package docs, parity tests, release
  tests, benchmark harnesses, and reproducible workflow harnesses under `repro/`.
- Added diagnostics and options so callers can see whether a run used the Rust path or the legacy
  fallback.

## Why It Matters

The largest gain comes from removing the old "fit one lavaan model per SNP" loop from supported
`userGWAS()` models. That loop was spending most of its time in repeated R/lavaan orchestration, not
in irreducible SEM mathematics.

On a remote 16-CPU panda pod, the public p-factor practical model scaled as follows:

| Path | Workload | Runtime |
|---|---:|---:|
| old R/lavaan, fitted slope from 1k-50k SNPs | projected 1M SNPs | 3.363 h |
| Rust-backed `userGWAS()` | measured 1M SNPs | 12.380 s |

That is a projected `978x` model-fitting speedup on the same class of machine. The old-path fit is
based on a linear regression with `R^2 = 0.999996`; the full old 1M run was not executed because the
bounded scaling run already made the cost clear.

Other measured gains:

| Workflow | Benchmark | Old | New | Speedup |
|---|---|---:|---:|---:|
| `sumstats()` | public p-factor data, 1,291,369 aligned SNPs, remote pod | 40.821 s | 9.816 s | 4.2x |
| `sumstats()` | synthetic 100k SNP / 2-trait prep benchmark | 1.322 s | 0.139 s | 9.5x |
| `sumstats()` | synthetic 100k SNP / 4-trait batched prep benchmark, clean-clone median | 2.546 s | 0.079 s | 32.2x |
| `munge()` | synthetic 100k SNP / 2-trait prep benchmark | 1.901 s | 0.433 s | 4.4x |
| `munge()` | synthetic 100k SNP / 4-trait batched prep benchmark, clean-clone median | 3.834 s | 0.376 s | 10.2x |

## What Was Reproduced

The branch includes three workflow harnesses:

1. `repro/grotzinger_2019_nhb.R`
   - Reproduces the published four-trait no-SNP common-factor example exactly from the public LDSC
     matrices.
   - Benchmarks the corresponding SNP-effect workflow shape with old R/lavaan versus Rust-backed R.

2. `repro/grotzinger_2025_nature.R`
   - Recreates the public 14-disorder workflow shape from the reported summary information.
   - This is a same-input performance/equivalence harness, not a full raw-input paper rerun, because
     the complete univariate inputs and full sampling covariance matrix are not public.

3. `repro/pfactor_practical_1m.R`
   - Runs the public GenomicSEM p-factor practical shape using public SCZ, BIP, and MDD summary
     statistics plus `GenomicSEMPractical.RData`.
   - Builds an aligned real-data 1M-SNP table, checks bounded old-vs-new equivalence, and runs the
     full Rust-backed 1M scan.

For the public p-factor benchmark, bounded old-vs-new checks matched to about `1e-6` for
`userGWAS()` and near machine precision for `sumstats()`.

## Current Scope

The Rust acceleration is intentionally selective.

- Accelerated: supported `munge()`, `sumstats()`, LDSC setup kernels, `commonfactorGWAS()`, and
  `userGWAS()` fast-path model classes.
- Preserved: existing R APIs, old implementations, and fallback semantics for unsupported syntax.
- Not claimed: a complete replacement of every lavaan model or every GenomicSEM function.

Use `options(GenomicSEM.fast_diagnostics=TRUE)` or inspect the `GenomicSEM.fast_path` result
attribute to confirm which path a run used.

## Evidence and Details

- Full benchmark history: [`benches/BENCHMARK_LOG.md`](benches/BENCHMARK_LOG.md)
- Reproduction harnesses and caveats: [`repro/README.md`](repro/README.md)
- User-facing patch notes: [`PATCHNOTES.md`](PATCHNOTES.md)

## Release Work Still Worth Doing

- Add clean-checkout CI for Linux, macOS, Windows, `R CMD check`, and `cargo test`.
- Expand golden compatibility fixtures across supported fast-path syntax and fallback cases.
- Decide which fast paths should be default-on versus explicitly experimental in the first public
  release.
- Re-run the public benchmark suite from a clean clone and publish one concise release benchmark
  report.
