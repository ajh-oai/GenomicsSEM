# lavaan_rust progress log

## 2026-05-04 12:34 PDT

### Implemented

- Added the nested experimental R package `lavaanrust/`.
- Added an `extendr` Rust backend with `fit_one_factor_dwls()`.
- Added a minimal S4 compatibility object, `lavaan_rust_fit`, and the first
  suffixed compatibility surface:
  - `sem_rust()`
  - `lavaan_rust()`
  - `lavInspect_rust()`
  - `inspect_rust()`
  - `parTable_rust()`
  - `fitted_rust()`
  - `resid_rust()`
  - defined-parameter helper stubs that currently strict-error
- Added `commonfactor_rust()` without copying the original function body:
  it clones the existing `commonfactor()` closure and rebinds only the lavaan
  symbols in a child environment.
- Added a second native closed-form DWLS slice for observed-covariance models,
  covering the common-factor CFI/null model and its parameter-table refit.
- Switched `_rust()` policy to strict mode: unsupported paths now error instead
  of silently delegating to upstream lavaan.
- Added frozen one-factor fixtures from lavaan 0.6-21 and test coverage for:
  - the rust-backed one-factor DWLS slice
  - strict errors for unsupported syntax

### Current support boundary

`sem_rust()` currently owns:

- one-factor covariance-only DWLS models used by the primary and standardized
  `commonfactor()` fits
- the observed-covariance DWLS family used by the common-factor CFI/null model
  and its parameter-table refit

Unsupported syntax now errors deliberately; `_rust()` wrappers are meant to
measure the Rust surface that actually exists, not hide missing coverage by
delegating back to lavaan.

### Validation

- `lavaanrust/tests/testthat`: 12 passing tests.
- `R CMD build lavaanrust` followed by
  `R CMD check lavaanrust_0.0.0.9000.tar.gz --no-manual`: `Status: OK`.
- Synthetic `commonfactor()` smoke comparison:
  - max unstandardized estimate difference: `3.74e-09`
  - max sandwich SE difference: `3.72e-09`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct one-factor `sem()` fit | `0.031630 s` | `0.000828 s` | `38.20x` |
| end-to-end `commonfactor()` smoke, first packet | `0.164460 s` | `0.064800 s` | `2.54x` |
| direct one-factor `sem()` fit after null-model slice | `0.017726 s` | `0.000522 s` | `33.96x` |
| end-to-end `commonfactor()` smoke after null-model slice | `0.092230 s` | `0.006320 s` | `14.59x` |

The direct-fit gain is the clean measurement of the Rust one-factor backend.
The second end-to-end benchmark shows the impact of removing the remaining
lavaan work from `commonfactor_rust()` itself.

### Next packet

1. Add strict base-model reuse support through `lavaan_rust()` for GWAS loops.
2. Broaden `sem_rust()` from current covariance-only families toward the
   `commonfactorGWAS()` model slice.
3. Add `commonfactorGWAS_rust()` once the model-reuse path is native.

## 2026-05-04 12:59 PDT

### Implemented

- Added a strict native `commonfactorGWAS()` backend family:
  - marker-scaled factor/SNP DWLS first-stage fit
  - slot-style base-model reuse through `lavaan_rust()`
  - Q-model refit with direct SNP effects and trait residual variances free
- Extended `lavaan_rust_fit` with the compatibility slots needed by the current
  GenomicSEM reuse path: `Options`, `Data`, and `Model`.
- Added `commonfactorGWAS_rust()` using the original GenomicSEM implementation
  with rust-bound helper closures for `.commonfactorGWAS_main()` and
  `.rearrange()`.
- Kept the wrapper boundary strict: `commonfactorGWAS_rust()` currently requires
  `parallel = FALSE`; unsupported parallel execution errors instead of drifting
  into the original worker path.
- Added a narrow `class()` compatibility shim inside the rust-bound wrapper
  environment so the unchanged upstream `class(fit)[1] == "lavaan"` checks see
  rust fit objects as successful lavaan-compatible fits.

### Current support boundary

`sem_rust()` now owns:

- one-factor covariance-only DWLS models used by `commonfactor()`
- the observed-covariance DWLS family used by the common-factor CFI/null model
  and its parameter-table refit
- the marker-scaled one-factor SNP model used by `commonfactorGWAS()`
- the `commonfactorGWAS()` Q-model parameter-table refit

`lavaan_rust()` now owns the strict slot-reuse path for the supported
`commonfactorGWAS()` DWLS base model. Unsupported reuse patterns still error and
do not delegate to lavaan.

`commonfactorGWAS_rust()` currently supports the sequential path only. The
existing `foreach` worker setup needs an explicit rust-bound export strategy
before the parallel path can be considered supported.

### Validation

- `lavaanrust/tests/testthat`: 20 passing tests.
- Frozen lavaan 0.6-21 `commonfactorGWAS()` fixture:
  - first-stage parameter estimates match within `2e-06`
  - Q-model refit parameter estimates match within `2e-06`
- Synthetic sequential GenomicSEM smoke:
  - `commonfactorGWAS()` vs `commonfactorGWAS_rust()` factor-SNP estimate max
    absolute difference: `1.57e-08`
  - Q-statistic absolute difference: `4.82e-08`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct `commonfactorGWAS()` first-stage `sem()` fit | `0.018000 s` | `0.001000 s` | `18.00x` |
| end-to-end sequential `commonfactorGWAS()` smoke, 1 SNP | `0.067000 s` | `0.005000 s` | `13.40x` |

### Next packet

1. Decide whether to push the experiment deeper into `commonfactorGWAS()`
   parallel worker support or move next to `usermodel_rust()`.
2. Broaden parser and constraint support only where another GenomicSEM wrapper
   actually needs it.

## 2026-05-04 13:20 PDT

### Implemented

- Added the first strict `usermodel_rust()` slice while keeping the original
  GenomicSEM `usermodel()` body intact.
- Added support for simple one-factor user models such as
  `F1 =~ A + B + C`:
  - marker scaling when `std.lv = FALSE`
  - latent-variance scaling when `std.lv = TRUE`
- Added the minimal `standardizedSolution_rust()` compatibility surface needed
  by the unchanged `usermodel()` standardized-output path.
- Split the implicit `std.lv = TRUE` constructor from the explicit
  `F1 ~~ 1*F1` constructor so parameter-table row order matches lavaan for both
  syntactic forms.
- Tightened the class shim used by rust-bound wrappers so exact upstream checks
  such as `class(x) != "lavaan"` continue to behave like the original code.

### Current support boundary

`usermodel_rust()` currently supports simple one-factor DWLS models with no:

- labels or equality constraints
- inequality constraints
- user-defined parameters via `:=`
- regressions
- multiple latent factors

Unsupported user-model syntax still errors through `sem_rust()` and does not
fall back to lavaan.

### Validation

- `lavaanrust/tests/testthat`: 28 passing tests.
- Frozen lavaan 0.6-21 fixtures for:
  - marker-scaled one-factor user models
  - `std.lv = TRUE` one-factor user models
  - standardized-solution rows consumed by `usermodel()`
- Synthetic `usermodel()` smoke, simple marker-scaled one-factor model:
  - max unstandardized estimate difference: `7.17e-09`
  - max sandwich SE difference: `1.60e-09`
  - max `STD_Genotype` difference: `1.11e-08`
  - max `STD_All` difference: `5.97e-09`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct user-model `sem()` fit | `0.017000 s` | `0.001000 s` | `17.00x` |
| end-to-end simple `usermodel()` smoke | `0.067000 s` | `0.005000 s` | `13.40x` |

### Next packet

1. Expand the user-model parser one feature at a time, starting with the
   smallest syntax family that unlocks a real GenomicSEM workflow beyond the
   current one-factor case.
2. Likely next targets are either constrained one-factor models or the
   user-defined-parameter machinery needed for `:=`.

## 2026-05-04 13:30 PDT

### Implemented

- Added the first strict `userGWAS_rust()` slice while preserving the original
  GenomicSEM `userGWAS()` body.
- Reused the existing marker-scaled one-factor SNP optimizer for the
  unrestricted first-stage `F1 =~ ...` plus `F1 ~ SNP` fit.
- Added a new native fixed-measurement DWLS solver for the default
  `fix_measurement = TRUE` path:
  - fixed loadings from the no-SNP measurement model
  - free trait residual variances
  - free latent residual variance
  - free factor-SNP regression
  - free SNP variance
- Added strict native slot reuse for the fixed-measurement base model through
  `lavaan_rust()`, which is the path exercised inside the per-SNP loop.
- Kept the wrapper boundary explicit. `userGWAS_rust()` currently requires:
  - `parallel = FALSE`
  - `fix_measurement = TRUE`
  - `Q_SNP = FALSE`
  - `estimation = "DWLS"`
  - `TWAS = FALSE`

### Current support boundary

The new `userGWAS_rust()` slice supports simple one-factor models such as:

```r
F1 =~ A + B + C
F1 ~ SNP
```

Unsupported user-GWAS paths still error and do not fall back to lavaan. The
current packet intentionally does not yet cover unconstrained measurement
refits, `Q_SNP`, TWAS, parallel workers, multiple latent factors, or more
general user-model syntax.

### Validation

- `lavaanrust/tests/testthat`: 37 passing tests.
- Frozen lavaan 0.6-21 fixtures for:
  - the unrestricted one-factor user-GWAS first stage
  - the fixed-measurement parameter-table refit
  - fixed-measurement slot reuse
- Synthetic sequential GenomicSEM smoke, 1 SNP:
  - max estimate difference: `5.20e-09`
  - max sandwich SE difference: `3.34e-10`
  - max chi-square difference: `2.06e-09`
- `R CMD build lavaanrust` followed by
  `R CMD check lavaanrust_0.0.0.9000.tar.gz --no-manual`: clean after adding
  the generated Rd page for the new exported solver.

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct fixed-measurement `sem()` refit | `0.015040 s` | `0.000244 s` | `61.64x` |
| end-to-end sequential `userGWAS()` smoke, 1 SNP | `0.105400 s` | `0.013300 s` | `7.92x` |

### Next packet

1. Decide whether to broaden `userGWAS_rust()` next into `Q_SNP = TRUE`, the
   unconstrained `fix_measurement = FALSE` path, or parallel worker support.
2. If the goal is still the pure isolated-backend experiment, the most
   informative next step is probably a second real user-model family rather
   than changing outer GenomicSEM orchestration.

## 2026-05-04 13:41 PDT

### Implemented

- Expanded `userGWAS_rust()` from the first default packet to the full current
  one-factor sequential DWLS slice:
  - both `fix_measurement = TRUE/FALSE`
  - both `Q_SNP = TRUE/FALSE`
- Added native `lavaan_rust()` slot reuse for the unrestricted
  `user_gwas_dwls` model family used when `fix_measurement = FALSE`.
- Kept the supported matrix explicit in the `userGWAS_rust()` source comments
  so the strict wrapper boundary is visible where the behavior is enforced.

### Current support boundary

`userGWAS_rust()` now supports simple one-factor SNP models such as:

```r
F1 =~ A + B + C
F1 ~ SNP
```

with:

- `parallel = FALSE`
- `estimation = "DWLS"`
- `TWAS = FALSE`
- either `fix_measurement = TRUE/FALSE`
- either `Q_SNP = TRUE/FALSE`

Still unsupported:

- `parallel = TRUE`
- `TWAS = TRUE`
- ML estimation
- multi-factor or more general user-model syntax

### Validation

- `lavaanrust/tests/testthat`: 38 passing tests.
- Added native reuse coverage for the unrestricted `user_gwas_dwls` base model.
- Synthetic sequential GenomicSEM smokes, 1 SNP:
  - `fix_measurement = FALSE`, `Q_SNP = FALSE`:
    - max estimate difference: `1.04e-07`
    - max sandwich SE difference: `4.89e-08`
    - max chi-square difference: `1.16e-13`
  - `fix_measurement = FALSE`, `Q_SNP = TRUE`:
    - max `Q_SNP` difference: `3.44e-08`
    - max `Q_SNP` p-value difference: `1.65e-08`
  - `fix_measurement = TRUE`, `Q_SNP = TRUE`:
    - max `Q_SNP` difference: `2.06e-09`
    - max `Q_SNP` p-value difference: `9.86e-10`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct unrestricted `userGWAS()` `sem()` fit | `0.021660 s` | `0.000732 s` | `29.59x` |
| end-to-end sequential, `fix_measurement = FALSE`, `Q_SNP = FALSE` | `0.062500 s` | `0.011800 s` | `5.30x` |
| end-to-end sequential, `fix_measurement = FALSE`, `Q_SNP = TRUE` | `0.057200 s` | `0.005300 s` | `10.79x` |
| end-to-end sequential, `fix_measurement = TRUE`, `Q_SNP = FALSE` | `0.100500 s` | `0.006600 s` | `15.23x` |
| end-to-end sequential, `fix_measurement = TRUE`, `Q_SNP = TRUE` | `0.103800 s` | `0.007100 s` | `14.62x` |

### Next packet

1. Add parallel worker support for the now-complete current one-factor
   `userGWAS_rust()` slice.
2. Then decide whether to broaden `usermodel_rust()` syntax or move to a
   second real user-GWAS model family.

## 2026-05-04 15:18 PDT

### Implemented

- Added explicit parallel worker plumbing for backend-swapped GWAS wrappers:
  - `commonfactorGWAS()` and `userGWAS()` now bind the worker helper before
    entering `foreach`
  - worker package exports are selected from wrapper-local settings, so
    `_rust()` variants load `lavaanrust` workers without changing the original
    serial logic
- Enabled strict native parallel execution for:
  - `commonfactorGWAS_rust()`
  - `userGWAS_rust()`
- Added `tools/bench_lavaan_rust_parallel.R`, a reproducible benchmark harness
  that:
  - generates deterministic synthetic LDSC/GWAS inputs
  - checks sequential-vs-parallel equality for the rust-backed wrappers
  - benchmarks lavaan and rust-backed wrappers over multiple core counts

### Current support boundary

`commonfactorGWAS_rust()` now supports the current one-factor DWLS model family
with either `parallel = TRUE/FALSE`.

`userGWAS_rust()` now supports simple one-factor SNP models such as:

```r
F1 =~ A + B + C
F1 ~ SNP
```

with:

- `parallel = TRUE/FALSE`
- `estimation = "DWLS"`
- `TWAS = FALSE`
- `std.lv = FALSE`
- either `fix_measurement = TRUE/FALSE`
- either `Q_SNP = TRUE/FALSE`

Still unsupported:

- `TWAS = TRUE`
- `std.lv = TRUE`
- ML estimation
- multi-factor or more general user-model syntax

### Validation

- Local 4-SNP and benchmark-harness smoke checks matched exactly between
  rust-backed sequential and parallel output on the checked numeric columns.
- Remote panda/flex 16-CPU benchmark pod, synthetic `n_snp = 1000`, `k = 5`
  fixture:
  - `commonfactorGWAS_rust()` sequential vs 2-worker parallel max absolute
    difference: `0`
  - `userGWAS_rust()` sequential vs 2-worker parallel max absolute difference:
    `0`

### Remote benchmarks

Synthetic 1,000-SNP, 5-trait fixture on a panda/flex 16-CPU pod,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
2026-05-04:

| Workflow | Backend | sequential | 2 cores | 4 cores | 8 cores | 16 cores |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `commonfactorGWAS()` | lavaan | `167.523 s` | `85.818 s` | `45.065 s` | `24.809 s` | `18.668 s` |
| `commonfactorGWAS_rust()` | rust-backed | `6.227 s` | `4.432 s` | `2.091 s` | `1.509 s` | `1.365 s` |
| `userGWAS()` | lavaan | `62.443 s` | `62.423 s` | `28.090 s` | `17.820 s` | `14.530 s` |
| `userGWAS_rust()` | rust-backed | `5.439 s` | `4.581 s` | `2.704 s` | `2.138 s` | `2.361 s` |

Selected comparisons:

- `commonfactorGWAS_rust()` is `26.90x` faster than lavaan sequentially and
  `13.68x` faster at 16 workers.
- `userGWAS_rust()` is `11.48x` faster than lavaan sequentially and `6.15x`
  faster at 16 workers.
- The common-factor rust path continues to improve through 16 workers on this
  fixture; the user-GWAS rust path is fastest at 8 workers here, so worker
  overhead is already visible once the per-SNP fit cost is small enough.

### Next packet

1. Broaden `usermodel_rust()` syntax or move to a second real user-GWAS model
   family.
2. If parallel performance remains a priority later, revisit task granularity
   and batching for the very fast rust-backed loops rather than only adding
   more workers.

## 2026-05-04 16:05 PDT

### Implemented

- Expanded `tools/bench_lavaan_rust_parallel.R` from the default
  fixed-measurement / `Q_SNP = TRUE` benchmark into the full current supported
  `userGWAS_rust()` matrix:
  - `fix_measurement = TRUE/FALSE`
  - `Q_SNP = TRUE/FALSE`
- Added full-matrix equivalence checks before timing:
  - lavaan vs rust-backed sequential output for every supported user-GWAS cell
  - rust-backed sequential vs 2-worker parallel output for every supported
    user-GWAS cell

### Validation

Remote panda/flex 16-CPU pod, synthetic `n_snp = 1000`, `k = 5` fixture:

| Workflow | `fix_measurement` | `Q_SNP` | Comparison | Max absolute difference |
| --- | --- | --- | --- | ---: |
| `commonfactorGWAS()` | `NA` | `NA` | rust sequential vs parallel | `0` |
| `userGWAS()` | `TRUE` | `TRUE` | lavaan vs rust sequential | `1.18e-07` |
| `userGWAS()` | `TRUE` | `TRUE` | rust sequential vs parallel | `0` |
| `userGWAS()` | `FALSE` | `TRUE` | lavaan vs rust sequential | `1.38e-07` |
| `userGWAS()` | `FALSE` | `TRUE` | rust sequential vs parallel | `0` |
| `userGWAS()` | `TRUE` | `FALSE` | lavaan vs rust sequential | `1.18e-07` |
| `userGWAS()` | `TRUE` | `FALSE` | rust sequential vs parallel | `0` |
| `userGWAS()` | `FALSE` | `FALSE` | lavaan vs rust sequential | `1.38e-07` |
| `userGWAS()` | `FALSE` | `FALSE` | rust sequential vs parallel | `0` |

### Remote benchmarks

Synthetic 1,000-SNP, 5-trait fixture on the same panda/flex 16-CPU pod,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
2026-05-04:

| Workflow | Backend | `fix_measurement` | `Q_SNP` | sequential | 2 cores | 4 cores | 8 cores | 16 cores |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `commonfactorGWAS()` | lavaan | `NA` | `NA` | `167.093 s` | `85.649 s` | `44.774 s` | `24.852 s` | `18.423 s` |
| `commonfactorGWAS_rust()` | rust-backed | `NA` | `NA` | `6.250 s` | `4.941 s` | `2.227 s` | `1.504 s` | `1.463 s` |
| `userGWAS()` | lavaan | `TRUE` | `TRUE` | `62.121 s` | `62.418 s` | `28.771 s` | `17.526 s` | `14.421 s` |
| `userGWAS_rust()` | rust-backed | `TRUE` | `TRUE` | `5.411 s` | `5.484 s` | `2.819 s` | `2.176 s` | `2.283 s` |
| `userGWAS()` | lavaan | `FALSE` | `TRUE` | `87.534 s` | `77.286 s` | `35.251 s` | `20.281 s` | `14.901 s` |
| `userGWAS_rust()` | rust-backed | `FALSE` | `TRUE` | `6.589 s` | `6.389 s` | `2.928 s` | `2.438 s` | `2.283 s` |
| `userGWAS()` | lavaan | `TRUE` | `FALSE` | `61.323 s` | `58.999 s` | `28.431 s` | `18.290 s` | `13.912 s` |
| `userGWAS_rust()` | rust-backed | `TRUE` | `FALSE` | `4.472 s` | `4.008 s` | `2.455 s` | `2.122 s` | `1.937 s` |
| `userGWAS()` | lavaan | `FALSE` | `FALSE` | `87.180 s` | `81.775 s` | `35.065 s` | `20.161 s` | `16.304 s` |
| `userGWAS_rust()` | rust-backed | `FALSE` | `FALSE` | `5.341 s` | `4.045 s` | `2.651 s` | `2.269 s` | `2.036 s` |

Selected comparisons:

- Sequential speedups span `11.48x` to `16.32x` across the full supported
  `userGWAS_rust()` matrix.
- At 16 workers, speedups span `6.32x` to `8.01x` across that matrix.
- `Q_SNP = FALSE` lowers rust-backed runtime materially in both measurement
  modes, especially for `fix_measurement = TRUE`.
- The best rust-backed worker count is not constant:
  - `fix_measurement = TRUE`, `Q_SNP = TRUE`: best at 8 workers
  - `fix_measurement = FALSE`, `Q_SNP = TRUE`: best at 16 workers
  - `fix_measurement = TRUE`, `Q_SNP = FALSE`: best at 16 workers
  - `fix_measurement = FALSE`, `Q_SNP = FALSE`: best at 16 workers

### Next packet

1. Broaden `usermodel_rust()` syntax or move to a second real user-GWAS model
   family.
2. If more parallel work is desired later, target batching/task granularity for
   the already-fast rust loops rather than expecting uniform gains from simply
   raising `cores`.

## 2026-05-04 16:13 PDT

### Tightening

- Added an explicit `std.lv = TRUE` rejection in `userGWAS_rust()`.
- This was already outside the intended native slice, but it was not previously
  guarded:
  - `fix_measurement = TRUE, std.lv = TRUE` failed mid-run
  - `fix_measurement = FALSE, std.lv = TRUE` completed but was not equivalent to
    lavaan
- The wrapper now keeps the experiment's strict contract: unsupported paths
  fail immediately instead of producing mixed or incorrect results.

## 2026-05-04 17:41 PDT

### Generic compiler scaffold

- Added `LAVAAN_FAST_PLAN.md` to make the next phase explicit: `lavaan_fast`
  should begin as a generic compiler layer rather than as another family of
  hand-written fast paths.
- Added a parameter-table compiler in `lavaanrust/R/compiler.R` that lowers the
  current supported subset (`=~`, `~`, `~~`) into a RAM representation:
  - directed matrix `A`
  - residual covariance matrix `S`
  - observed-variable selector `F`
- Added generic reconstruction of:
  - implied covariance `Sigma = F (I - A)^-1 S (I - A)^-T F^T`
  - analytic `d vech(Sigma) / d theta` Jacobians
- Kept unsupported operators strict-errors at the compiler boundary; `:=` still
  does not silently pass through.

### Validation

- Added compiler tests for both current user-GWAS shapes:
  - unrestricted one-factor model
  - fixed-measurement one-factor model
- The compiler now exactly reproduces the specialized rust-backed fit objects
  for both:
  - implied covariance matrices
  - Jacobian matrices exposed through `lavInspect_rust(..., "delta")`
- Local package validation after the compiler packet:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `46` passing tests

### Next packet

1. Route one existing family through the generic compiler without changing its
   public wrapper contract.
2. Then widen the syntax/compiler surface toward the real blocker set for
   broader `userGWAS_rust()` coverage:
   - labels and equality reuse
   - direct effects
   - residual covariances
   - eventually defined parameters and nonlinear constraints

## 2026-05-04 18:50 PDT

### Native compiler evaluator

- Added `evaluate_ram_surfaces()` in Rust and the R bridge
  `.lavaan_fast_implied_surfaces_rust()`.
- The native evaluator consumes the compiler's row-wise RAM encoding and
  returns both:
  - implied observed covariance
  - analytic Jacobian over `vech(Sigma)`
- Added fixture checks proving the native evaluator matches the specialized
  rust-backed `userGWAS` surfaces for both current model shapes.

### Validation and profiling

- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `50` passing tests
- Direct evaluator timing on the unrestricted user-GWAS fixture:
  - R reference evaluator, 10,000 calls: `2.274 s`
  - Rust evaluator, 10,000 calls: `1.004 s`
- I also tried routing the unrestricted `userGWAS` constructor through the
  R-side compiler surfaces directly. That was behaviorally exact but slower:
  - baseline specialized constructor path, 1,000 tiny fits: `0.842 s`
  - R compiler-backed constructor path, 1,000 tiny fits: `1.134 s`
- I did not keep that hot-path change. The next production-quality step is a
  native generic optimizer that reuses compiler-backed surfaces during fitting,
  instead of paying an extra post-fit reconstruction cost after a specialized
  optimizer has already computed them.

### Next packet

1. Add a generic native DWLS optimizer over the RAM IR.
2. Route one existing supported family through that optimizer only if the
   benchmark remains competitive with the specialized kernel.

## 2026-05-04 18:53 PDT

### Generic RAM optimizer

- Added `fit_ram_dwls()` in Rust and the R bridge
  `.lavaan_fast_fit_dwls_rust()`.
- The generic optimizer now:
  - consumes the same compiled RAM rows as the native evaluator
  - reuses generic implied covariance and Jacobian surfaces during fitting
  - clamps free diagonal covariance parameters positive
  - returns estimates, implied covariance, Jacobian, naive SEs, fit objective,
    SRMR, convergence, and iteration count
- Added a fixed-measurement `userGWAS` fixture test showing the generic optimizer
  reproduces the existing specialized rust-backed fit.

### Validation and benchmark

- Local package validation:
  - `R CMD INSTALL lavaanrust`
  - `Rscript -e 'testthat::test_dir("lavaanrust/tests/testthat")'`
  - result: `53` passing tests
- Direct low-level benchmark on the fixed-measurement user-GWAS fixture:
  - specialized kernel, 1,000 fits: `0.014 s`
  - generic RAM optimizer, 1,000 fits: `0.109 s`
- Conclusion:
  - the generic optimizer is now a viable compiler-backed coverage path
  - it is not yet competitive with the tiny analytic kernels, so those should
    remain the production hot paths for already-supported models

### Next packet

1. Broaden the compiler syntax surface enough to unlock a real currently
   unsupported model shape.
2. Use the generic optimizer there first, where the comparison is against
   unsupported behavior rather than against an already-optimal specialized
   kernel.
