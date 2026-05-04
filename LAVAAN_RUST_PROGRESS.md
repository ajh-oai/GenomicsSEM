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
