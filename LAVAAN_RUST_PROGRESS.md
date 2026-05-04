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
