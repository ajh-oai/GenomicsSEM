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
  - defined-parameter helper stubs that currently delegate to lavaan
- Added `commonfactor_rust()` without copying the original function body:
  it clones the existing `commonfactor()` closure and rebinds only the lavaan
  symbols in a child environment.
- Added a second native closed-form DWLS slice for observed-covariance models,
  covering the common-factor CFI/null model and its parameter-table refit.
- Switched `_rust()` policy to strict mode: unsupported paths now error instead
  of silently delegating to upstream lavaan.
- Added frozen one-factor fixtures from lavaan 0.6-21 and test coverage for:
  - the rust-backed one-factor DWLS slice
  - fallback to upstream lavaan for unsupported syntax

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
