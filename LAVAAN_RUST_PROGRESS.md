# lavaan_rust progress log

## 2026-05-04

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
- Added frozen one-factor fixtures from lavaan 0.6-21 and test coverage for:
  - the rust-backed one-factor DWLS slice
  - fallback to upstream lavaan for unsupported syntax

### Current support boundary

`sem_rust()` currently owns only the covariance-only one-factor DWLS model used
by the primary and standardized `commonfactor()` fits. Unsupported syntax falls
back to upstream lavaan intentionally, so `commonfactor_rust()` can already run
end-to-end while the parser and compatibility surface are broadened.

The CFI/null-model fits inside `commonfactor_rust()` still use lavaan today.

### Validation

- `lavaanrust/tests/testthat`: 7 passing tests.
- Synthetic `commonfactor()` smoke comparison:
  - max unstandardized estimate difference: `3.74e-09`
  - max sandwich SE difference: `3.72e-09`

### Local benchmarks

Synthetic 3-trait fixture, macOS laptop, 2026-05-04:

| Benchmark | lavaan | rust-backed | Speedup |
| --- | ---: | ---: | ---: |
| direct one-factor `sem()` fit | `0.031630 s` | `0.000828 s` | `38.20x` |
| end-to-end `commonfactor()` smoke | `0.164460 s` | `0.064800 s` | `2.54x` |

The larger direct-fit gain is the clean measurement of the first Rust backend
slice. The smaller end-to-end gain is expected because the null-model/CFI path
still delegates to lavaan.

### Next packet

1. Replace the CFI/null-model fallback path by adding a small generic
   covariance-model representation that can handle fixed/free `~~` rows and
   parameter-table refits.
2. Make `resid_rust()` real for supported fit objects.
3. Add `commonfactorGWAS_rust()` once `parTable` mutation plus refit is native.
