# lavaan_fast generic compiler plan

## Why start here

The specialized `lavaanrust` solvers have already shown that replacing lavaan
under unchanged GenomicSEM wrappers can produce large gains. The next risk is
coverage: if every new model family needs its own parser, parameter-table
builder, and solver, the backend will become a collection of bespoke fast paths.

`lavaan_fast` should therefore start as a generic compiler layer. The first
compiler boundary is the lavaan-style parameter table, not raw syntax:

- GenomicSEM already consumes and mutates parameter tables.
- A parameter table carries free/fixed status and equality structure directly.
- Syntax parsing and model compilation are separable problems.
- Compiling parameter tables first lets us reuse the existing fixtures while we
  test the generic numerical representation.

## Initial representation

The first compiler packet lowers a supported parameter table into RAM form:

- `A`: directed paths, including both `~` regressions and `=~` loadings
- `S`: residual covariance matrix from `~~` rows
- `F`: selector from the full variable system to observed variables

The implied observed covariance is:

```text
Sigma = F (I - A)^-1 S (I - A)^-T F^T
```

This representation is broad enough to cover the current one-factor
`userGWAS()` families and is the natural bridge toward multi-factor path models,
direct SNP effects, residual covariances, and equality-constrained parameter
tables.

## First compiler packet

The first compiler packet is deliberately non-invasive:

1. Compile existing fitted parameter tables into a generic RAM intermediate
   representation.
2. Recompute implied covariance matrices from that IR.
3. Recompute analytic Jacobians over `vech(Sigma)` from that IR.
4. Prove exact agreement against the current rust-backed fitted objects for:
   - unrestricted one-factor `userGWAS()`
   - fixed-measurement one-factor `userGWAS()`
5. Keep unsupported operators such as `:=` strict-errors for now.

No existing fit path is switched over in this packet. The point is to establish
the compiler seam before moving solvers onto it.

## Next compiler packets

1. Add a native generic DWLS optimizer over the same RAM IR, so the hot path
   can reuse generic surfaces without recomputing them after a specialized fit.
2. Re-express one current family through that native compiler-backed optimizer.
3. Add parser support for a larger syntax subset:
   - fixed coefficients
   - labels
   - residual covariances
   - direct effects
4. Add parameter-table equality reuse, then decide whether inequality
   constraints and `:=` belong in the same generic layer or a later symbolic
   layer.

The decision point after that is whether `lavaan_fast` remains a compiler that
feeds a few specialized kernels, or becomes the front end for a fully generic
SEM optimizer.
