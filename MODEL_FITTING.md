# Rust Model-Fitting Fast Paths

This note describes the two native model-fitting paths added in this branch. They target different
problems:

- `commonfactorGWAS()` has a fixed generated model shape, so it uses a specialized analytic solver.
- `userGWAS()` accepts user-supplied model syntax, so it uses a compiled generic SEM solver for the
  supported subset of models.

## `commonfactorGWAS()`: Specialized Analytic Solver

`commonfactorGWAS()` is unusually amenable to specialization because the function itself generates
the model rather than accepting arbitrary lavaan syntax. For `k` traits, the generated model is:

```r
F1 =~ T1 + T2 + ... + Tk
F1 ~ SNP
T1 ~ 0*SNP
...
Tk ~ 0*SNP
```

The Rust solver implements that exact model directly instead of routing each SNP through a general
SEM engine.

### Parameterization

Let `X` denote the SNP, `F` the common factor, and `Y_j` the `j`th observed trait:

```text
F   = b * X + zeta
Y_j = lambda_j * F + epsilon_j
lambda_1 = 1
```

The main model has `2k + 2` free parameters:

```text
lambda_2 ... lambda_k
b
theta_1 ... theta_k
psi
var_x
```

where:

- `theta_j = Var(epsilon_j)`
- `psi = Var(zeta)`
- `var_x = Var(X)`

### Implied Covariance and Analytic Jacobian

The model-implied covariance entries have closed forms:

```text
Var(X)        = var_x
Cov(X, Y_j)   = var_x * b * lambda_j
Cov(Y_i, Y_j) = lambda_i * lambda_j * (psi + b^2 * var_x)
Var(Y_j)      = lambda_j^2 * (psi + b^2 * var_x) + theta_j
```

The native path computes both the implied lower-triangle covariance vector and its Jacobian with
respect to the free parameters analytically. This is implemented in
`commonfactor_implied_delta()` in `crates/genomicssem-core/src/lib.rs`.

That is the main difference from the generic `userGWAS()` path: the common-factor solver does not
need to build generic RAM matrices or estimate derivatives numerically.

### Optimization

The solver minimizes the same DWLS objective as the legacy path:

```text
f(p) = sum_r w_r * (s_r - sigma_r(p))^2
```

where:

- `s` is `vech(S_full)` for the current SNP;
- `sigma(p)` is the model-implied covariance vector;
- `w_r` is the diagonal DWLS weight derived from `V_full`.

It uses Gauss-Newton iterations:

```text
A    = Delta' W Delta
g    = Delta' W (s - sigma)
step = solve(A, g)
```

with a short backtracking line search. After convergence, the solver computes the same
sandwich-corrected covariance used by the original GenomicSEM code:

```text
bread = (Delta' W Delta)^-1
meat  = Delta' W V W Delta
cov   = bread * meat * bread
```

The corrected standard error returned for the GWAS result is the standard error of the `b`
parameter, corresponding to `F1 ~ SNP`.

### Starting Values

The fast path uses a domain-specific initializer from the LDSC trait covariance matrix:

1. take the leading eigenvector of `S_LD`;
2. scale it so the first loading is `1`;
3. derive an initial factor variance from the leading eigenvalue;
4. derive residual variances from the diagonal remainder;
5. initialize the SNP-to-factor effect at `0`.

The R helper for this is `.commonfactor_fast_start_from_cov()`.

### Q Statistic

The legacy `commonfactorGWAS()` implementation fits a second model to calculate SNP heterogeneity.
The Rust path preserves that behavior with a second specialized solver.

For the follow-up model, the first-stage common-factor parameters are held fixed and trait-specific
direct SNP effects are added:

```text
Y_j = lambda_j * F + gamma_j * X + epsilon_j
```

The free parameters are:

```text
gamma_1 ... gamma_k
theta_1 ... theta_k
```

The solver again uses analytic implied-covariance formulas and analytic derivatives, implemented in
`commonfactor_q_implied_delta()`. After fitting, it computes:

```text
Q = gamma' Var(gamma)^-1 gamma
```

which matches the role of the follow-up lavaan Q model in the original implementation.

### Batched Execution

For each SNP, the native routine:

1. builds `S_full`, `V_snp`, and `V_full`;
2. fits the main common-factor model;
3. fits the Q follow-up model;
4. returns the factor effect, corrected standard error, and Q statistic.

`fit_commonfactor_batch()` runs that per-SNP routine across the full scan using Rayon worker
threads. The result is fast because it removes repeated lavaan parsing, repeated generic SEM setup,
and R-level per-SNP orchestration while preserving the original workflow outputs.

## `userGWAS()`: Generic Compiled Solver

`userGWAS()` cannot use the same specialization because its model is user-provided. The supported
fast path instead:

1. parses supported lavaan-like rows in R;
2. compiles them into RAM-style `B` and `Psi` matrices plus parameter maps and bounds;
3. evaluates the implied covariance as:

```text
Sigma = (I - B)^-1 Psi (I - B)^-T
```

4. fits the model with the same DWLS objective in Rust;
5. batches the per-SNP scan in Rust.

The generic path trades some per-model efficiency for modeling flexibility. Unsupported syntax still
falls back to the original R/lavaan path.
