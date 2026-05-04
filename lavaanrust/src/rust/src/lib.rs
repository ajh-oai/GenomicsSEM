use extendr_api::prelude::*;
use nalgebra::{DMatrix, DVector, SymmetricEigen};

fn from_col_major(values: &[f64], nrows: usize, ncols: usize) -> DMatrix<f64> {
    DMatrix::from_column_slice(nrows, ncols, values)
}

fn to_col_major(matrix: &DMatrix<f64>) -> Vec<f64> {
    matrix.as_slice().to_vec()
}

fn vech(matrix: &DMatrix<f64>) -> DVector<f64> {
    let n = matrix.nrows();
    let mut values = Vec::with_capacity(n * (n + 1) / 2);

    for col in 0..n {
        for row in col..n {
            values.push(matrix[(row, col)]);
        }
    }

    DVector::from_vec(values)
}

fn implied_one_factor(loadings: &DVector<f64>, residuals: &DVector<f64>) -> DMatrix<f64> {
    let mut implied = loadings * loadings.transpose();

    for idx in 0..residuals.len() {
        implied[(idx, idx)] += residuals[idx];
    }

    implied
}

fn delta_one_factor(loadings: &DVector<f64>) -> DMatrix<f64> {
    let n = loadings.len();
    let n_stats = n * (n + 1) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, 2 * n);
    let mut row = 0;

    for col in 0..n {
        for obs_row in col..n {
            if obs_row == col {
                delta[(row, col)] = 2.0 * loadings[col];
                delta[(row, n + col)] = 1.0;
            } else {
                delta[(row, col)] = loadings[obs_row];
                delta[(row, obs_row)] = loadings[col];
            }

            row += 1;
        }
    }

    delta
}

fn initial_one_factor(sample_cov: &DMatrix<f64>) -> (DVector<f64>, DVector<f64>) {
    let eigen = SymmetricEigen::new(sample_cov.clone());
    let mut best_idx = 0;

    for idx in 1..eigen.eigenvalues.len() {
        if eigen.eigenvalues[idx] > eigen.eigenvalues[best_idx] {
            best_idx = idx;
        }
    }

    let scale = eigen.eigenvalues[best_idx].max(1e-8).sqrt();
    let mut loadings = eigen.eigenvectors.column(best_idx).into_owned() * scale;

    if loadings.sum() < 0.0 {
        loadings *= -1.0;
    }

    let mut residuals = DVector::<f64>::zeros(sample_cov.nrows());
    for idx in 0..sample_cov.nrows() {
        residuals[idx] = (sample_cov[(idx, idx)] - loadings[idx] * loadings[idx]).max(1e-8);
    }

    (loadings, residuals)
}

fn srmr(sample_cov: &DMatrix<f64>, implied: &DMatrix<f64>) -> f64 {
    let observed_sd: Vec<f64> = (0..sample_cov.nrows())
        .map(|idx| sample_cov[(idx, idx)].abs().sqrt())
        .collect();
    let mut sum_sq = 0.0;
    let mut count = 0usize;

    for col in 0..sample_cov.ncols() {
        for row in col..sample_cov.nrows() {
            let denom = observed_sd[row] * observed_sd[col];
            if denom > 0.0 {
                let resid = (sample_cov[(row, col)] - implied[(row, col)]) / denom;
                sum_sq += resid * resid;
                count += 1;
            }
        }
    }

    if count == 0 {
        0.0
    } else {
        (sum_sq / count as f64).sqrt()
    }
}

fn from_vech(values: &DVector<f64>, n: usize) -> DMatrix<f64> {
    let mut matrix = DMatrix::<f64>::zeros(n, n);
    let mut idx = 0;

    for col in 0..n {
        for row in col..n {
            matrix[(row, col)] = values[idx];
            matrix[(col, row)] = values[idx];
            idx += 1;
        }
    }

    matrix
}

/// Fit the covariance-only one-factor DWLS slice used by GenomicSEM's
/// `commonfactor()` model.
///
/// `sample_cov` and `wls_v` are flattened column-major R matrices.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param n Number of observed variables.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_one_factor_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    n: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if n <= 0 {
        return Err(Error::Other("n must be positive".into()));
    }

    let n = n as usize;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = from_col_major(&wls_values, n_stats, n_stats);
    let observed = vech(&sample_cov);
    let (mut loadings, mut residuals) = initial_one_factor(&sample_cov);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_one_factor(&loadings, &residuals);
        let residual = &observed - vech(&implied);
        let delta = delta_one_factor(&loadings);
        let gradient = delta.transpose() * &wls_v * &residual;
        let hessian = delta.transpose() * &wls_v * &delta;
        let mut regularized = hessian.clone();

        for idx in 0..regularized.nrows() {
            regularized[(idx, idx)] += damping;
        }

        let Some(step) = regularized.lu().solve(&gradient) else {
            damping *= 10.0;
            continue;
        };

        let old_objective = 0.5 * residual.dot(&(&wls_v * &residual));
        let mut next_loadings = loadings.clone();
        let mut next_residuals = residuals.clone();

        for idx in 0..n {
            next_loadings[idx] += step[idx];
            next_residuals[idx] = (next_residuals[idx] + step[n + idx]).max(1e-10);
        }

        let next_implied = implied_one_factor(&next_loadings, &next_residuals);
        let next_residual = &observed - vech(&next_implied);
        let next_objective = 0.5 * next_residual.dot(&(&wls_v * &next_residual));

        if next_objective <= old_objective {
            loadings = next_loadings;
            residuals = next_residuals;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    if loadings.sum() < 0.0 {
        loadings *= -1.0;
    }

    let implied = implied_one_factor(&loadings, &residuals);
    let residual = &observed - vech(&implied);
    let delta = delta_one_factor(&loadings);
    let bread = (delta.transpose() * &wls_v * &delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(2 * n, 2 * n));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = 0.5 * residual.dot(&(&wls_v * &residual));

    Ok(list!(
        loadings = loadings.as_slice().to_vec(),
        residuals = residuals.as_slice().to_vec(),
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Fit a covariance-only observed-variable DWLS model where the free parameters
/// are selected entries of `vech(Sigma)`.
///
/// This covers the common-factor null model and the parameter-table refit used
/// to estimate its residual-covariance Q statistic.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param free_mask Integer mask over `vech(Sigma)`.
/// @param fixed_values Fixed values over `vech(Sigma)`.
/// @param n Number of observed variables.
/// @export
#[extendr]
fn fit_observed_covariance_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    free_mask: Integers,
    fixed_values: Doubles,
    n: i32,
) -> std::result::Result<List, Error> {
    if n <= 0 {
        return Err(Error::Other("n must be positive".into()));
    }

    let n = n as usize;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();
    let free_mask_values = free_mask.iter().map(|value| value.0 != 0).collect::<Vec<_>>();
    let fixed_values = fixed_values.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    if free_mask_values.len() != n_stats || fixed_values.len() != n_stats {
        return Err(Error::Other(
            "free_mask and fixed_values must match vech(sample_cov)".into(),
        ));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = from_col_major(&wls_values, n_stats, n_stats);
    let observed = vech(&sample_cov);
    let fixed = DVector::from_vec(fixed_values);
    let free_indices = free_mask_values
        .iter()
        .enumerate()
        .filter_map(|(idx, is_free)| is_free.then_some(idx))
        .collect::<Vec<_>>();
    let mut delta = DMatrix::<f64>::zeros(n_stats, free_indices.len());

    for (col, row) in free_indices.iter().enumerate() {
        delta[(*row, col)] = 1.0;
    }

    let target = observed.clone() - fixed.clone();
    let hessian = delta.transpose() * &wls_v * &delta;
    let rhs = delta.transpose() * &wls_v * &target;
    let Some(params) = hessian.clone().lu().solve(&rhs) else {
        return Err(Error::Other("observed covariance solve failed".into()));
    };

    let implied_vec = fixed + &delta * &params;
    let implied = from_vech(&implied_vec, n);
    let residual = observed - implied_vec;
    let bread = hessian
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(free_indices.len(), free_indices.len()));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = 0.5 * residual.dot(&(&wls_v * &residual));

    Ok(list!(
        estimates = params.as_slice().to_vec(),
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = true
    ))
}

extendr_module! {
    mod lavaanrust;
    fn fit_one_factor_dwls;
    fn fit_observed_covariance_dwls;
}
