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

fn implied_commonfactor_gwas(
    loadings: &DVector<f64>,
    gamma: f64,
    residuals: &DVector<f64>,
    psi: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let mut implied = DMatrix::<f64>::zeros(k + 1, k + 1);
    let factor_total = gamma * gamma * phi + psi;

    implied[(0, 0)] = phi;

    for idx in 0..k {
        let cov = phi * loadings[idx] * gamma;
        implied[(idx + 1, 0)] = cov;
        implied[(0, idx + 1)] = cov;
    }

    for col in 0..k {
        for row in col..k {
            let mut cov = loadings[row] * loadings[col] * factor_total;
            if row == col {
                cov += residuals[row];
            }

            implied[(row + 1, col + 1)] = cov;
            implied[(col + 1, row + 1)] = cov;
        }
    }

    implied
}

fn delta_commonfactor_gwas(
    loadings: &DVector<f64>,
    gamma: f64,
    psi: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let n_stats = (k + 1) * (k + 2) / 2;
    let n_params = 2 * k + 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, n_params);
    let factor_total = gamma * gamma * phi + psi;
    let gamma_col = k - 1;
    let theta_start = gamma_col + 1;
    let psi_col = theta_start + k;
    let phi_col = psi_col + 1;
    let mut stat_row = 0;

    for col in 0..=k {
        for row in col..=k {
            if col == 0 && row == 0 {
                delta[(stat_row, phi_col)] = 1.0;
            } else if col == 0 || row == 0 {
                let trait_idx = row.max(col) - 1;

                if trait_idx > 0 {
                    delta[(stat_row, trait_idx - 1)] = phi * gamma;
                }

                delta[(stat_row, gamma_col)] = phi * loadings[trait_idx];
                delta[(stat_row, phi_col)] = loadings[trait_idx] * gamma;
            } else {
                let trait_row = row - 1;
                let trait_col = col - 1;

                if trait_row > 0 {
                    delta[(stat_row, trait_row - 1)] += loadings[trait_col] * factor_total;
                }

                if trait_col > 0 {
                    delta[(stat_row, trait_col - 1)] += loadings[trait_row] * factor_total;
                }

                delta[(stat_row, gamma_col)] =
                    loadings[trait_row] * loadings[trait_col] * 2.0 * gamma * phi;
                delta[(stat_row, psi_col)] = loadings[trait_row] * loadings[trait_col];
                delta[(stat_row, phi_col)] =
                    loadings[trait_row] * loadings[trait_col] * gamma * gamma;

                if trait_row == trait_col {
                    delta[(stat_row, theta_start + trait_row)] = 1.0;
                }
            }

            stat_row += 1;
        }
    }

    delta
}

fn implied_commonfactor_gwas_q(
    loadings: &DVector<f64>,
    gamma: f64,
    direct: &DVector<f64>,
    residuals: &DVector<f64>,
    psi: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let mut implied = DMatrix::<f64>::zeros(k + 1, k + 1);
    let beta = loadings * gamma + direct;

    implied[(0, 0)] = phi;

    for idx in 0..k {
        let cov = phi * beta[idx];
        implied[(idx + 1, 0)] = cov;
        implied[(0, idx + 1)] = cov;
    }

    for col in 0..k {
        for row in col..k {
            let mut cov = loadings[row] * loadings[col] * psi + phi * beta[row] * beta[col];
            if row == col {
                cov += residuals[row];
            }

            implied[(row + 1, col + 1)] = cov;
            implied[(col + 1, row + 1)] = cov;
        }
    }

    implied
}

fn delta_commonfactor_gwas_q(
    loadings: &DVector<f64>,
    gamma: f64,
    direct: &DVector<f64>,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let n_stats = (k + 1) * (k + 2) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, 2 * k);
    let beta = loadings * gamma + direct;
    let mut stat_row = 0;

    for col in 0..=k {
        for row in col..=k {
            if col == 0 && row == 0 {
                // The predictor variance is fixed in the Q-model refit.
            } else if col == 0 || row == 0 {
                let trait_idx = row.max(col) - 1;
                delta[(stat_row, trait_idx)] = phi;
            } else {
                let trait_row = row - 1;
                let trait_col = col - 1;

                delta[(stat_row, trait_row)] += phi * beta[trait_col];
                delta[(stat_row, trait_col)] += phi * beta[trait_row];

                if trait_row == trait_col {
                    delta[(stat_row, k + trait_row)] = 1.0;
                }
            }

            stat_row += 1;
        }
    }

    delta
}

fn delta_user_gwas_fixed_measurement(
    loadings: &DVector<f64>,
    gamma: f64,
    phi: f64,
) -> DMatrix<f64> {
    let k = loadings.len();
    let n_stats = (k + 1) * (k + 2) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, k + 3);
    let psi_col = k;
    let gamma_col = k + 1;
    let phi_col = k + 2;
    let mut stat_row = 0;

    for col in 0..=k {
        for row in col..=k {
            if col == 0 && row == 0 {
                delta[(stat_row, phi_col)] = 1.0;
            } else if col == 0 || row == 0 {
                let trait_idx = row.max(col) - 1;
                delta[(stat_row, gamma_col)] = phi * loadings[trait_idx];
                delta[(stat_row, phi_col)] = loadings[trait_idx] * gamma;
            } else {
                let trait_row = row - 1;
                let trait_col = col - 1;

                if trait_row == trait_col {
                    delta[(stat_row, trait_row)] = 1.0;
                }

                delta[(stat_row, psi_col)] = loadings[trait_row] * loadings[trait_col];
                delta[(stat_row, gamma_col)] =
                    loadings[trait_row] * loadings[trait_col] * 2.0 * gamma * phi;
                delta[(stat_row, phi_col)] =
                    loadings[trait_row] * loadings[trait_col] * gamma * gamma;
            }

            stat_row += 1;
        }
    }

    delta
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

fn initial_commonfactor_gwas(
    sample_cov: &DMatrix<f64>,
) -> (DVector<f64>, f64, DVector<f64>, f64, f64) {
    let k = sample_cov.nrows() - 1;
    let trait_cov = sample_cov.view((1, 1), (k, k)).into_owned();
    let (std_loadings, residuals) = initial_one_factor(&trait_cov);
    let marker = if std_loadings[0].abs() < 1e-8 {
        1.0
    } else {
        std_loadings[0]
    };
    let mut loadings = std_loadings / marker;
    loadings[0] = 1.0;
    let phi = sample_cov[(0, 0)].max(1e-8);
    let gamma = sample_cov[(1, 0)] / phi;
    let factor_total = marker * marker;
    let psi = (factor_total - gamma * gamma * phi).max(1e-8);

    (loadings, gamma, residuals, psi, phi)
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

fn select_observed(matrix: &DMatrix<f64>, observed_indices: &[usize]) -> DMatrix<f64> {
    let n_observed = observed_indices.len();
    let mut selected = DMatrix::<f64>::zeros(n_observed, n_observed);

    for (out_row, source_row) in observed_indices.iter().enumerate() {
        for (out_col, source_col) in observed_indices.iter().enumerate() {
            selected[(out_row, out_col)] = matrix[(*source_row, *source_col)];
        }
    }

    selected
}

fn implied_ram(
    directed: &DMatrix<f64>,
    covariance: &DMatrix<f64>,
    observed_indices: &[usize],
) -> std::result::Result<(DMatrix<f64>, DMatrix<f64>), Error> {
    let identity = DMatrix::<f64>::identity(directed.nrows(), directed.ncols());
    let Some(inverse) = (identity - directed).try_inverse() else {
        return Err(Error::Other("RAM solve failed because I - A is singular".into()));
    };
    let full_implied = &inverse * covariance * inverse.transpose();

    Ok((select_observed(&full_implied, observed_indices), inverse))
}

fn validate_ram_indices(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    fixed_values: &[f64],
    observed_index: &[i32],
    n_variables: usize,
    n_free: usize,
) -> std::result::Result<Vec<usize>, Error> {
    let n_rows = lhs_index.len();
    if rhs_index.len() != n_rows
        || op_code.len() != n_rows
        || free_index.len() != n_rows
        || fixed_values.len() != n_rows
    {
        return Err(Error::Other(
            "RAM row vectors must all have the same length".into(),
        ));
    }

    if observed_index.is_empty() {
        return Err(Error::Other("observed_index must be non-empty".into()));
    }

    for row_idx in 0..n_rows {
        let lhs = lhs_index[row_idx];
        let rhs = rhs_index[row_idx];
        if lhs <= 0 || rhs <= 0 || lhs as usize > n_variables || rhs as usize > n_variables {
            return Err(Error::Other("RAM row index is out of bounds".into()));
        }

        let free = free_index[row_idx];
        if free < 0 || free as usize > n_free {
            return Err(Error::Other("free_index is out of bounds".into()));
        }

        if !matches!(op_code[row_idx], 1..=3) {
            return Err(Error::Other("unsupported RAM operator code".into()));
        }
    }

    observed_index
        .iter()
        .map(|index| {
            if *index <= 0 || *index as usize > n_variables {
                Err(Error::Other("observed_index is out of bounds".into()))
            } else {
                Ok((*index as usize) - 1)
            }
        })
        .collect::<std::result::Result<Vec<_>, _>>()
}

fn ram_matrices_from_rows(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    fixed_values: &[f64],
    free_values: &[f64],
    n_variables: usize,
) -> (DMatrix<f64>, DMatrix<f64>) {
    let mut directed = DMatrix::<f64>::zeros(n_variables, n_variables);
    let mut covariance = DMatrix::<f64>::zeros(n_variables, n_variables);

    for row_idx in 0..lhs_index.len() {
        let lhs = lhs_index[row_idx] as usize - 1;
        let rhs = rhs_index[row_idx] as usize - 1;
        let free = free_index[row_idx];
        let value = if free > 0 {
            free_values[free as usize - 1]
        } else {
            fixed_values[row_idx]
        };

        match op_code[row_idx] {
            1 => directed[(rhs, lhs)] = value,
            2 => directed[(lhs, rhs)] = value,
            3 => {
                covariance[(lhs, rhs)] = value;
                covariance[(rhs, lhs)] = value;
            }
            _ => unreachable!(),
        }
    }

    (directed, covariance)
}

fn ram_delta_from_rows(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    covariance: &DMatrix<f64>,
    inverse: &DMatrix<f64>,
    observed_indices: &[usize],
    n_variables: usize,
    n_free: usize,
) -> DMatrix<f64> {
    let n_observed = observed_indices.len();
    let n_stats = n_observed * (n_observed + 1) / 2;
    let mut delta = DMatrix::<f64>::zeros(n_stats, n_free);

    for free_position in 1..=n_free {
        let mut d_directed = DMatrix::<f64>::zeros(n_variables, n_variables);
        let mut d_covariance = DMatrix::<f64>::zeros(n_variables, n_variables);

        for row_idx in 0..lhs_index.len() {
            if free_index[row_idx] != free_position as i32 {
                continue;
            }

            let lhs = lhs_index[row_idx] as usize - 1;
            let rhs = rhs_index[row_idx] as usize - 1;

            match op_code[row_idx] {
                1 => d_directed[(rhs, lhs)] += 1.0,
                2 => d_directed[(lhs, rhs)] += 1.0,
                3 => {
                    d_covariance[(lhs, rhs)] += 1.0;
                    if lhs != rhs {
                        d_covariance[(rhs, lhs)] += 1.0;
                    }
                }
                _ => unreachable!(),
            }
        }

        let d_inverse = inverse * d_directed * inverse;
        let d_full_implied = &d_inverse * covariance * inverse.transpose()
            + inverse * d_covariance * inverse.transpose()
            + inverse * covariance * d_inverse.transpose();
        let d_implied = select_observed(&d_full_implied, observed_indices);
        let d_vech = vech(&d_implied);

        for row_idx in 0..n_stats {
            delta[(row_idx, free_position - 1)] = d_vech[row_idx];
        }
    }

    delta
}

fn ram_surfaces_from_rows(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    fixed_values: &[f64],
    free_values: &[f64],
    observed_indices: &[usize],
    n_variables: usize,
) -> std::result::Result<(DMatrix<f64>, DMatrix<f64>), Error> {
    let (directed, covariance) = ram_matrices_from_rows(
        lhs_index,
        rhs_index,
        op_code,
        free_index,
        fixed_values,
        free_values,
        n_variables,
    );
    let (implied, inverse) = implied_ram(&directed, &covariance, observed_indices)?;
    let delta = ram_delta_from_rows(
        lhs_index,
        rhs_index,
        op_code,
        free_index,
        &covariance,
        &inverse,
        observed_indices,
        n_variables,
        free_values.len(),
    );

    Ok((implied, delta))
}

fn covariance_variance_free_mask(
    lhs_index: &[i32],
    rhs_index: &[i32],
    op_code: &[i32],
    free_index: &[i32],
    n_free: usize,
) -> Vec<bool> {
    let mut mask = vec![false; n_free];

    for row_idx in 0..lhs_index.len() {
        let free = free_index[row_idx];
        if op_code[row_idx] == 3 && lhs_index[row_idx] == rhs_index[row_idx] && free > 0 {
            mask[free as usize - 1] = true;
        }
    }

    mask
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
    let free_mask_values = free_mask
        .iter()
        .map(|value| value.0 != 0)
        .collect::<Vec<_>>();
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

/// Fit the marker-scaled common-factor GWAS DWLS slice used by
/// GenomicSEM's `commonfactorGWAS()` model.
///
/// Variable order is expected to be `SNP, trait_1, ..., trait_k`.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param k Number of trait indicators.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_commonfactor_gwas_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    k: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if k <= 0 {
        return Err(Error::Other("k must be positive".into()));
    }

    let k = k as usize;
    let n = k + 1;
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
    let (mut loadings, mut gamma, mut residuals, mut psi, mut phi) =
        initial_commonfactor_gwas(&sample_cov);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
        let residual = &observed - vech(&implied);
        let delta = delta_commonfactor_gwas(&loadings, gamma, psi, phi);
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

        for idx in 1..k {
            next_loadings[idx] += step[idx - 1];
        }

        let gamma_idx = k - 1;
        let theta_start = gamma_idx + 1;
        let next_gamma = gamma + step[gamma_idx];

        for idx in 0..k {
            next_residuals[idx] = (next_residuals[idx] + step[theta_start + idx]).max(1e-10);
        }

        let next_psi = (psi + step[theta_start + k]).max(1e-10);
        let next_phi = (phi + step[theta_start + k + 1]).max(1e-10);
        let next_implied = implied_commonfactor_gwas(
            &next_loadings,
            next_gamma,
            &next_residuals,
            next_psi,
            next_phi,
        );
        let next_residual = &observed - vech(&next_implied);
        let next_objective = 0.5 * next_residual.dot(&(&wls_v * &next_residual));

        if next_objective <= old_objective {
            loadings = next_loadings;
            gamma = next_gamma;
            residuals = next_residuals;
            psi = next_psi;
            phi = next_phi;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
    let residual = &observed - vech(&implied);
    let delta = delta_commonfactor_gwas(&loadings, gamma, psi, phi);
    let bread = (delta.transpose() * &wls_v * &delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(2 * k + 2, 2 * k + 2));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = 0.5 * residual.dot(&(&wls_v * &residual));

    Ok(list!(
        loadings = loadings.as_slice().to_vec(),
        gamma = gamma,
        residuals = residuals.as_slice().to_vec(),
        psi = psi,
        phi = phi,
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Fit the common-factor GWAS Q-model refit where direct SNP effects and trait
/// residual variances are free and all first-stage quantities are fixed.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param loadings Fixed marker-scaled trait loadings.
/// @param gamma Fixed factor regression coefficient.
/// @param direct Initial direct SNP effects.
/// @param residuals Initial trait residual variances.
/// @param psi Fixed factor residual variance.
/// @param phi Fixed SNP variance.
/// @param k Number of trait indicators.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_commonfactor_gwas_q_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    loadings: Doubles,
    gamma: f64,
    direct: Doubles,
    residuals: Doubles,
    psi: f64,
    phi: f64,
    k: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if k <= 0 {
        return Err(Error::Other("k must be positive".into()));
    }

    let k = k as usize;
    let n = k + 1;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();
    let loadings = loadings.iter().map(|value| value.0).collect::<Vec<_>>();
    let direct = direct.iter().map(|value| value.0).collect::<Vec<_>>();
    let residuals = residuals.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    if loadings.len() != k || direct.len() != k || residuals.len() != k {
        return Err(Error::Other(
            "loadings, direct, and residuals must each have length k".into(),
        ));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = from_col_major(&wls_values, n_stats, n_stats);
    let observed = vech(&sample_cov);
    let loadings = DVector::from_vec(loadings);
    let mut direct = DVector::from_vec(direct);
    let mut residuals = DVector::from_vec(residuals);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_commonfactor_gwas_q(&loadings, gamma, &direct, &residuals, psi, phi);
        let residual = &observed - vech(&implied);
        let delta = delta_commonfactor_gwas_q(&loadings, gamma, &direct, phi);
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
        let mut next_direct = direct.clone();
        let mut next_residuals = residuals.clone();

        for idx in 0..k {
            next_direct[idx] += step[idx];
            next_residuals[idx] = (next_residuals[idx] + step[k + idx]).max(1e-10);
        }

        let next_implied =
            implied_commonfactor_gwas_q(&loadings, gamma, &next_direct, &next_residuals, psi, phi);
        let next_residual = &observed - vech(&next_implied);
        let next_objective = 0.5 * next_residual.dot(&(&wls_v * &next_residual));

        if next_objective <= old_objective {
            direct = next_direct;
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

    let implied = implied_commonfactor_gwas_q(&loadings, gamma, &direct, &residuals, psi, phi);
    let residual = &observed - vech(&implied);
    let delta = delta_commonfactor_gwas_q(&loadings, gamma, &direct, phi);
    let bread = (delta.transpose() * &wls_v * &delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(2 * k, 2 * k));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = 0.5 * residual.dot(&(&wls_v * &residual));

    Ok(list!(
        direct = direct.as_slice().to_vec(),
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

/// Fit the fixed-measurement one-factor GWAS DWLS slice used by the default
/// `userGWAS(..., fix_measurement = TRUE)` path.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param loadings Fixed marker-scaled trait loadings.
/// @param residuals Initial trait residual variances.
/// @param psi Initial factor residual variance.
/// @param gamma Initial factor regression coefficient.
/// @param phi Initial SNP variance.
/// @param k Number of trait indicators.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_user_gwas_fixed_measurement_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    loadings: Doubles,
    residuals: Doubles,
    psi: f64,
    gamma: f64,
    phi: f64,
    k: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if k <= 0 {
        return Err(Error::Other("k must be positive".into()));
    }

    let k = k as usize;
    let n = k + 1;
    let n_stats = n * (n + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();
    let loadings = loadings.iter().map(|value| value.0).collect::<Vec<_>>();
    let residuals = residuals.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n * n {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    if loadings.len() != k || residuals.len() != k {
        return Err(Error::Other(
            "loadings and residuals must each have length k".into(),
        ));
    }

    let sample_cov = from_col_major(&sample_cov_values, n, n);
    let wls_v = from_col_major(&wls_values, n_stats, n_stats);
    let observed = vech(&sample_cov);
    let loadings = DVector::from_vec(loadings);
    let mut residuals = DVector::from_vec(residuals);
    let mut psi = psi.max(1e-10);
    let mut gamma = gamma;
    let mut phi = phi.max(1e-10);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;

        let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
        let residual = &observed - vech(&implied);
        let delta = delta_user_gwas_fixed_measurement(&loadings, gamma, phi);
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
        let mut next_residuals = residuals.clone();

        for idx in 0..k {
            next_residuals[idx] = (next_residuals[idx] + step[idx]).max(1e-10);
        }

        let next_psi = (psi + step[k]).max(1e-10);
        let next_gamma = gamma + step[k + 1];
        let next_phi = (phi + step[k + 2]).max(1e-10);
        let next_implied = implied_commonfactor_gwas(
            &loadings,
            next_gamma,
            &next_residuals,
            next_psi,
            next_phi,
        );
        let next_residual = &observed - vech(&next_implied);
        let next_objective = 0.5 * next_residual.dot(&(&wls_v * &next_residual));

        if next_objective <= old_objective {
            residuals = next_residuals;
            psi = next_psi;
            gamma = next_gamma;
            phi = next_phi;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    let implied = implied_commonfactor_gwas(&loadings, gamma, &residuals, psi, phi);
    let residual = &observed - vech(&implied);
    let delta = delta_user_gwas_fixed_measurement(&loadings, gamma, phi);
    let bread = (delta.transpose() * &wls_v * &delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(k + 3, k + 3));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = 0.5 * residual.dot(&(&wls_v * &residual));

    Ok(list!(
        residuals = residuals.as_slice().to_vec(),
        psi = psi,
        gamma = gamma,
        phi = phi,
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

/// Evaluate generic RAM-model implied covariance and Jacobian surfaces.
///
/// Row indices are 1-based to match R. Operator codes are:
/// `1 = =~`, `2 = ~`, `3 = ~~`.
/// @param lhs_index Parameter-table lhs row indices.
/// @param rhs_index Parameter-table rhs row indices.
/// @param op_code Parameter-table operator codes.
/// @param free_index Free-parameter indices, with `0` for fixed rows.
/// @param fixed_values Fixed row values, ignored for free rows.
/// @param free_values Current free-parameter values.
/// @param observed_index Observed-variable indices in the full RAM system.
/// @param n_variables Number of variables in the full RAM system.
/// @export
#[extendr]
fn evaluate_ram_surfaces(
    lhs_index: Integers,
    rhs_index: Integers,
    op_code: Integers,
    free_index: Integers,
    fixed_values: Doubles,
    free_values: Doubles,
    observed_index: Integers,
    n_variables: i32,
) -> std::result::Result<List, Error> {
    if n_variables <= 0 {
        return Err(Error::Other("n_variables must be positive".into()));
    }

    let n_variables = n_variables as usize;
    let lhs_index = lhs_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let rhs_index = rhs_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let op_code = op_code.iter().map(|value| value.0).collect::<Vec<_>>();
    let free_index = free_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let fixed_values = fixed_values.iter().map(|value| value.0).collect::<Vec<_>>();
    let free_values = free_values.iter().map(|value| value.0).collect::<Vec<_>>();
    let observed_index = observed_index
        .iter()
        .map(|value| value.0)
        .collect::<Vec<_>>();

    let observed_indices = validate_ram_indices(
        &lhs_index,
        &rhs_index,
        &op_code,
        &free_index,
        &fixed_values,
        &observed_index,
        n_variables,
        free_values.len(),
    )?;
    let (implied, delta) = ram_surfaces_from_rows(
        &lhs_index,
        &rhs_index,
        &op_code,
        &free_index,
        &fixed_values,
        &free_values,
        &observed_indices,
        n_variables,
    )?;

    Ok(list!(
        implied = to_col_major(&implied),
        delta = to_col_major(&delta)
    ))
}

/// Fit a generic RAM-model DWLS slice from compiled row data.
///
/// This is the generic optimizer counterpart to `evaluate_ram_surfaces()`.
/// Free diagonal covariance parameters are constrained to remain positive.
/// @param sample_cov Flattened observed covariance matrix.
/// @param wls_v Flattened DWLS weight matrix.
/// @param lhs_index Parameter-table lhs row indices.
/// @param rhs_index Parameter-table rhs row indices.
/// @param op_code Parameter-table operator codes.
/// @param free_index Free-parameter indices, with `0` for fixed rows.
/// @param fixed_values Fixed row values, ignored for free rows.
/// @param free_values Initial free-parameter values.
/// @param observed_index Observed-variable indices in the full RAM system.
/// @param n_variables Number of variables in the full RAM system.
/// @param max_iter Maximum optimizer iterations.
/// @param tol Convergence tolerance.
/// @export
#[extendr]
fn fit_ram_dwls(
    sample_cov: Doubles,
    wls_v: Doubles,
    lhs_index: Integers,
    rhs_index: Integers,
    op_code: Integers,
    free_index: Integers,
    fixed_values: Doubles,
    free_values: Doubles,
    observed_index: Integers,
    n_variables: i32,
    max_iter: i32,
    tol: f64,
) -> std::result::Result<List, Error> {
    if n_variables <= 0 {
        return Err(Error::Other("n_variables must be positive".into()));
    }

    let n_variables = n_variables as usize;
    let lhs_index = lhs_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let rhs_index = rhs_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let op_code = op_code.iter().map(|value| value.0).collect::<Vec<_>>();
    let free_index = free_index.iter().map(|value| value.0).collect::<Vec<_>>();
    let fixed_values = fixed_values.iter().map(|value| value.0).collect::<Vec<_>>();
    let initial_free_values = free_values.iter().map(|value| value.0).collect::<Vec<_>>();
    let observed_index = observed_index
        .iter()
        .map(|value| value.0)
        .collect::<Vec<_>>();
    let observed_indices = validate_ram_indices(
        &lhs_index,
        &rhs_index,
        &op_code,
        &free_index,
        &fixed_values,
        &observed_index,
        n_variables,
        initial_free_values.len(),
    )?;

    let n_observed = observed_indices.len();
    let n_stats = n_observed * (n_observed + 1) / 2;
    let sample_cov_values = sample_cov.iter().map(|value| value.0).collect::<Vec<_>>();
    let wls_values = wls_v.iter().map(|value| value.0).collect::<Vec<_>>();

    if sample_cov_values.len() != n_observed * n_observed {
        return Err(Error::Other("sample_cov has the wrong size".into()));
    }

    if wls_values.len() != n_stats * n_stats {
        return Err(Error::Other("wls_v has the wrong size".into()));
    }

    let sample_cov = from_col_major(&sample_cov_values, n_observed, n_observed);
    let wls_v = from_col_major(&wls_values, n_stats, n_stats);
    let observed = vech(&sample_cov);
    let variance_mask = covariance_variance_free_mask(
        &lhs_index,
        &rhs_index,
        &op_code,
        &free_index,
        initial_free_values.len(),
    );
    let mut params = DVector::from_vec(initial_free_values);
    let mut damping = 1e-6;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..max_iter.max(1) {
        iterations = iter + 1;
        let params_vec = params.as_slice();
        let (implied, delta) = ram_surfaces_from_rows(
            &lhs_index,
            &rhs_index,
            &op_code,
            &free_index,
            &fixed_values,
            params_vec,
            &observed_indices,
            n_variables,
        )?;
        let residual = &observed - vech(&implied);
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
        let mut next_params = &params + &step;

        for (idx, should_clamp) in variance_mask.iter().enumerate() {
            if *should_clamp {
                next_params[idx] = next_params[idx].max(1e-10);
            }
        }

        let (next_implied, _) = ram_surfaces_from_rows(
            &lhs_index,
            &rhs_index,
            &op_code,
            &free_index,
            &fixed_values,
            next_params.as_slice(),
            &observed_indices,
            n_variables,
        )?;
        let next_residual = &observed - vech(&next_implied);
        let next_objective = 0.5 * next_residual.dot(&(&wls_v * &next_residual));

        if next_objective <= old_objective {
            params = next_params;
            damping = (damping / 3.0).max(1e-12);

            if step.amax() < tol {
                converged = true;
                break;
            }
        } else {
            damping *= 10.0;
        }
    }

    let (implied, delta) = ram_surfaces_from_rows(
        &lhs_index,
        &rhs_index,
        &op_code,
        &free_index,
        &fixed_values,
        params.as_slice(),
        &observed_indices,
        n_variables,
    )?;
    let residual = &observed - vech(&implied);
    let bread = (delta.transpose() * &wls_v * &delta)
        .try_inverse()
        .unwrap_or_else(|| DMatrix::<f64>::zeros(params.len(), params.len()));
    let naive_se = bread.diagonal().map(|value| value.max(0.0).sqrt());
    let objective = 0.5 * residual.dot(&(&wls_v * &residual));

    Ok(list!(
        estimates = params.as_slice().to_vec(),
        implied = to_col_major(&implied),
        delta = to_col_major(&delta),
        naive_se = naive_se.as_slice().to_vec(),
        objective = objective,
        srmr = srmr(&sample_cov, &implied),
        converged = converged,
        iterations = iterations
    ))
}

extendr_module! {
    mod lavaanrust;
    fn fit_one_factor_dwls;
    fn fit_observed_covariance_dwls;
    fn fit_commonfactor_gwas_dwls;
    fn fit_commonfactor_gwas_q_dwls;
    fn fit_user_gwas_fixed_measurement_dwls;
    fn evaluate_ram_surfaces;
    fn fit_ram_dwls;
}
