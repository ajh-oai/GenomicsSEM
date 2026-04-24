#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GenomicControl {
    Conserv,
    Standard,
    None,
}

impl GenomicControl {
    pub fn from_code(code: i32) -> Option<Self> {
        match code {
            0 => Some(Self::Conserv),
            1 => Some(Self::Standard),
            2 => Some(Self::None),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum KernelError {
    BadDimensions,
    BadIndex,
    BadGenomicControl,
    NullPointer,
    ThreadPoolBuild,
    Singular,
}

pub type KernelResult<T> = Result<T, KernelError>;

#[inline]
fn idx(row: usize, col: usize, nrow: usize) -> usize {
    row + col * nrow
}

#[inline]
fn ridx(row: usize, col: usize, ncol: usize) -> usize {
    row * ncol + col
}

fn solve_linear(a: &[f64], b: &[f64], n: usize) -> KernelResult<Vec<f64>> {
    if a.len() != n * n || b.len() != n {
        return Err(KernelError::BadDimensions);
    }

    let mut aug = vec![0.0; n * (n + 1)];
    for row in 0..n {
        for col in 0..n {
            aug[ridx(row, col, n + 1)] = a[ridx(row, col, n)];
        }
        aug[ridx(row, n, n + 1)] = b[row];
    }

    for col in 0..n {
        let mut pivot = col;
        let mut pivot_abs = aug[ridx(col, col, n + 1)].abs();
        for row in (col + 1)..n {
            let value_abs = aug[ridx(row, col, n + 1)].abs();
            if value_abs > pivot_abs {
                pivot = row;
                pivot_abs = value_abs;
            }
        }

        if !pivot_abs.is_finite() || pivot_abs <= f64::EPSILON {
            return Err(KernelError::Singular);
        }

        if pivot != col {
            for j in col..=n {
                aug.swap(ridx(col, j, n + 1), ridx(pivot, j, n + 1));
            }
        }

        for row in (col + 1)..n {
            let factor = aug[ridx(row, col, n + 1)] / aug[ridx(col, col, n + 1)];
            aug[ridx(row, col, n + 1)] = 0.0;
            for j in (col + 1)..=n {
                aug[ridx(row, j, n + 1)] -= factor * aug[ridx(col, j, n + 1)];
            }
        }
    }

    let mut x = vec![0.0; n];
    for row_rev in 0..n {
        let row = n - 1 - row_rev;
        let mut sum = aug[ridx(row, n, n + 1)];
        for col in (row + 1)..n {
            sum -= aug[ridx(row, col, n + 1)] * x[col];
        }
        x[row] = sum / aug[ridx(row, row, n + 1)];
    }

    Ok(x)
}

fn invert_matrix(a: &[f64], n: usize) -> KernelResult<Vec<f64>> {
    if a.len() != n * n {
        return Err(KernelError::BadDimensions);
    }

    let mut inv = vec![0.0; n * n];
    for col in 0..n {
        let mut rhs = vec![0.0; n];
        rhs[col] = 1.0;
        let sol = solve_linear(a, &rhs, n)?;
        for row in 0..n {
            inv[ridx(row, col, n)] = sol[row];
        }
    }
    Ok(inv)
}

fn vech_index(row: usize, col: usize, n: usize) -> usize {
    debug_assert!(row >= col);
    col * n - (col * (col.saturating_sub(1))) / 2 + (row - col)
}

fn commonfactor_implied_delta(k: usize, p: &[f64], sigma: &mut [f64], delta: &mut [f64]) {
    let n = k + 1;
    let q = 2 * k + 2;

    sigma.fill(0.0);
    delta.fill(0.0);

    let b = p[k - 1];
    let psi = p[2 * k];
    let var_x = p[2 * k + 1];
    let phi = psi + b * b * var_x;

    let lambda = |trait_zero: usize, p: &[f64]| -> f64 {
        if trait_zero == 0 {
            1.0
        } else {
            p[trait_zero - 1]
        }
    };

    sigma[vech_index(0, 0, n)] = var_x;
    delta[ridx(vech_index(0, 0, n), 2 * k + 1, q)] = 1.0;

    for j in 0..k {
        let lam_j = lambda(j, p);
        let row = j + 1;
        let out_idx = vech_index(row, 0, n);
        sigma[out_idx] = var_x * b * lam_j;
        if j > 0 {
            delta[ridx(out_idx, j - 1, q)] = var_x * b;
        }
        delta[ridx(out_idx, k - 1, q)] = var_x * lam_j;
        delta[ridx(out_idx, 2 * k + 1, q)] = b * lam_j;
    }

    for col_trait in 0..k {
        let lam_col = lambda(col_trait, p);
        for row_trait in col_trait..k {
            let lam_row = lambda(row_trait, p);
            let out_idx = vech_index(row_trait + 1, col_trait + 1, n);
            sigma[out_idx] = lam_row * lam_col * phi;
            if row_trait == col_trait {
                sigma[out_idx] += p[k + row_trait];
                delta[ridx(out_idx, k + row_trait, q)] = 1.0;
            }

            if row_trait > 0 {
                delta[ridx(out_idx, row_trait - 1, q)] += lam_col * phi;
            }
            if col_trait > 0 && col_trait != row_trait {
                delta[ridx(out_idx, col_trait - 1, q)] += lam_row * phi;
            }

            delta[ridx(out_idx, k - 1, q)] = 2.0 * b * var_x * lam_row * lam_col;
            delta[ridx(out_idx, 2 * k, q)] = lam_row * lam_col;
            delta[ridx(out_idx, 2 * k + 1, q)] = b * b * lam_row * lam_col;
        }
    }
}

fn commonfactor_q_implied_delta(
    k: usize,
    fixed: &[f64],
    p: &[f64],
    sigma: &mut [f64],
    delta: &mut [f64],
) {
    let n = k + 1;
    let q = 2 * k;

    sigma.fill(0.0);
    delta.fill(0.0);

    let b = fixed[k];
    let psi = fixed[k + 1];
    let var_x = fixed[k + 2];

    sigma[vech_index(0, 0, n)] = var_x;

    for j in 0..k {
        let lam_j = fixed[j];
        let gamma_j = p[j];
        let a_j = lam_j * b + gamma_j;
        let out_idx = vech_index(j + 1, 0, n);

        sigma[out_idx] = var_x * a_j;
        delta[ridx(out_idx, j, q)] = var_x;
    }

    for col_trait in 0..k {
        let lam_col = fixed[col_trait];
        let gamma_col = p[col_trait];
        let a_col = lam_col * b + gamma_col;

        for row_trait in col_trait..k {
            let lam_row = fixed[row_trait];
            let gamma_row = p[row_trait];
            let a_row = lam_row * b + gamma_row;
            let out_idx = vech_index(row_trait + 1, col_trait + 1, n);

            sigma[out_idx] = lam_row * lam_col * psi + a_row * a_col * var_x;

            if row_trait == col_trait {
                sigma[out_idx] += p[k + row_trait];
                delta[ridx(out_idx, row_trait, q)] += 2.0 * a_row * var_x;
                delta[ridx(out_idx, k + row_trait, q)] = 1.0;
            } else {
                delta[ridx(out_idx, row_trait, q)] += a_col * var_x;
                delta[ridx(out_idx, col_trait, q)] += a_row * var_x;
            }
        }
    }
}

fn vech_from_matrix_col_major(matrix: &[f64], n: usize, out: &mut [f64]) {
    let mut out_i = 0;
    for col in 0..n {
        for row in col..n {
            out[out_i] = matrix[idx(row, col, n)];
            out_i += 1;
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn fit_commonfactor_main(
    k: usize,
    s_full: &[f64],
    s_nrow: usize,
    s_ncol: usize,
    v_full_reorder: &[f64],
    v_nrow: usize,
    v_ncol: usize,
    w_diag: &[f64],
    start: &[f64],
    max_iter: usize,
    tol: f64,
    out: &mut [f64],
) -> KernelResult<()> {
    let n = k + 1;
    let m = n * (n + 1) / 2;
    let q = 2 * k + 2;
    let out_len = 2 * q + 3;

    if k == 0
        || s_nrow != n
        || s_ncol != n
        || s_full.len() != n * n
        || v_nrow != m
        || v_ncol != m
        || v_full_reorder.len() != m * m
        || w_diag.len() != m
        || start.len() != q
        || out.len() != out_len
    {
        return Err(KernelError::BadDimensions);
    }

    let mut y = vec![0.0; m];
    vech_from_matrix_col_major(s_full, n, &mut y);

    let mut p = start.to_vec();
    let mut sigma = vec![0.0; m];
    let mut delta = vec![0.0; m * q];
    let mut converged = false;
    let mut iterations = 0usize;

    for iter in 0..max_iter {
        iterations = iter + 1;
        commonfactor_implied_delta(k, &p, &mut sigma, &mut delta);

        let mut residual = vec![0.0; m];
        let mut obj = 0.0;
        for r in 0..m {
            residual[r] = y[r] - sigma[r];
            obj += w_diag[r] * residual[r] * residual[r];
        }

        let mut a = vec![0.0; q * q];
        let mut g = vec![0.0; q];
        for row in 0..m {
            let wr = w_diag[row] * residual[row];
            let w = w_diag[row];
            for col_p in 0..q {
                let d_col = delta[ridx(row, col_p, q)];
                g[col_p] += d_col * wr;
                for row_p in 0..q {
                    a[ridx(row_p, col_p, q)] += delta[ridx(row, row_p, q)] * w * d_col;
                }
            }
        }

        let step = solve_linear(&a, &g, q)?;
        let mut max_scaled_step = 0.0_f64;
        let mut alpha = 1.0_f64;
        let mut accepted = false;
        let mut candidate = p.clone();

        for _ in 0..24 {
            for j in 0..q {
                candidate[j] = p[j] + alpha * step[j];
            }

            let mut sigma_candidate = vec![0.0; m];
            let mut delta_candidate = vec![0.0; m * q];
            commonfactor_implied_delta(k, &candidate, &mut sigma_candidate, &mut delta_candidate);
            let mut candidate_obj = 0.0;
            for row in 0..m {
                let r = y[row] - sigma_candidate[row];
                candidate_obj += w_diag[row] * r * r;
            }

            if candidate_obj.is_finite() && candidate_obj <= obj {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if !accepted {
            alpha = 1.0e-4;
            for j in 0..q {
                candidate[j] = p[j] + alpha * step[j];
            }
        }

        for j in 0..q {
            let scaled = (alpha * step[j]).abs() / (1.0 + p[j].abs());
            if scaled > max_scaled_step {
                max_scaled_step = scaled;
            }
        }
        p.copy_from_slice(&candidate);

        if max_scaled_step < tol {
            converged = true;
            break;
        }
    }

    commonfactor_implied_delta(k, &p, &mut sigma, &mut delta);

    let mut obj = 0.0;
    for row in 0..m {
        let r = y[row] - sigma[row];
        obj += w_diag[row] * r * r;
    }

    let mut a = vec![0.0; q * q];
    for row in 0..m {
        let w = w_diag[row];
        for col_p in 0..q {
            let d_col = delta[ridx(row, col_p, q)];
            for row_p in 0..q {
                a[ridx(row_p, col_p, q)] += delta[ridx(row, row_p, q)] * w * d_col;
            }
        }
    }

    let bread = invert_matrix(&a, q)?;

    let mut meat = vec![0.0; q * q];
    for a_col in 0..q {
        for b_col in 0..q {
            let mut sum = 0.0;
            for r1 in 0..m {
                let wr1 = w_diag[r1] * delta[ridx(r1, a_col, q)];
                for r2 in 0..m {
                    let wr2 = w_diag[r2] * delta[ridx(r2, b_col, q)];
                    sum += wr1 * v_full_reorder[idx(r1, r2, m)] * wr2;
                }
            }
            meat[ridx(a_col, b_col, q)] = sum;
        }
    }

    let mut tmp = vec![0.0; q * q];
    for row in 0..q {
        for col in 0..q {
            let mut sum = 0.0;
            for mid in 0..q {
                sum += bread[ridx(row, mid, q)] * meat[ridx(mid, col, q)];
            }
            tmp[ridx(row, col, q)] = sum;
        }
    }

    let mut cov = vec![0.0; q * q];
    for row in 0..q {
        for col in 0..q {
            let mut sum = 0.0;
            for mid in 0..q {
                sum += tmp[ridx(row, mid, q)] * bread[ridx(mid, col, q)];
            }
            cov[ridx(row, col, q)] = sum;
        }
    }

    out[..q].copy_from_slice(&p);
    for j in 0..q {
        out[q + j] = cov[ridx(j, j, q)].max(0.0).sqrt();
    }
    out[2 * q] = obj;
    out[2 * q + 1] = if converged { 1.0 } else { 0.0 };
    out[2 * q + 2] = iterations as f64;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn fit_commonfactor_q(
    k: usize,
    s_full: &[f64],
    s_nrow: usize,
    s_ncol: usize,
    v_full_reorder: &[f64],
    v_nrow: usize,
    v_ncol: usize,
    w_diag: &[f64],
    fixed: &[f64],
    start: &[f64],
    max_iter: usize,
    tol: f64,
    out: &mut [f64],
) -> KernelResult<()> {
    let n = k + 1;
    let m = n * (n + 1) / 2;
    let q = 2 * k;
    let out_len = q + k * k + 3;

    if k == 0
        || s_nrow != n
        || s_ncol != n
        || s_full.len() != n * n
        || v_nrow != m
        || v_ncol != m
        || v_full_reorder.len() != m * m
        || w_diag.len() != m
        || fixed.len() != k + 3
        || start.len() != q
        || out.len() != out_len
    {
        return Err(KernelError::BadDimensions);
    }

    let mut y = vec![0.0; m];
    vech_from_matrix_col_major(s_full, n, &mut y);

    let mut p = start.to_vec();
    let mut sigma = vec![0.0; m];
    let mut delta = vec![0.0; m * q];
    let mut converged = false;
    let mut iterations = 0usize;

    for iter in 0..max_iter {
        iterations = iter + 1;
        commonfactor_q_implied_delta(k, fixed, &p, &mut sigma, &mut delta);

        let mut residual = vec![0.0; m];
        let mut obj = 0.0;
        for r in 0..m {
            residual[r] = y[r] - sigma[r];
            obj += w_diag[r] * residual[r] * residual[r];
        }

        let mut a = vec![0.0; q * q];
        let mut g = vec![0.0; q];
        for row in 0..m {
            let wr = w_diag[row] * residual[row];
            let w = w_diag[row];
            for col_p in 0..q {
                let d_col = delta[ridx(row, col_p, q)];
                g[col_p] += d_col * wr;
                for row_p in 0..q {
                    a[ridx(row_p, col_p, q)] += delta[ridx(row, row_p, q)] * w * d_col;
                }
            }
        }

        let step = solve_linear(&a, &g, q)?;
        let mut max_scaled_step = 0.0_f64;
        let mut alpha = 1.0_f64;
        let mut accepted = false;
        let mut candidate = p.clone();

        for _ in 0..24 {
            for j in 0..q {
                candidate[j] = p[j] + alpha * step[j];
            }

            let mut sigma_candidate = vec![0.0; m];
            let mut delta_candidate = vec![0.0; m * q];
            commonfactor_q_implied_delta(
                k,
                fixed,
                &candidate,
                &mut sigma_candidate,
                &mut delta_candidate,
            );

            let mut candidate_obj = 0.0;
            for row in 0..m {
                let r = y[row] - sigma_candidate[row];
                candidate_obj += w_diag[row] * r * r;
            }

            if candidate_obj.is_finite() && candidate_obj <= obj {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if !accepted {
            alpha = 1.0e-4;
            for j in 0..q {
                candidate[j] = p[j] + alpha * step[j];
            }
        }

        for j in 0..q {
            let scaled = (alpha * step[j]).abs() / (1.0 + p[j].abs());
            if scaled > max_scaled_step {
                max_scaled_step = scaled;
            }
        }
        p.copy_from_slice(&candidate);

        if max_scaled_step < tol {
            converged = true;
            break;
        }
    }

    commonfactor_q_implied_delta(k, fixed, &p, &mut sigma, &mut delta);

    let mut obj = 0.0;
    for row in 0..m {
        let r = y[row] - sigma[row];
        obj += w_diag[row] * r * r;
    }

    let mut a = vec![0.0; q * q];
    for row in 0..m {
        let w = w_diag[row];
        for col_p in 0..q {
            let d_col = delta[ridx(row, col_p, q)];
            for row_p in 0..q {
                a[ridx(row_p, col_p, q)] += delta[ridx(row, row_p, q)] * w * d_col;
            }
        }
    }

    let bread = invert_matrix(&a, q)?;

    let mut meat = vec![0.0; q * q];
    for a_col in 0..q {
        for b_col in 0..q {
            let mut sum = 0.0;
            for r1 in 0..m {
                let wr1 = w_diag[r1] * delta[ridx(r1, a_col, q)];
                for r2 in 0..m {
                    let wr2 = w_diag[r2] * delta[ridx(r2, b_col, q)];
                    sum += wr1 * v_full_reorder[idx(r1, r2, m)] * wr2;
                }
            }
            meat[ridx(a_col, b_col, q)] = sum;
        }
    }

    let mut tmp = vec![0.0; q * q];
    for row in 0..q {
        for col in 0..q {
            let mut sum = 0.0;
            for mid in 0..q {
                sum += bread[ridx(row, mid, q)] * meat[ridx(mid, col, q)];
            }
            tmp[ridx(row, col, q)] = sum;
        }
    }

    let mut cov = vec![0.0; q * q];
    for row in 0..q {
        for col in 0..q {
            let mut sum = 0.0;
            for mid in 0..q {
                sum += tmp[ridx(row, mid, q)] * bread[ridx(mid, col, q)];
            }
            cov[ridx(row, col, q)] = sum;
        }
    }

    out[..q].copy_from_slice(&p);
    for col in 0..k {
        for row in 0..k {
            out[q + idx(row, col, k)] = cov[ridx(row, col, q)];
        }
    }
    out[q + k * k] = obj;
    out[q + k * k + 1] = if converged { 1.0 } else { 0.0 };
    out[q + k * k + 2] = iterations as f64;

    Ok(())
}

pub fn fill_v_snp(
    se_snp: &[f64],
    se_nrow: usize,
    se_ncol: usize,
    i_zero: usize,
    i_ld: &[f64],
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    var_snp: &[f64],
    coords: &[i32],
    coords_nrow: usize,
    coords_ncol: usize,
    k: usize,
    gc: GenomicControl,
    out: &mut [f64],
) -> KernelResult<()> {
    if se_ncol < k
        || i_ld_nrow < k
        || i_ld_ncol < k
        || coords_ncol < 2
        || out.len() != k * k
        || se_snp.len() != se_nrow * se_ncol
        || i_ld.len() != i_ld_nrow * i_ld_ncol
        || coords.len() != coords_nrow * coords_ncol
    {
        return Err(KernelError::BadDimensions);
    }

    if i_zero >= se_nrow || i_zero >= var_snp.len() {
        return Err(KernelError::BadIndex);
    }

    out.fill(0.0);
    for d in 0..k {
        out[idx(d, d, k)] = 1.0;
    }

    let var = var_snp[i_zero];
    let var_sq = var * var;

    for p in 0..coords_nrow {
        let x_one = coords[idx(p, 0, coords_nrow)];
        let y_one = coords[idx(p, 1, coords_nrow)];
        if x_one <= 0 || y_one <= 0 {
            return Err(KernelError::BadIndex);
        }

        let x = (x_one - 1) as usize;
        let y = (y_one - 1) as usize;
        if x >= k || y >= k {
            return Err(KernelError::BadIndex);
        }

        let se_x = se_snp[idx(i_zero, x, se_nrow)];
        let se_y = se_snp[idx(i_zero, y, se_nrow)];
        let i_xy = i_ld[idx(x, y, i_ld_nrow)];

        let value = if x == y {
            let i_xx = i_ld[idx(x, x, i_ld_nrow)];
            match gc {
                GenomicControl::Conserv => (se_x * i_xx * var).powi(2),
                GenomicControl::Standard => (se_x * i_xx.sqrt() * var).powi(2),
                GenomicControl::None => (se_x * var).powi(2),
            }
        } else {
            match gc {
                GenomicControl::Conserv => {
                    let i_xx = i_ld[idx(x, x, i_ld_nrow)];
                    let i_yy = i_ld[idx(y, y, i_ld_nrow)];
                    se_y * se_x * i_xy * i_xx * i_yy * var_sq
                }
                GenomicControl::Standard => {
                    let i_xx = i_ld[idx(x, x, i_ld_nrow)].sqrt();
                    let i_yy = i_ld[idx(y, y, i_ld_nrow)].sqrt();
                    se_y * se_x * i_xy * i_xx * i_yy * var_sq
                }
                GenomicControl::None => se_y * se_x * i_xy * var_sq,
            }
        };

        out[idx(x, y, k)] = value;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn fill_v_snp_batch(
    se_snp: &[f64],
    se_nrow: usize,
    se_ncol: usize,
    i_ld: &[f64],
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    var_snp: &[f64],
    coords: &[i32],
    coords_nrow: usize,
    coords_ncol: usize,
    k: usize,
    gc: GenomicControl,
    n_threads: usize,
    out: &mut [f64],
) -> KernelResult<()> {
    if out.len() != se_nrow * k * k {
        return Err(KernelError::BadDimensions);
    }

    let run = || {
        use rayon::prelude::*;

        out.par_chunks_mut(k * k)
            .enumerate()
            .try_for_each(|(i_zero, out_one)| {
                fill_v_snp(
                    se_snp,
                    se_nrow,
                    se_ncol,
                    i_zero,
                    i_ld,
                    i_ld_nrow,
                    i_ld_ncol,
                    var_snp,
                    coords,
                    coords_nrow,
                    coords_ncol,
                    k,
                    gc,
                    out_one,
                )
            })
    };

    if n_threads <= 1 {
        for (i_zero, out_one) in out.chunks_mut(k * k).enumerate() {
            fill_v_snp(
                se_snp,
                se_nrow,
                se_ncol,
                i_zero,
                i_ld,
                i_ld_nrow,
                i_ld_ncol,
                var_snp,
                coords,
                coords_nrow,
                coords_ncol,
                k,
                gc,
                out_one,
            )?;
        }
        return Ok(());
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|_| KernelError::ThreadPoolBuild)?
        .install(run)
}

pub fn fill_v_full(
    k: usize,
    v_ld: &[f64],
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    var_snp_se2: f64,
    v_snp: &[f64],
    v_snp_nrow: usize,
    v_snp_ncol: usize,
    out: &mut [f64],
) -> KernelResult<()> {
    let n = (k + 1) * (k + 2) / 2;
    let ld_start = k + 1;

    if out.len() != n * n
        || v_ld_nrow != n - ld_start
        || v_ld_ncol != n - ld_start
        || v_ld.len() != v_ld_nrow * v_ld_ncol
        || v_snp_nrow != k
        || v_snp_ncol != k
        || v_snp.len() != k * k
    {
        return Err(KernelError::BadDimensions);
    }

    out.fill(0.0);
    for d in 0..n {
        out[idx(d, d, n)] = 1.0;
    }

    out[idx(0, 0, n)] = var_snp_se2;

    for col in 0..k {
        for row in 0..k {
            out[idx(row + 1, col + 1, n)] = v_snp[idx(row, col, k)];
        }
    }

    for col in 0..v_ld_ncol {
        for row in 0..v_ld_nrow {
            out[idx(row + ld_start, col + ld_start, n)] = v_ld[idx(row, col, v_ld_nrow)];
        }
    }

    Ok(())
}

pub fn fill_s_full(
    k: usize,
    s_ld: &[f64],
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    var_snp: &[f64],
    beta_snp: &[f64],
    beta_nrow: usize,
    beta_ncol: usize,
    i_zero: usize,
    out: &mut [f64],
) -> KernelResult<()> {
    let n = k + 1;

    if s_ld_nrow != k
        || s_ld_ncol != k
        || s_ld.len() != k * k
        || beta_ncol < k
        || beta_snp.len() != beta_nrow * beta_ncol
        || out.len() != n * n
    {
        return Err(KernelError::BadDimensions);
    }

    if i_zero >= beta_nrow || i_zero >= var_snp.len() {
        return Err(KernelError::BadIndex);
    }

    out.fill(0.0);
    for d in 0..n {
        out[idx(d, d, n)] = 1.0;
    }

    for col in 0..k {
        for row in 0..k {
            out[idx(row + 1, col + 1, n)] = s_ld[idx(row, col, k)];
        }
    }

    let var = var_snp[i_zero];
    out[idx(0, 0, n)] = var;
    for p in 0..k {
        let cov = var * beta_snp[idx(i_zero, p, beta_nrow)];
        out[idx(p + 1, 0, n)] = cov;
        out[idx(0, p + 1, n)] = cov;
    }

    Ok(())
}

pub fn fill_z_pre(
    beta_snp: &[f64],
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: &[f64],
    se_nrow: usize,
    se_ncol: usize,
    i_ld: &[f64],
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    i_zero: usize,
    gc: GenomicControl,
    out: &mut [f64],
) -> KernelResult<()> {
    let k = out.len();
    if beta_ncol < k
        || se_ncol < k
        || i_ld_nrow < k
        || i_ld_ncol < k
        || beta_snp.len() != beta_nrow * beta_ncol
        || se_snp.len() != se_nrow * se_ncol
        || i_ld.len() != i_ld_nrow * i_ld_ncol
    {
        return Err(KernelError::BadDimensions);
    }

    if i_zero >= beta_nrow || i_zero >= se_nrow {
        return Err(KernelError::BadIndex);
    }

    for p in 0..k {
        let beta = beta_snp[idx(i_zero, p, beta_nrow)];
        let se = se_snp[idx(i_zero, p, se_nrow)];
        let diag = i_ld[idx(p, p, i_ld_nrow)];
        out[p] = match gc {
            GenomicControl::Conserv => beta / (se * diag),
            GenomicControl::Standard => beta / (se * diag.sqrt()),
            GenomicControl::None => beta / se,
        };
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fills_v_snp_in_r_column_major_order() {
        let se_snp = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let i_ld = vec![1.2, 0.4, 0.4, 1.5];
        let var_snp = vec![0.7, 0.8, 0.9];
        let coords = vec![1, 2, 1, 2, 1, 1, 2, 2];
        let mut out = vec![0.0; 4];

        fill_v_snp(
            &se_snp,
            3,
            2,
            1,
            &i_ld,
            2,
            2,
            &var_snp,
            &coords,
            4,
            2,
            2,
            GenomicControl::Standard,
            &mut out,
        )
        .unwrap();

        assert_eq!(out[0], (0.2_f64 * 1.2_f64.sqrt() * 0.8_f64).powi(2));
        assert_eq!(out[3], (0.5_f64 * 1.5_f64.sqrt() * 0.8_f64).powi(2));
        assert_eq!(
            out[1],
            0.5_f64 * 0.2_f64 * 0.4_f64 * 1.2_f64.sqrt() * 1.5_f64.sqrt() * 0.8_f64.powi(2)
        );
    }

    #[test]
    fn fills_s_full_with_snp_row_and_column() {
        let s_ld = vec![1.0, 0.2, 0.2, 1.0];
        let var_snp = vec![0.5];
        let beta_snp = vec![0.1, -0.2];
        let mut out = vec![0.0; 9];

        fill_s_full(2, &s_ld, 2, 2, &var_snp, &beta_snp, 1, 2, 0, &mut out).unwrap();

        assert_eq!(out[idx(0, 0, 3)], 0.5);
        assert_eq!(out[idx(1, 0, 3)], 0.05);
        assert_eq!(out[idx(0, 2, 3)], -0.1);
        assert_eq!(out[idx(2, 2, 3)], 1.0);
    }

    #[test]
    fn batch_matches_single_v_snp() {
        let se_snp = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let i_ld = vec![1.2, 0.4, 0.4, 1.5];
        let var_snp = vec![0.7, 0.8, 0.9];
        let coords = vec![1, 2, 1, 2, 1, 1, 2, 2];
        let mut batch = vec![0.0; 3 * 2 * 2];

        fill_v_snp_batch(
            &se_snp,
            3,
            2,
            &i_ld,
            2,
            2,
            &var_snp,
            &coords,
            4,
            2,
            2,
            GenomicControl::Conserv,
            2,
            &mut batch,
        )
        .unwrap();

        for i_zero in 0..3 {
            let mut one = vec![0.0; 4];
            fill_v_snp(
                &se_snp,
                3,
                2,
                i_zero,
                &i_ld,
                2,
                2,
                &var_snp,
                &coords,
                4,
                2,
                2,
                GenomicControl::Conserv,
                &mut one,
            )
            .unwrap();
            assert_eq!(&batch[(i_zero * 4)..((i_zero + 1) * 4)], one.as_slice());
        }
    }
}
