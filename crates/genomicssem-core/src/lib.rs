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

pub const COMMONFACTOR_BATCH_OUT_COLS: usize = 7;
pub const MUNGE_QC_COUNT_LEN: usize = 8;
pub const SUMSTATS_QC_COUNT_LEN: usize = 10;

mod prep;
pub use prep::{
    munge_fused, munge_fused_batch, sumstats_fused, sumstats_fused_batch, MungeFusedInput,
    MungeFusedReport, SumstatsFusedBatchOutput, SumstatsFusedInput, SumstatsFusedOutput,
    SumstatsFusedReport, MUNGE_FUSED_COLS, SUMSTATS_FUSED_COLS,
};

#[inline]
fn idx(row: usize, col: usize, nrow: usize) -> usize {
    row + col * nrow
}

#[inline]
fn ridx(row: usize, col: usize, ncol: usize) -> usize {
    row * ncol + col
}

fn median_finite(values: &mut Vec<f64>) -> f64 {
    values.retain(|x| x.is_finite());
    if values.is_empty() {
        return f64::NAN;
    }

    values.sort_by(|a, b| a.total_cmp(b));
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        0.5 * (values[mid - 1] + values[mid])
    }
}

fn normal_quantile(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const P_LOW: f64 = 0.02425;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }

    if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        return (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0);
    }

    if p <= P_HIGH {
        let q = p - 0.5;
        let r = q * q;
        return (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0);
    }

    let q = (-2.0 * (1.0 - p).ln()).sqrt();
    -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
        / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
}

#[inline]
fn p_to_z_abs(p: f64) -> f64 {
    -normal_quantile(0.5 * p)
}

#[inline]
fn valid_allele(x: i32) -> bool {
    (1..=4).contains(&x)
}

#[inline]
fn allele_flip(a1_ref: i32, a1_file: i32, a2_file: i32) -> bool {
    a1_ref != a1_file && a1_ref == a2_file
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

fn quadratic_form_inverse(a: &[f64], x: &[f64], n: usize) -> KernelResult<f64> {
    if a.len() != n * n || x.len() != n {
        return Err(KernelError::BadDimensions);
    }

    let sol = solve_linear(a, x, n)?;
    let mut out = 0.0;
    for i in 0..n {
        out += x[i] * sol[i];
    }
    Ok(out)
}

pub struct MungeQcOutput<'a> {
    pub keep: &'a mut [i32],
    pub z: &'a mut [f64],
    pub counts: &'a mut [i32],
}

pub fn munge_qc(
    a1_ref: &[i32],
    a2_ref: &[i32],
    a1_file: &[i32],
    a2_file: &[i32],
    effect: &[f64],
    p: &[f64],
    info: Option<&[f64]>,
    maf: Option<&[f64]>,
    info_filter: f64,
    maf_filter: f64,
    out: &mut MungeQcOutput<'_>,
) -> KernelResult<usize> {
    let n = effect.len();
    if a1_ref.len() != n
        || a2_ref.len() != n
        || a1_file.len() != n
        || a2_file.len() != n
        || p.len() != n
        || info.is_some_and(|x| x.len() != n)
        || maf.is_some_and(|x| x.len() != n)
        || out.keep.len() < n
        || out.z.len() < n
        || out.counts.len() < MUNGE_QC_COUNT_LEN
    {
        return Err(KernelError::BadDimensions);
    }

    out.counts[..MUNGE_QC_COUNT_LEN].fill(0);

    let mut median_values = Vec::with_capacity(n);
    for i in 0..n {
        if !p[i].is_finite() {
            continue;
        }
        if !effect[i].is_finite() {
            continue;
        }
        median_values.push(effect[i]);
    }
    let or_detected = median_finite(&mut median_values).round() == 1.0;
    out.counts[7] = i32::from(or_detected);

    let mut out_n = 0usize;
    for i in 0..n {
        let p_i = p[i];
        if !p_i.is_finite() {
            out.counts[0] += 1;
            continue;
        }

        let mut effect_i = effect[i];
        if !effect_i.is_finite() {
            out.counts[1] += 1;
            continue;
        }

        if or_detected {
            effect_i = effect_i.ln();
        }

        let a1r = a1_ref[i];
        let a2r = a2_ref[i];
        let a1f = a1_file[i];
        let a2f = a2_file[i];
        if !valid_allele(a1r) || !valid_allele(a2r) || !valid_allele(a1f) || !valid_allele(a2f) {
            out.counts[2] += 1;
            continue;
        }

        if allele_flip(a1r, a1f, a2f) {
            effect_i = -effect_i;
        }

        if a1r != a1f && a1r != a2f {
            out.counts[2] += 1;
            continue;
        }
        if a2r != a2f && a2r != a1f {
            out.counts[3] += 1;
            continue;
        }

        if p_i > 1.0 || p_i < 0.0 {
            out.counts[4] += 1;
        }

        if let Some(info) = info {
            if !info[i].is_finite() || info[i] < info_filter {
                out.counts[5] += 1;
                continue;
            }
        }

        if let Some(maf) = maf {
            if !maf[i].is_finite() || maf[i] < maf_filter {
                out.counts[6] += 1;
                continue;
            }
        }

        out.keep[out_n] = (i + 1) as i32;
        out.z[out_n] = effect_i.signum() * p_to_z_abs(p_i);
        out_n += 1;
    }

    Ok(out_n)
}

pub struct SumstatsQcOutput<'a> {
    pub keep: &'a mut [i32],
    pub beta: &'a mut [f64],
    pub se_out: &'a mut [f64],
    pub counts: &'a mut [i32],
}

#[allow(clippy::too_many_arguments)]
pub fn sumstats_qc(
    a1_ref: &[i32],
    a2_ref: &[i32],
    a1_file: &[i32],
    a2_file: &[i32],
    effect: &[f64],
    se: &[f64],
    p: &[f64],
    n_values: &[f64],
    maf_ref: &[f64],
    maf_file: Option<&[f64]>,
    info: Option<&[f64]>,
    info_filter: f64,
    ols: bool,
    beta_is_character: bool,
    linprob: bool,
    se_logit: bool,
    out: &mut SumstatsQcOutput<'_>,
) -> KernelResult<usize> {
    let n = effect.len();
    if a1_ref.len() != n
        || a2_ref.len() != n
        || a1_file.len() != n
        || a2_file.len() != n
        || se.len() != n
        || p.len() != n
        || n_values.len() != n
        || maf_ref.len() != n
        || maf_file.is_some_and(|x| x.len() != n)
        || info.is_some_and(|x| x.len() != n)
        || out.keep.len() < n
        || out.beta.len() < n
        || out.se_out.len() < n
        || out.counts.len() < SUMSTATS_QC_COUNT_LEN
    {
        return Err(KernelError::BadDimensions);
    }

    out.counts[..SUMSTATS_QC_COUNT_LEN].fill(0);

    let mut median_values = Vec::with_capacity(n);
    for i in 0..n {
        if !p[i].is_finite() || !effect[i].is_finite() {
            continue;
        }
        let maf = maf_file.map(|x| x[i]).unwrap_or(maf_ref[i]);
        let maf = if maf > 0.5 { 1.0 - maf } else { maf };
        if maf_file.is_some() && (maf == 0.0 || maf == 1.0 || !maf.is_finite()) {
            continue;
        }
        median_values.push(effect[i]);
    }
    let or_detected = median_finite(&mut median_values).round() == 1.0;
    out.counts[9] = i32::from(or_detected);

    let mut out_n = 0usize;
    for i in 0..n {
        let p_i = p[i];
        if !p_i.is_finite() {
            out.counts[0] += 1;
            continue;
        }

        let mut effect_i = effect[i];
        if !effect_i.is_finite() {
            out.counts[1] += 1;
            continue;
        }

        let maf_raw = maf_file.map(|x| x[i]).unwrap_or(maf_ref[i]);
        let maf = if maf_raw > 0.5 {
            1.0 - maf_raw
        } else {
            maf_raw
        };
        if maf_file.is_some() && (maf == 0.0 || maf == 1.0 || !maf.is_finite()) {
            out.counts[2] += 1;
            continue;
        }
        let var_snp = 2.0 * maf * (1.0 - maf);

        if or_detected {
            effect_i = effect_i.ln();
        }

        if effect_i == 0.0 || !effect_i.is_finite() {
            out.counts[3] += 1;
            continue;
        }

        let z = effect_i.signum() * p_to_z_abs(p_i);
        let mut se_i = se[i];
        if ols && !beta_is_character {
            if !n_values[i].is_finite() || !var_snp.is_finite() || var_snp <= 0.0 {
                continue;
            }
            effect_i = z / (n_values[i] * var_snp).sqrt();
        }
        if linprob {
            if !n_values[i].is_finite() || !var_snp.is_finite() || var_snp <= 0.0 {
                continue;
            }
            effect_i = z / ((n_values[i] / 4.0) * var_snp).sqrt();
            se_i = 1.0 / ((n_values[i] / 4.0) * var_snp).sqrt();
        }

        let a1r = a1_ref[i];
        let a2r = a2_ref[i];
        let a1f = a1_file[i];
        let a2f = a2_file[i];
        if !valid_allele(a1r) || !valid_allele(a2r) || !valid_allele(a1f) || !valid_allele(a2f) {
            out.counts[4] += 1;
            continue;
        }
        if allele_flip(a1r, a1f, a2f) {
            effect_i = -effect_i;
        }
        if a1r != a1f && a1r != a2f {
            out.counts[4] += 1;
            continue;
        }
        if a2r != a2f && a2r != a1f {
            out.counts[5] += 1;
            continue;
        }

        if p_i > 1.0 || p_i < 0.0 {
            out.counts[6] += 1;
        }

        if let Some(info) = info {
            if !info[i].is_finite() || info[i] < info_filter {
                out.counts[7] += 1;
                continue;
            }
        }

        let (beta_out, se_out) = if ols {
            (effect_i, (effect_i / z).abs())
        } else if linprob {
            let denom =
                ((effect_i * effect_i) * var_snp + std::f64::consts::PI.powi(2) / 3.0).sqrt();
            (effect_i / denom, se_i / denom)
        } else if se_logit {
            let denom =
                ((effect_i * effect_i) * var_snp + std::f64::consts::PI.powi(2) / 3.0).sqrt();
            (effect_i / denom, se_i / denom)
        } else {
            let denom =
                ((effect_i * effect_i) * var_snp + std::f64::consts::PI.powi(2) / 3.0).sqrt();
            (effect_i / denom, (se_i / effect_i.exp()) / denom)
        };

        if !beta_out.is_finite() || !se_out.is_finite() {
            continue;
        }
        if linprob && (beta_out == 0.0 || se_out == 0.0) {
            out.counts[8] += 1;
            continue;
        }

        out.keep[out_n] = (i + 1) as i32;
        out.beta[out_n] = beta_out;
        out.se_out[out_n] = se_out;
        out_n += 1;
    }

    Ok(out_n)
}

fn ldsc_block_bounds(n_snps: usize, n_blocks: usize, block: usize) -> Option<(usize, usize)> {
    let from_one =
        (1.0 + (block as f64) * ((n_snps as f64) - 1.0) / (n_blocks as f64)).floor() as usize;
    let to_one = if block + 1 == n_blocks {
        n_snps
    } else {
        (1.0 + ((block + 1) as f64) * ((n_snps as f64) - 1.0) / (n_blocks as f64)).floor() as usize
            - 1
    };

    if from_one == 0 || to_one < from_one {
        return None;
    }

    Some((from_one - 1, to_one))
}

fn compute_ldsc_block_products(
    weighted_ld: &[f64],
    n_snps: usize,
    n_annot: usize,
    weighted_chi: &[f64],
    start: usize,
    end: usize,
    xty: &mut [f64],
    xtx: &mut [f64],
) {
    xty.fill(0.0);
    xtx.fill(0.0);

    for col in 0..n_annot {
        let x_col = &weighted_ld[idx(start, col, n_snps)..idx(end, col, n_snps)];

        let mut xy = 0.0;
        for row_offset in 0..x_col.len() {
            xy += x_col[row_offset] * weighted_chi[start + row_offset];
        }
        xty[col] = xy;

        for col2 in col..n_annot {
            let x_col2 = &weighted_ld[idx(start, col2, n_snps)..idx(end, col2, n_snps)];
            let mut xx = 0.0;
            for row_offset in 0..x_col.len() {
                xx += x_col[row_offset] * x_col2[row_offset];
            }
            xtx[idx(col, col2, n_annot)] = xx;
            if col2 != col {
                xtx[idx(col2, col, n_annot)] = xx;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn ldsc_block_products(
    weighted_ld: &[f64],
    n_snps: usize,
    n_annot: usize,
    weighted_chi: &[f64],
    n_blocks: usize,
    xty_block: &mut [f64],
    xtx_block: &mut [f64],
    xty: &mut [f64],
    xtx: &mut [f64],
    n_threads: usize,
) -> KernelResult<()> {
    if weighted_ld.len() != n_snps * n_annot
        || weighted_chi.len() != n_snps
        || xty_block.len() != n_blocks * n_annot
        || xtx_block.len() != n_blocks * n_annot * n_annot
        || xty.len() != n_annot
        || xtx.len() != n_annot * n_annot
        || n_blocks == 0
        || n_snps == 0
        || n_annot == 0
    {
        return Err(KernelError::BadDimensions);
    }

    xty_block.fill(0.0);
    xtx_block.fill(0.0);
    xty.fill(0.0);
    xtx.fill(0.0);

    if n_threads <= 1 {
        let mut block_xty = vec![0.0; n_annot];
        let mut block_xtx = vec![0.0; n_annot * n_annot];
        for block in 0..n_blocks {
            let Some((start, end)) = ldsc_block_bounds(n_snps, n_blocks, block) else {
                continue;
            };
            compute_ldsc_block_products(
                weighted_ld,
                n_snps,
                n_annot,
                weighted_chi,
                start,
                end,
                &mut block_xty,
                &mut block_xtx,
            );
            copy_ldsc_block(
                block, n_blocks, n_annot, &block_xty, &block_xtx, xty_block, xtx_block, xty, xtx,
            );
        }
        return Ok(());
    }

    let run = || {
        use rayon::prelude::*;

        (0..n_blocks)
            .into_par_iter()
            .map(|block| {
                let mut block_xty = vec![0.0; n_annot];
                let mut block_xtx = vec![0.0; n_annot * n_annot];
                if let Some((start, end)) = ldsc_block_bounds(n_snps, n_blocks, block) {
                    compute_ldsc_block_products(
                        weighted_ld,
                        n_snps,
                        n_annot,
                        weighted_chi,
                        start,
                        end,
                        &mut block_xty,
                        &mut block_xtx,
                    );
                }
                (block, block_xty, block_xtx)
            })
            .collect::<Vec<_>>()
    };

    let blocks = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|_| KernelError::ThreadPoolBuild)?
        .install(run);

    for (block, block_xty, block_xtx) in blocks {
        copy_ldsc_block(
            block, n_blocks, n_annot, &block_xty, &block_xtx, xty_block, xtx_block, xty, xtx,
        );
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn copy_ldsc_block(
    block: usize,
    n_blocks: usize,
    n_annot: usize,
    block_xty: &[f64],
    block_xtx: &[f64],
    xty_block: &mut [f64],
    xtx_block: &mut [f64],
    xty: &mut [f64],
    xtx: &mut [f64],
) {
    for col in 0..n_annot {
        let value = block_xty[col];
        xty_block[idx(block, col, n_blocks)] = value;
        xty[col] += value;
    }

    for col2 in 0..n_annot {
        for col in 0..n_annot {
            let value = block_xtx[idx(col, col2, n_annot)];
            xtx_block[idx(block * n_annot + col, col2, n_blocks * n_annot)] = value;
            xtx[idx(col, col2, n_annot)] += value;
        }
    }
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

fn validate_parameter_bounds(q: usize, lower: &[f64], upper: &[f64]) -> KernelResult<()> {
    if lower.len() != q || upper.len() != q {
        return Err(KernelError::BadDimensions);
    }

    for j in 0..q {
        if lower[j].is_nan() || upper[j].is_nan() || lower[j] >= upper[j] {
            return Err(KernelError::BadDimensions);
        }
    }

    Ok(())
}

fn project_parameters_to_bounds(p: &mut [f64], lower: &[f64], upper: &[f64]) -> KernelResult<()> {
    validate_parameter_bounds(p.len(), lower, upper)?;

    for j in 0..p.len() {
        if p[j].is_nan() {
            return Err(KernelError::BadDimensions);
        }
        if p[j] < lower[j] {
            p[j] = lower[j];
        } else if p[j] > upper[j] {
            p[j] = upper[j];
        }
        if !p[j].is_finite() {
            return Err(KernelError::BadDimensions);
        }
    }

    Ok(())
}

fn generic_sem_implied(
    obs_n: usize,
    total_n: usize,
    b_fixed_col_major: &[f64],
    psi_fixed_col_major: &[f64],
    b_free_col_major: &[i32],
    psi_free_col_major: &[i32],
    p: &[f64],
    sigma_vech: &mut [f64],
    sigma_obs_col_major: Option<&mut [f64]>,
) -> KernelResult<()> {
    let q = p.len();
    let total_sq = total_n * total_n;
    if obs_n == 0
        || obs_n > total_n
        || b_fixed_col_major.len() != total_sq
        || psi_fixed_col_major.len() != total_sq
        || b_free_col_major.len() != total_sq
        || psi_free_col_major.len() != total_sq
        || sigma_vech.len() != obs_n * (obs_n + 1) / 2
    {
        return Err(KernelError::BadDimensions);
    }

    let mut b = vec![0.0; total_sq];
    let mut psi = vec![0.0; total_sq];
    for col in 0..total_n {
        for row in 0..total_n {
            let source = idx(row, col, total_n);
            let dest = ridx(row, col, total_n);

            let b_free = b_free_col_major[source];
            b[dest] = if b_free > 0 {
                let p_idx = (b_free - 1) as usize;
                if p_idx >= q {
                    return Err(KernelError::BadIndex);
                }
                p[p_idx]
            } else {
                b_fixed_col_major[source]
            };

            let psi_free = psi_free_col_major[source];
            psi[dest] = if psi_free > 0 {
                let p_idx = (psi_free - 1) as usize;
                if p_idx >= q {
                    return Err(KernelError::BadIndex);
                }
                p[p_idx]
            } else {
                psi_fixed_col_major[source]
            };
        }
    }

    let mut a = vec![0.0; total_sq];
    for row in 0..total_n {
        for col in 0..total_n {
            a[ridx(row, col, total_n)] =
                if row == col { 1.0 } else { 0.0 } - b[ridx(row, col, total_n)];
        }
    }

    let inv_a = invert_matrix(&a, total_n)?;

    let mut tmp = vec![0.0; total_sq];
    for row in 0..total_n {
        for col in 0..total_n {
            let mut sum = 0.0;
            for mid in 0..total_n {
                sum += inv_a[ridx(row, mid, total_n)] * psi[ridx(mid, col, total_n)];
            }
            tmp[ridx(row, col, total_n)] = sum;
        }
    }

    let mut cov = vec![0.0; total_sq];
    for row in 0..total_n {
        for col in 0..total_n {
            let mut sum = 0.0;
            for mid in 0..total_n {
                sum += tmp[ridx(row, mid, total_n)] * inv_a[ridx(col, mid, total_n)];
            }
            cov[ridx(row, col, total_n)] = sum;
        }
    }

    let mut out_i = 0;
    for col in 0..obs_n {
        for row in col..obs_n {
            sigma_vech[out_i] = cov[ridx(row, col, total_n)];
            out_i += 1;
        }
    }

    if let Some(out) = sigma_obs_col_major {
        if out.len() != obs_n * obs_n {
            return Err(KernelError::BadDimensions);
        }
        for col in 0..obs_n {
            for row in 0..obs_n {
                out[idx(row, col, obs_n)] = cov[ridx(row, col, total_n)];
            }
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn generic_sem_delta_numeric(
    obs_n: usize,
    total_n: usize,
    b_fixed: &[f64],
    psi_fixed: &[f64],
    b_free: &[i32],
    psi_free: &[i32],
    p: &[f64],
    sigma: &[f64],
    delta: &mut [f64],
) -> KernelResult<()> {
    let m = obs_n * (obs_n + 1) / 2;
    let q = p.len();
    if sigma.len() != m || delta.len() != m * q {
        return Err(KernelError::BadDimensions);
    }

    let mut p_step = p.to_vec();
    let mut sigma_step = vec![0.0; m];
    for param in 0..q {
        let step = 1.0e-6 * (1.0 + p[param].abs());
        p_step[param] = p[param] + step;
        generic_sem_implied(
            obs_n,
            total_n,
            b_fixed,
            psi_fixed,
            b_free,
            psi_free,
            &p_step,
            &mut sigma_step,
            None,
        )?;

        for row in 0..m {
            delta[ridx(row, param, q)] = (sigma_step[row] - sigma[row]) / step;
        }
        p_step[param] = p[param];
    }

    Ok(())
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

#[allow(clippy::too_many_arguments)]
fn fit_commonfactor_snp_fast(
    k: usize,
    s_ld: &[f64],
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    v_ld: &[f64],
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    i_ld: &[f64],
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    beta_snp: &[f64],
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: &[f64],
    se_nrow: usize,
    se_ncol: usize,
    var_snp: &[f64],
    coords: &[i32],
    coords_nrow: usize,
    coords_ncol: usize,
    var_snp_se2: f64,
    gc: GenomicControl,
    start: &[f64],
    max_iter_main: usize,
    max_iter_q: usize,
    tol: f64,
    i_zero: usize,
    out: &mut [f64],
) -> KernelResult<()> {
    if out.len() != COMMONFACTOR_BATCH_OUT_COLS {
        return Err(KernelError::BadDimensions);
    }

    let n = k + 1;
    let m = n * (n + 1) / 2;
    let main_q = 2 * k + 2;
    let q_q = 2 * k;

    let mut s_full = vec![0.0; n * n];
    fill_s_full(
        k,
        s_ld,
        s_ld_nrow,
        s_ld_ncol,
        var_snp,
        beta_snp,
        beta_nrow,
        beta_ncol,
        i_zero,
        &mut s_full,
    )?;

    let mut v_snp = vec![0.0; k * k];
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
        &mut v_snp,
    )?;

    let mut v_full = vec![0.0; m * m];
    fill_v_full(
        k,
        v_ld,
        v_ld_nrow,
        v_ld_ncol,
        var_snp_se2,
        &v_snp,
        k,
        k,
        &mut v_full,
    )?;

    let mut w_diag = vec![0.0; m];
    for d in 0..m {
        let value = v_full[idx(d, d, m)];
        if !value.is_finite() || value == 0.0 {
            return Err(KernelError::Singular);
        }
        w_diag[d] = 1.0 / value;
    }

    let mut main_out = vec![0.0; 2 * main_q + 3];
    fit_commonfactor_main(
        k,
        &s_full,
        n,
        n,
        &v_full,
        m,
        m,
        &w_diag,
        start,
        max_iter_main,
        tol,
        &mut main_out,
    )?;

    let main_converged = main_out[2 * main_q + 1] == 1.0;
    if !main_converged || !main_out[..(2 * main_q)].iter().all(|v| v.is_finite()) {
        return Err(KernelError::Singular);
    }

    let mut fixed = vec![0.0; k + 3];
    fixed[0] = 1.0;
    for j in 1..k {
        fixed[j] = main_out[j - 1];
    }
    fixed[k] = main_out[k - 1];
    fixed[k + 1] = main_out[2 * k];
    fixed[k + 2] = main_out[2 * k + 1];

    let b = fixed[k];
    let psi = fixed[k + 1];
    let var_x = fixed[k + 2];

    let mut q_start = vec![0.0; q_q];
    let mut q_start_finite = var_x.is_finite() && var_x != 0.0;
    for j in 0..k {
        let snp_trait_cov = s_full[idx(0, j + 1, n)];
        q_start[j] = (snp_trait_cov - var_x * b * fixed[j]) / var_x;
        q_start[k + j] = s_full[idx(j + 1, j + 1, n)]
            - fixed[j] * fixed[j] * psi
            - (fixed[j] * b + q_start[j]).powi(2) * var_x;
        q_start_finite &= q_start[j].is_finite() && q_start[k + j].is_finite();
    }
    if !q_start_finite {
        for j in 0..k {
            q_start[j] = 0.0;
            q_start[k + j] = main_out[k + j];
        }
    }

    let mut q_out = vec![0.0; q_q + k * k + 3];
    fit_commonfactor_q(
        k, &s_full, n, n, &v_full, m, m, &w_diag, &fixed, &q_start, max_iter_q, tol, &mut q_out,
    )?;

    let q_converged = q_out[q_q + k * k + 1] == 1.0;
    if !q_converged || !q_out[..(q_q + k * k)].iter().all(|v| v.is_finite()) {
        return Err(KernelError::Singular);
    }

    let mut gamma_cov = vec![0.0; k * k];
    for col in 0..k {
        for row in 0..k {
            gamma_cov[ridx(row, col, k)] = q_out[q_q + idx(row, col, k)];
        }
    }

    let gamma = &q_out[..k];
    let q_stat = quadratic_form_inverse(&gamma_cov, gamma, k)?;
    if !q_stat.is_finite() {
        return Err(KernelError::Singular);
    }

    out[0] = main_out[k - 1];
    out[1] = main_out[main_q + k - 1];
    out[2] = q_stat;
    out[3] = 1.0;
    out[4] = 1.0;
    out[5] = main_out[2 * main_q + 2];
    out[6] = q_out[q_q + k * k + 2];

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn fit_commonfactor_batch(
    k: usize,
    s_ld: &[f64],
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    v_ld: &[f64],
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    i_ld: &[f64],
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    beta_snp: &[f64],
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: &[f64],
    se_nrow: usize,
    se_ncol: usize,
    var_snp: &[f64],
    coords: &[i32],
    coords_nrow: usize,
    coords_ncol: usize,
    var_snp_se2: f64,
    gc: GenomicControl,
    start: &[f64],
    max_iter_main: usize,
    max_iter_q: usize,
    tol: f64,
    n_threads: usize,
    out: &mut [f64],
) -> KernelResult<()> {
    if k == 0
        || beta_nrow != se_nrow
        || beta_ncol < k
        || se_ncol < k
        || var_snp.len() != beta_nrow
        || s_ld.len() != s_ld_nrow * s_ld_ncol
        || v_ld.len() != v_ld_nrow * v_ld_ncol
        || i_ld.len() != i_ld_nrow * i_ld_ncol
        || beta_snp.len() != beta_nrow * beta_ncol
        || se_snp.len() != se_nrow * se_ncol
        || coords.len() != coords_nrow * coords_ncol
        || start.len() != 2 * k + 2
        || out.len() != beta_nrow * COMMONFACTOR_BATCH_OUT_COLS
    {
        return Err(KernelError::BadDimensions);
    }

    let run = || {
        use rayon::prelude::*;

        out.par_chunks_mut(COMMONFACTOR_BATCH_OUT_COLS)
            .enumerate()
            .for_each(|(i_zero, out_one)| {
                let result = fit_commonfactor_snp_fast(
                    k,
                    s_ld,
                    s_ld_nrow,
                    s_ld_ncol,
                    v_ld,
                    v_ld_nrow,
                    v_ld_ncol,
                    i_ld,
                    i_ld_nrow,
                    i_ld_ncol,
                    beta_snp,
                    beta_nrow,
                    beta_ncol,
                    se_snp,
                    se_nrow,
                    se_ncol,
                    var_snp,
                    coords,
                    coords_nrow,
                    coords_ncol,
                    var_snp_se2,
                    gc,
                    start,
                    max_iter_main,
                    max_iter_q,
                    tol,
                    i_zero,
                    out_one,
                );
                if result.is_err() {
                    out_one.fill(f64::NAN);
                    out_one[3] = 0.0;
                    out_one[4] = 0.0;
                }
            })
    };

    if n_threads <= 1 {
        for (i_zero, out_one) in out.chunks_mut(COMMONFACTOR_BATCH_OUT_COLS).enumerate() {
            let result = fit_commonfactor_snp_fast(
                k,
                s_ld,
                s_ld_nrow,
                s_ld_ncol,
                v_ld,
                v_ld_nrow,
                v_ld_ncol,
                i_ld,
                i_ld_nrow,
                i_ld_ncol,
                beta_snp,
                beta_nrow,
                beta_ncol,
                se_snp,
                se_nrow,
                se_ncol,
                var_snp,
                coords,
                coords_nrow,
                coords_ncol,
                var_snp_se2,
                gc,
                start,
                max_iter_main,
                max_iter_q,
                tol,
                i_zero,
                out_one,
            );
            if result.is_err() {
                out_one.fill(f64::NAN);
                out_one[3] = 0.0;
                out_one[4] = 0.0;
            }
        }
        return Ok(());
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|_| KernelError::ThreadPoolBuild)?
        .install(run);

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn fit_generic_sem(
    obs_n: usize,
    total_n: usize,
    s_full: &[f64],
    s_nrow: usize,
    s_ncol: usize,
    v_full_reorder: &[f64],
    v_nrow: usize,
    v_ncol: usize,
    w_diag: &[f64],
    b_fixed: &[f64],
    psi_fixed: &[f64],
    b_free: &[i32],
    psi_free: &[i32],
    start: &[f64],
    lower: &[f64],
    upper: &[f64],
    max_iter: usize,
    tol: f64,
    out: &mut [f64],
) -> KernelResult<()> {
    let m = obs_n * (obs_n + 1) / 2;
    let q = start.len();
    let out_len = 2 * q + obs_n * obs_n + 3;
    let total_sq = total_n * total_n;

    if obs_n == 0
        || total_n < obs_n
        || q == 0
        || s_nrow != obs_n
        || s_ncol != obs_n
        || s_full.len() != obs_n * obs_n
        || v_nrow != m
        || v_ncol != m
        || v_full_reorder.len() != m * m
        || w_diag.len() != m
        || b_fixed.len() != total_sq
        || psi_fixed.len() != total_sq
        || b_free.len() != total_sq
        || psi_free.len() != total_sq
        || lower.len() != q
        || upper.len() != q
        || out.len() != out_len
    {
        return Err(KernelError::BadDimensions);
    }
    validate_parameter_bounds(q, lower, upper)?;

    let mut y = vec![0.0; m];
    vech_from_matrix_col_major(s_full, obs_n, &mut y);

    let mut p = start.to_vec();
    project_parameters_to_bounds(&mut p, lower, upper)?;
    let mut sigma = vec![0.0; m];
    let mut delta = vec![0.0; m * q];
    let mut converged = false;
    let mut iterations = 0usize;

    for iter in 0..max_iter {
        iterations = iter + 1;
        generic_sem_implied(
            obs_n, total_n, b_fixed, psi_fixed, b_free, psi_free, &p, &mut sigma, None,
        )?;
        generic_sem_delta_numeric(
            obs_n, total_n, b_fixed, psi_fixed, b_free, psi_free, &p, &sigma, &mut delta,
        )?;

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
            project_parameters_to_bounds(&mut candidate, lower, upper)?;

            let mut sigma_candidate = vec![0.0; m];
            generic_sem_implied(
                obs_n,
                total_n,
                b_fixed,
                psi_fixed,
                b_free,
                psi_free,
                &candidate,
                &mut sigma_candidate,
                None,
            )?;

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
            project_parameters_to_bounds(&mut candidate, lower, upper)?;
        }

        for j in 0..q {
            let scaled = (candidate[j] - p[j]).abs() / (1.0 + p[j].abs());
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

    let mut implied_obs = vec![0.0; obs_n * obs_n];
    generic_sem_implied(
        obs_n,
        total_n,
        b_fixed,
        psi_fixed,
        b_free,
        psi_free,
        &p,
        &mut sigma,
        Some(&mut implied_obs),
    )?;
    generic_sem_delta_numeric(
        obs_n, total_n, b_fixed, psi_fixed, b_free, psi_free, &p, &sigma, &mut delta,
    )?;

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
    out[(2 * q)..(2 * q + obs_n * obs_n)].copy_from_slice(&implied_obs);
    out[2 * q + obs_n * obs_n] = obj;
    out[2 * q + obs_n * obs_n + 1] = if converged { 1.0 } else { 0.0 };
    out[2 * q + obs_n * obs_n + 2] = iterations as f64;

    Ok(())
}

fn reorder_square_col_major(
    input: &[f64],
    input_n: usize,
    order_zero: &[usize],
    out: &mut [f64],
) -> KernelResult<()> {
    let n = order_zero.len();
    if input.len() != input_n * input_n || out.len() != n * n {
        return Err(KernelError::BadDimensions);
    }

    for col in 0..n {
        let source_col = order_zero[col];
        if source_col >= input_n {
            return Err(KernelError::BadIndex);
        }
        for row in 0..n {
            let source_row = order_zero[row];
            if source_row >= input_n {
                return Err(KernelError::BadIndex);
            }
            out[idx(row, col, n)] = input[idx(source_row, source_col, input_n)];
        }
    }

    Ok(())
}

fn vech_residual_original_order(
    s_original: &[f64],
    implied_spec: &[f64],
    obs_n: usize,
    spec_to_original: &[usize],
    out: &mut [f64],
) -> KernelResult<()> {
    let m = obs_n * (obs_n + 1) / 2;
    if s_original.len() != obs_n * obs_n
        || implied_spec.len() != obs_n * obs_n
        || spec_to_original.len() != obs_n
        || out.len() != m
    {
        return Err(KernelError::BadDimensions);
    }

    let mut original_to_spec = vec![usize::MAX; obs_n];
    for (spec_i, &original_i) in spec_to_original.iter().enumerate() {
        if original_i >= obs_n || original_to_spec[original_i] != usize::MAX {
            return Err(KernelError::BadIndex);
        }
        original_to_spec[original_i] = spec_i;
    }

    let mut out_i = 0;
    for col in 0..obs_n {
        let spec_col = original_to_spec[col];
        for row in col..obs_n {
            let spec_row = original_to_spec[row];
            out[out_i] =
                s_original[idx(row, col, obs_n)] - implied_spec[idx(spec_row, spec_col, obs_n)];
            out_i += 1;
        }
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn block_diagonal_gwas_q(
    eta: &[f64],
    k: usize,
    var_snp_se2: f64,
    v_snp: &[f64],
    v_ld_inv_row_major: &[f64],
) -> KernelResult<f64> {
    let obs_n = k + 1;
    let m = obs_n * (obs_n + 1) / 2;
    let ld_m = k * (k + 1) / 2;
    if eta.len() != m
        || v_snp.len() != k * k
        || v_ld_inv_row_major.len() != ld_m * ld_m
        || !var_snp_se2.is_finite()
        || var_snp_se2 == 0.0
    {
        return Err(KernelError::BadDimensions);
    }

    let mut out = eta[0] * eta[0] / var_snp_se2;

    let eta_snp = &eta[1..(k + 1)];
    let mut v_snp_row_major = vec![0.0; k * k];
    for col in 0..k {
        for row in 0..k {
            v_snp_row_major[ridx(row, col, k)] = v_snp[idx(row, col, k)];
        }
    }
    out += quadratic_form_inverse(&v_snp_row_major, eta_snp, k)?;

    let eta_ld = &eta[(k + 1)..];
    for row in 0..ld_m {
        let mut row_sum = 0.0;
        for col in 0..ld_m {
            row_sum += v_ld_inv_row_major[ridx(row, col, ld_m)] * eta_ld[col];
        }
        out += eta_ld[row] * row_sum;
    }

    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn fit_generic_sem_snp_fast(
    obs_n: usize,
    total_n: usize,
    s_ld: &[f64],
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    v_ld: &[f64],
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    v_ld_inv_row_major: &[f64],
    i_ld: &[f64],
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    beta_snp: &[f64],
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: &[f64],
    se_nrow: usize,
    se_ncol: usize,
    var_snp: &[f64],
    coords: &[i32],
    coords_nrow: usize,
    coords_ncol: usize,
    var_snp_se2: f64,
    gc: GenomicControl,
    order_zero: &[usize],
    spec_to_original: &[usize],
    b_fixed: &[f64],
    psi_fixed: &[f64],
    b_free: &[i32],
    psi_free: &[i32],
    start: &[f64],
    lower: &[f64],
    upper: &[f64],
    q_snp_indices_zero: &[i32],
    q_snp_nrow: usize,
    q_snp_ncol: usize,
    q_snp_lengths: &[i32],
    max_iter: usize,
    tol: f64,
    i_zero: usize,
    out: &mut [f64],
) -> KernelResult<()> {
    let k = obs_n - 1;
    let m = obs_n * (obs_n + 1) / 2;
    let q = start.len();
    let out_cols = 2 * q + 1 + q_snp_ncol + 2;
    if out.len() != out_cols || obs_n == 0 || spec_to_original.len() != obs_n {
        return Err(KernelError::BadDimensions);
    }

    let mut s_original = vec![0.0; obs_n * obs_n];
    fill_s_full(
        k,
        s_ld,
        s_ld_nrow,
        s_ld_ncol,
        var_snp,
        beta_snp,
        beta_nrow,
        beta_ncol,
        i_zero,
        &mut s_original,
    )?;

    let mut s_spec = vec![0.0; obs_n * obs_n];
    reorder_square_col_major(&s_original, obs_n, spec_to_original, &mut s_spec)?;

    let mut v_snp = vec![0.0; k * k];
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
        &mut v_snp,
    )?;

    let mut v_full = vec![0.0; m * m];
    fill_v_full(
        k,
        v_ld,
        v_ld_nrow,
        v_ld_ncol,
        var_snp_se2,
        &v_snp,
        k,
        k,
        &mut v_full,
    )?;

    let mut v_reorder = vec![0.0; m * m];
    reorder_square_col_major(&v_full, m, order_zero, &mut v_reorder)?;

    let mut w_diag = vec![0.0; m];
    for d in 0..m {
        let value = v_reorder[idx(d, d, m)];
        if !value.is_finite() || value == 0.0 {
            return Err(KernelError::Singular);
        }
        w_diag[d] = 1.0 / value;
    }

    let mut fit_out = vec![0.0; 2 * q + obs_n * obs_n + 3];
    fit_generic_sem(
        obs_n,
        total_n,
        &s_spec,
        obs_n,
        obs_n,
        &v_reorder,
        m,
        m,
        &w_diag,
        b_fixed,
        psi_fixed,
        b_free,
        psi_free,
        start,
        lower,
        upper,
        max_iter,
        tol,
        &mut fit_out,
    )?;

    let converged = fit_out[2 * q + obs_n * obs_n + 1] == 1.0;
    if !converged
        || !fit_out[..(2 * q + obs_n * obs_n)]
            .iter()
            .all(|v| v.is_finite())
    {
        return Err(KernelError::Singular);
    }

    let implied = &fit_out[(2 * q)..(2 * q + obs_n * obs_n)];
    let mut eta = vec![0.0; m];
    vech_residual_original_order(&s_original, implied, obs_n, spec_to_original, &mut eta)?;
    let chisq = block_diagonal_gwas_q(&eta, k, var_snp_se2, &v_snp, v_ld_inv_row_major)?;
    if !chisq.is_finite() {
        return Err(KernelError::Singular);
    }

    out[..q].copy_from_slice(&fit_out[..q]);
    out[q..(2 * q)].copy_from_slice(&fit_out[q..(2 * q)]);
    out[2 * q] = chisq;

    for q_idx in 0..q_snp_ncol {
        let len = if q_idx < q_snp_lengths.len() {
            q_snp_lengths[q_idx]
        } else {
            0
        };
        let out_idx = 2 * q + 1 + q_idx;
        if len <= 0 {
            out[out_idx] = f64::NAN;
            continue;
        }

        let len = len as usize;
        if len > q_snp_nrow {
            return Err(KernelError::BadDimensions);
        }

        let mut eta_snp = vec![0.0; len];
        let mut v_sub = vec![0.0; len * len];
        for col in 0..len {
            let trait_col_i32 = q_snp_indices_zero[idx(col, q_idx, q_snp_nrow)];
            if trait_col_i32 < 0 {
                out[out_idx] = f64::NAN;
                continue;
            }
            let trait_col = trait_col_i32 as usize;
            if trait_col >= k {
                return Err(KernelError::BadIndex);
            }

            let original_col = trait_col + 1;
            let eta_idx = vech_index(original_col, 0, obs_n);
            eta_snp[col] = eta[eta_idx];

            for row in 0..len {
                let trait_row_i32 = q_snp_indices_zero[idx(row, q_idx, q_snp_nrow)];
                if trait_row_i32 < 0 {
                    out[out_idx] = f64::NAN;
                    continue;
                }
                let trait_row = trait_row_i32 as usize;
                if trait_row >= k {
                    return Err(KernelError::BadIndex);
                }
                v_sub[ridx(row, col, len)] = v_snp[idx(trait_row, trait_col, k)];
            }
        }

        out[out_idx] = quadratic_form_inverse(&v_sub, &eta_snp, len).unwrap_or(f64::NAN);
    }

    out[2 * q + 1 + q_snp_ncol] = 1.0;
    out[2 * q + 1 + q_snp_ncol + 1] = fit_out[2 * q + obs_n * obs_n + 2];

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn fit_generic_sem_batch(
    obs_n: usize,
    total_n: usize,
    s_ld: &[f64],
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    v_ld: &[f64],
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    i_ld: &[f64],
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    beta_snp: &[f64],
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: &[f64],
    se_nrow: usize,
    se_ncol: usize,
    var_snp: &[f64],
    coords: &[i32],
    coords_nrow: usize,
    coords_ncol: usize,
    var_snp_se2: f64,
    gc: GenomicControl,
    order_one: &[i32],
    spec_to_original_one: &[i32],
    b_fixed: &[f64],
    psi_fixed: &[f64],
    b_free: &[i32],
    psi_free: &[i32],
    start: &[f64],
    lower: &[f64],
    upper: &[f64],
    q_snp_indices_one: &[i32],
    q_snp_nrow: usize,
    q_snp_ncol: usize,
    q_snp_lengths: &[i32],
    max_iter: usize,
    tol: f64,
    n_threads: usize,
    out: &mut [f64],
) -> KernelResult<()> {
    let k = obs_n.checked_sub(1).ok_or(KernelError::BadDimensions)?;
    let m = obs_n * (obs_n + 1) / 2;
    let q = start.len();
    let out_cols = 2 * q + 1 + q_snp_ncol + 2;
    let total_sq = total_n * total_n;

    if obs_n == 0
        || total_n < obs_n
        || beta_nrow != se_nrow
        || beta_ncol < k
        || se_ncol < k
        || var_snp.len() != beta_nrow
        || s_ld.len() != s_ld_nrow * s_ld_ncol
        || v_ld.len() != v_ld_nrow * v_ld_ncol
        || i_ld.len() != i_ld_nrow * i_ld_ncol
        || beta_snp.len() != beta_nrow * beta_ncol
        || se_snp.len() != se_nrow * se_ncol
        || coords.len() != coords_nrow * coords_ncol
        || order_one.len() != m
        || spec_to_original_one.len() != obs_n
        || b_fixed.len() != total_sq
        || psi_fixed.len() != total_sq
        || b_free.len() != total_sq
        || psi_free.len() != total_sq
        || lower.len() != q
        || upper.len() != q
        || q_snp_indices_one.len() != q_snp_nrow * q_snp_ncol
        || q_snp_lengths.len() != q_snp_ncol
        || out.len() != beta_nrow * out_cols
    {
        return Err(KernelError::BadDimensions);
    }
    validate_parameter_bounds(q, lower, upper)?;

    let mut order_zero = vec![0usize; m];
    for (i, value) in order_one.iter().enumerate() {
        if *value <= 0 {
            return Err(KernelError::BadIndex);
        }
        order_zero[i] = (*value - 1) as usize;
    }

    let mut spec_to_original = vec![0usize; obs_n];
    for (i, value) in spec_to_original_one.iter().enumerate() {
        if *value <= 0 {
            return Err(KernelError::BadIndex);
        }
        spec_to_original[i] = (*value - 1) as usize;
    }

    let mut q_snp_indices_zero = vec![0i32; q_snp_nrow * q_snp_ncol];
    for (i, value) in q_snp_indices_one.iter().enumerate() {
        q_snp_indices_zero[i] = if *value <= 0 { -1 } else { *value - 1 };
    }

    let ld_m = k * (k + 1) / 2;
    if v_ld_nrow != ld_m || v_ld_ncol != ld_m {
        return Err(KernelError::BadDimensions);
    }
    let mut v_ld_row_major = vec![0.0; ld_m * ld_m];
    for col in 0..ld_m {
        for row in 0..ld_m {
            v_ld_row_major[ridx(row, col, ld_m)] = v_ld[idx(row, col, ld_m)];
        }
    }
    let v_ld_inv_row_major = invert_matrix(&v_ld_row_major, ld_m)?;

    let run = || {
        use rayon::prelude::*;

        out.par_chunks_mut(out_cols)
            .enumerate()
            .for_each(|(i_zero, out_one)| {
                let result = fit_generic_sem_snp_fast(
                    obs_n,
                    total_n,
                    s_ld,
                    s_ld_nrow,
                    s_ld_ncol,
                    v_ld,
                    v_ld_nrow,
                    v_ld_ncol,
                    &v_ld_inv_row_major,
                    i_ld,
                    i_ld_nrow,
                    i_ld_ncol,
                    beta_snp,
                    beta_nrow,
                    beta_ncol,
                    se_snp,
                    se_nrow,
                    se_ncol,
                    var_snp,
                    coords,
                    coords_nrow,
                    coords_ncol,
                    var_snp_se2,
                    gc,
                    &order_zero,
                    &spec_to_original,
                    b_fixed,
                    psi_fixed,
                    b_free,
                    psi_free,
                    start,
                    lower,
                    upper,
                    &q_snp_indices_zero,
                    q_snp_nrow,
                    q_snp_ncol,
                    q_snp_lengths,
                    max_iter,
                    tol,
                    i_zero,
                    out_one,
                );
                if result.is_err() {
                    out_one.fill(f64::NAN);
                    out_one[2 * q + 1 + q_snp_ncol] = 0.0;
                }
            })
    };

    if n_threads <= 1 {
        for (i_zero, out_one) in out.chunks_mut(out_cols).enumerate() {
            let result = fit_generic_sem_snp_fast(
                obs_n,
                total_n,
                s_ld,
                s_ld_nrow,
                s_ld_ncol,
                v_ld,
                v_ld_nrow,
                v_ld_ncol,
                &v_ld_inv_row_major,
                i_ld,
                i_ld_nrow,
                i_ld_ncol,
                beta_snp,
                beta_nrow,
                beta_ncol,
                se_snp,
                se_nrow,
                se_ncol,
                var_snp,
                coords,
                coords_nrow,
                coords_ncol,
                var_snp_se2,
                gc,
                &order_zero,
                &spec_to_original,
                b_fixed,
                psi_fixed,
                b_free,
                psi_free,
                start,
                lower,
                upper,
                &q_snp_indices_zero,
                q_snp_nrow,
                q_snp_ncol,
                q_snp_lengths,
                max_iter,
                tol,
                i_zero,
                out_one,
            );
            if result.is_err() {
                out_one.fill(f64::NAN);
                out_one[2 * q + 1 + q_snp_ncol] = 0.0;
            }
        }
        return Ok(());
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|_| KernelError::ThreadPoolBuild)?
        .install(run);

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
