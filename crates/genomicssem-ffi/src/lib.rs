#![allow(unsafe_op_in_unsafe_fn)]

use genomicssem_core::{
    fill_s_full, fill_v_full, fill_v_snp, fill_v_snp_batch, fill_z_pre, fit_commonfactor_batch,
    fit_commonfactor_main, fit_commonfactor_q, fit_generic_sem, fit_generic_sem_batch,
    ldsc_block_products, munge_fused, munge_fused_batch, munge_qc, sumstats_fused,
    sumstats_fused_batch, sumstats_qc, GenomicControl, KernelError, MungeFusedInput, MungeQcOutput,
    SumstatsFusedInput, SumstatsFusedOutput, SumstatsQcOutput,
};
use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

const OK: i32 = 0;
const BAD_DIMENSIONS: i32 = 1;
const BAD_INDEX: i32 = 2;
const BAD_GC: i32 = 3;
const NULL_POINTER: i32 = 4;
const THREAD_POOL_BUILD: i32 = 5;
const SINGULAR: i32 = 6;

fn code(err: KernelError) -> i32 {
    match err {
        KernelError::BadDimensions => BAD_DIMENSIONS,
        KernelError::BadIndex => BAD_INDEX,
        KernelError::BadGenomicControl => BAD_GC,
        KernelError::NullPointer => NULL_POINTER,
        KernelError::ThreadPoolBuild => THREAD_POOL_BUILD,
        KernelError::Singular => SINGULAR,
    }
}

unsafe fn checked_slice<'a, T>(ptr: *const T, len: usize) -> Result<&'a [T], KernelError> {
    if ptr.is_null() && len != 0 {
        return Err(KernelError::NullPointer);
    }
    Ok(slice::from_raw_parts(ptr, len))
}

unsafe fn checked_slice_mut<'a, T>(ptr: *mut T, len: usize) -> Result<&'a mut [T], KernelError> {
    if ptr.is_null() && len != 0 {
        return Err(KernelError::NullPointer);
    }
    Ok(slice::from_raw_parts_mut(ptr, len))
}

unsafe fn checked_cstr<'a>(ptr: *const c_char) -> Result<&'a str, KernelError> {
    if ptr.is_null() {
        return Err(KernelError::NullPointer);
    }
    CStr::from_ptr(ptr)
        .to_str()
        .map_err(|_| KernelError::BadDimensions)
}

unsafe fn checked_cstrs<'a>(
    ptr: *const *const c_char,
    len: usize,
) -> Result<Vec<&'a str>, KernelError> {
    if ptr.is_null() && len != 0 {
        return Err(KernelError::NullPointer);
    }

    let raw = slice::from_raw_parts(ptr, len);
    let mut out = Vec::with_capacity(len);
    for item in raw {
        out.push(checked_cstr(*item)?);
    }
    Ok(out)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fill_v_snp(
    se_snp: *const f64,
    se_nrow: usize,
    se_ncol: usize,
    i_zero: usize,
    i_ld: *const f64,
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    var_snp: *const f64,
    var_snp_len: usize,
    coords: *const i32,
    coords_nrow: usize,
    coords_ncol: usize,
    k: usize,
    gc_code: i32,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let Some(gc) = GenomicControl::from_code(gc_code) else {
        return BAD_GC;
    };

    let result = (|| {
        let se_snp = checked_slice(se_snp, se_nrow * se_ncol)?;
        let i_ld = checked_slice(i_ld, i_ld_nrow * i_ld_ncol)?;
        let var_snp = checked_slice(var_snp, var_snp_len)?;
        let coords = checked_slice(coords, coords_nrow * coords_ncol)?;
        let out = checked_slice_mut(out, out_len)?;

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
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fill_v_full(
    k: usize,
    v_ld: *const f64,
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    var_snp_se2: f64,
    v_snp: *const f64,
    v_snp_nrow: usize,
    v_snp_ncol: usize,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let result = (|| {
        let v_ld = checked_slice(v_ld, v_ld_nrow * v_ld_ncol)?;
        let v_snp = checked_slice(v_snp, v_snp_nrow * v_snp_ncol)?;
        let out = checked_slice_mut(out, out_len)?;

        fill_v_full(
            k,
            v_ld,
            v_ld_nrow,
            v_ld_ncol,
            var_snp_se2,
            v_snp,
            v_snp_nrow,
            v_snp_ncol,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fill_v_snp_batch(
    se_snp: *const f64,
    se_nrow: usize,
    se_ncol: usize,
    i_ld: *const f64,
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    var_snp: *const f64,
    var_snp_len: usize,
    coords: *const i32,
    coords_nrow: usize,
    coords_ncol: usize,
    k: usize,
    gc_code: i32,
    n_threads: usize,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let Some(gc) = GenomicControl::from_code(gc_code) else {
        return BAD_GC;
    };

    let result = (|| {
        let se_snp = checked_slice(se_snp, se_nrow * se_ncol)?;
        let i_ld = checked_slice(i_ld, i_ld_nrow * i_ld_ncol)?;
        let var_snp = checked_slice(var_snp, var_snp_len)?;
        let coords = checked_slice(coords, coords_nrow * coords_ncol)?;
        let out = checked_slice_mut(out, out_len)?;

        fill_v_snp_batch(
            se_snp,
            se_nrow,
            se_ncol,
            i_ld,
            i_ld_nrow,
            i_ld_ncol,
            var_snp,
            coords,
            coords_nrow,
            coords_ncol,
            k,
            gc,
            n_threads,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fill_s_full(
    k: usize,
    s_ld: *const f64,
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    var_snp: *const f64,
    var_snp_len: usize,
    beta_snp: *const f64,
    beta_nrow: usize,
    beta_ncol: usize,
    i_zero: usize,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let result = (|| {
        let s_ld = checked_slice(s_ld, s_ld_nrow * s_ld_ncol)?;
        let var_snp = checked_slice(var_snp, var_snp_len)?;
        let beta_snp = checked_slice(beta_snp, beta_nrow * beta_ncol)?;
        let out = checked_slice_mut(out, out_len)?;

        fill_s_full(
            k, s_ld, s_ld_nrow, s_ld_ncol, var_snp, beta_snp, beta_nrow, beta_ncol, i_zero, out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fill_z_pre(
    beta_snp: *const f64,
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: *const f64,
    se_nrow: usize,
    se_ncol: usize,
    i_ld: *const f64,
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    i_zero: usize,
    gc_code: i32,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let Some(gc) = GenomicControl::from_code(gc_code) else {
        return BAD_GC;
    };

    let result = (|| {
        let beta_snp = checked_slice(beta_snp, beta_nrow * beta_ncol)?;
        let se_snp = checked_slice(se_snp, se_nrow * se_ncol)?;
        let i_ld = checked_slice(i_ld, i_ld_nrow * i_ld_ncol)?;
        let out = checked_slice_mut(out, out_len)?;

        fill_z_pre(
            beta_snp, beta_nrow, beta_ncol, se_snp, se_nrow, se_ncol, i_ld, i_ld_nrow, i_ld_ncol,
            i_zero, gc, out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fit_commonfactor_main(
    k: usize,
    s_full: *const f64,
    s_nrow: usize,
    s_ncol: usize,
    v_full_reorder: *const f64,
    v_nrow: usize,
    v_ncol: usize,
    w_diag: *const f64,
    w_len: usize,
    start: *const f64,
    start_len: usize,
    max_iter: usize,
    tol: f64,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let result = (|| {
        let s_full = checked_slice(s_full, s_nrow * s_ncol)?;
        let v_full_reorder = checked_slice(v_full_reorder, v_nrow * v_ncol)?;
        let w_diag = checked_slice(w_diag, w_len)?;
        let start = checked_slice(start, start_len)?;
        let out = checked_slice_mut(out, out_len)?;

        fit_commonfactor_main(
            k,
            s_full,
            s_nrow,
            s_ncol,
            v_full_reorder,
            v_nrow,
            v_ncol,
            w_diag,
            start,
            max_iter,
            tol,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fit_commonfactor_q(
    k: usize,
    s_full: *const f64,
    s_nrow: usize,
    s_ncol: usize,
    v_full_reorder: *const f64,
    v_nrow: usize,
    v_ncol: usize,
    w_diag: *const f64,
    w_len: usize,
    fixed: *const f64,
    fixed_len: usize,
    start: *const f64,
    start_len: usize,
    max_iter: usize,
    tol: f64,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let result = (|| {
        let s_full = checked_slice(s_full, s_nrow * s_ncol)?;
        let v_full_reorder = checked_slice(v_full_reorder, v_nrow * v_ncol)?;
        let w_diag = checked_slice(w_diag, w_len)?;
        let fixed = checked_slice(fixed, fixed_len)?;
        let start = checked_slice(start, start_len)?;
        let out = checked_slice_mut(out, out_len)?;

        fit_commonfactor_q(
            k,
            s_full,
            s_nrow,
            s_ncol,
            v_full_reorder,
            v_nrow,
            v_ncol,
            w_diag,
            fixed,
            start,
            max_iter,
            tol,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fit_commonfactor_batch(
    k: usize,
    s_ld: *const f64,
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    v_ld: *const f64,
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    i_ld: *const f64,
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    beta_snp: *const f64,
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: *const f64,
    se_nrow: usize,
    se_ncol: usize,
    var_snp: *const f64,
    var_snp_len: usize,
    coords: *const i32,
    coords_nrow: usize,
    coords_ncol: usize,
    var_snp_se2: f64,
    gc_code: i32,
    start: *const f64,
    start_len: usize,
    max_iter_main: usize,
    max_iter_q: usize,
    tol: f64,
    n_threads: usize,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let Some(gc) = GenomicControl::from_code(gc_code) else {
        return BAD_GC;
    };

    let result = (|| {
        let s_ld = checked_slice(s_ld, s_ld_nrow * s_ld_ncol)?;
        let v_ld = checked_slice(v_ld, v_ld_nrow * v_ld_ncol)?;
        let i_ld = checked_slice(i_ld, i_ld_nrow * i_ld_ncol)?;
        let beta_snp = checked_slice(beta_snp, beta_nrow * beta_ncol)?;
        let se_snp = checked_slice(se_snp, se_nrow * se_ncol)?;
        let var_snp = checked_slice(var_snp, var_snp_len)?;
        let coords = checked_slice(coords, coords_nrow * coords_ncol)?;
        let start = checked_slice(start, start_len)?;
        let out = checked_slice_mut(out, out_len)?;

        fit_commonfactor_batch(
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
            n_threads,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fit_generic_sem(
    obs_n: usize,
    total_n: usize,
    s_full: *const f64,
    s_nrow: usize,
    s_ncol: usize,
    v_full_reorder: *const f64,
    v_nrow: usize,
    v_ncol: usize,
    w_diag: *const f64,
    w_len: usize,
    b_fixed: *const f64,
    b_fixed_len: usize,
    psi_fixed: *const f64,
    psi_fixed_len: usize,
    b_free: *const i32,
    b_free_len: usize,
    psi_free: *const i32,
    psi_free_len: usize,
    start: *const f64,
    start_len: usize,
    lower: *const f64,
    lower_len: usize,
    upper: *const f64,
    upper_len: usize,
    max_iter: usize,
    tol: f64,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let result = (|| {
        let s_full = checked_slice(s_full, s_nrow * s_ncol)?;
        let v_full_reorder = checked_slice(v_full_reorder, v_nrow * v_ncol)?;
        let w_diag = checked_slice(w_diag, w_len)?;
        let b_fixed = checked_slice(b_fixed, b_fixed_len)?;
        let psi_fixed = checked_slice(psi_fixed, psi_fixed_len)?;
        let b_free = checked_slice(b_free, b_free_len)?;
        let psi_free = checked_slice(psi_free, psi_free_len)?;
        let start = checked_slice(start, start_len)?;
        let lower = checked_slice(lower, lower_len)?;
        let upper = checked_slice(upper, upper_len)?;
        let out = checked_slice_mut(out, out_len)?;

        fit_generic_sem(
            obs_n,
            total_n,
            s_full,
            s_nrow,
            s_ncol,
            v_full_reorder,
            v_nrow,
            v_ncol,
            w_diag,
            b_fixed,
            psi_fixed,
            b_free,
            psi_free,
            start,
            lower,
            upper,
            max_iter,
            tol,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_fit_generic_sem_batch(
    obs_n: usize,
    total_n: usize,
    s_ld: *const f64,
    s_ld_nrow: usize,
    s_ld_ncol: usize,
    v_ld: *const f64,
    v_ld_nrow: usize,
    v_ld_ncol: usize,
    i_ld: *const f64,
    i_ld_nrow: usize,
    i_ld_ncol: usize,
    beta_snp: *const f64,
    beta_nrow: usize,
    beta_ncol: usize,
    se_snp: *const f64,
    se_nrow: usize,
    se_ncol: usize,
    var_snp: *const f64,
    var_snp_len: usize,
    coords: *const i32,
    coords_nrow: usize,
    coords_ncol: usize,
    var_snp_se2: f64,
    gc_code: i32,
    order: *const i32,
    order_len: usize,
    spec_to_original: *const i32,
    spec_to_original_len: usize,
    b_fixed: *const f64,
    b_fixed_len: usize,
    psi_fixed: *const f64,
    psi_fixed_len: usize,
    b_free: *const i32,
    b_free_len: usize,
    psi_free: *const i32,
    psi_free_len: usize,
    start: *const f64,
    start_len: usize,
    lower: *const f64,
    lower_len: usize,
    upper: *const f64,
    upper_len: usize,
    q_snp_indices: *const i32,
    q_snp_nrow: usize,
    q_snp_ncol: usize,
    q_snp_lengths: *const i32,
    q_snp_lengths_len: usize,
    max_iter: usize,
    tol: f64,
    n_threads: usize,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    let Some(gc) = GenomicControl::from_code(gc_code) else {
        return BAD_GC;
    };

    let result = (|| {
        let s_ld = checked_slice(s_ld, s_ld_nrow * s_ld_ncol)?;
        let v_ld = checked_slice(v_ld, v_ld_nrow * v_ld_ncol)?;
        let i_ld = checked_slice(i_ld, i_ld_nrow * i_ld_ncol)?;
        let beta_snp = checked_slice(beta_snp, beta_nrow * beta_ncol)?;
        let se_snp = checked_slice(se_snp, se_nrow * se_ncol)?;
        let var_snp = checked_slice(var_snp, var_snp_len)?;
        let coords = checked_slice(coords, coords_nrow * coords_ncol)?;
        let order = checked_slice(order, order_len)?;
        let spec_to_original = checked_slice(spec_to_original, spec_to_original_len)?;
        let b_fixed = checked_slice(b_fixed, b_fixed_len)?;
        let psi_fixed = checked_slice(psi_fixed, psi_fixed_len)?;
        let b_free = checked_slice(b_free, b_free_len)?;
        let psi_free = checked_slice(psi_free, psi_free_len)?;
        let start = checked_slice(start, start_len)?;
        let lower = checked_slice(lower, lower_len)?;
        let upper = checked_slice(upper, upper_len)?;
        let q_snp_indices = checked_slice(q_snp_indices, q_snp_nrow * q_snp_ncol)?;
        let q_snp_lengths = checked_slice(q_snp_lengths, q_snp_lengths_len)?;
        let out = checked_slice_mut(out, out_len)?;

        fit_generic_sem_batch(
            obs_n,
            total_n,
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
            order,
            spec_to_original,
            b_fixed,
            psi_fixed,
            b_free,
            psi_free,
            start,
            lower,
            upper,
            q_snp_indices,
            q_snp_nrow,
            q_snp_ncol,
            q_snp_lengths,
            max_iter,
            tol,
            n_threads,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_munge_qc(
    a1_ref: *const i32,
    a2_ref: *const i32,
    a1_file: *const i32,
    a2_file: *const i32,
    effect: *const f64,
    p: *const f64,
    n: usize,
    info: *const f64,
    info_len: usize,
    maf: *const f64,
    maf_len: usize,
    info_filter: f64,
    maf_filter: f64,
    keep: *mut i32,
    z: *mut f64,
    out_len: usize,
    counts: *mut i32,
    counts_len: usize,
    out_count: *mut usize,
) -> i32 {
    let result = (|| {
        let a1_ref = checked_slice(a1_ref, n)?;
        let a2_ref = checked_slice(a2_ref, n)?;
        let a1_file = checked_slice(a1_file, n)?;
        let a2_file = checked_slice(a2_file, n)?;
        let effect = checked_slice(effect, n)?;
        let p = checked_slice(p, n)?;
        let info = if info_len == 0 {
            None
        } else {
            Some(checked_slice(info, info_len)?)
        };
        let maf = if maf_len == 0 {
            None
        } else {
            Some(checked_slice(maf, maf_len)?)
        };
        let keep = checked_slice_mut(keep, out_len)?;
        let z = checked_slice_mut(z, out_len)?;
        let counts = checked_slice_mut(counts, counts_len)?;
        if out_count.is_null() {
            return Err(KernelError::NullPointer);
        }

        let mut out = MungeQcOutput { keep, z, counts };
        *out_count = munge_qc(
            a1_ref,
            a2_ref,
            a1_file,
            a2_file,
            effect,
            p,
            info,
            maf,
            info_filter,
            maf_filter,
            &mut out,
        )?;
        Ok(())
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_munge_fused(
    filename: *const c_char,
    output_path: *const c_char,
    ref_snp: *const *const c_char,
    ref_len: usize,
    ref_a1: *const i32,
    ref_a1_len: usize,
    ref_a2: *const i32,
    ref_a2_len: usize,
    col_indices: *const i32,
    col_indices_len: usize,
    provided_n: f64,
    n_multiplier: f64,
    info_filter: f64,
    maf_filter: f64,
    out_counts: *mut i32,
    out_counts_len: usize,
    rows_total: *mut usize,
    rows_joined: *mut usize,
    rows_written: *mut usize,
    unsupported: *mut i32,
) -> i32 {
    let result = (|| {
        if rows_total.is_null()
            || rows_joined.is_null()
            || rows_written.is_null()
            || unsupported.is_null()
        {
            return Err(KernelError::NullPointer);
        }
        let filename = checked_cstr(filename)?;
        let output_path = checked_cstr(output_path)?;
        let ref_snp = checked_cstrs(ref_snp, ref_len)?;
        let ref_a1 = checked_slice(ref_a1, ref_a1_len)?;
        let ref_a2 = checked_slice(ref_a2, ref_a2_len)?;
        let col_indices = checked_slice(col_indices, col_indices_len)?;
        let out_counts = checked_slice_mut(out_counts, out_counts_len)?;
        if out_counts_len < genomicssem_core::MUNGE_QC_COUNT_LEN {
            return Err(KernelError::BadDimensions);
        }

        let report = munge_fused(
            filename,
            output_path,
            &ref_snp,
            ref_a1,
            ref_a2,
            col_indices,
            provided_n.is_finite().then_some(provided_n),
            n_multiplier,
            info_filter,
            maf_filter,
        )?;

        out_counts[..genomicssem_core::MUNGE_QC_COUNT_LEN].copy_from_slice(&report.counts);
        *rows_total = report.rows_total;
        *rows_joined = report.rows_joined;
        *rows_written = report.rows_written;
        *unsupported = i32::from(report.unsupported);
        Ok(())
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_munge_fused_batch(
    filenames: *const *const c_char,
    output_paths: *const *const c_char,
    file_count: usize,
    ref_snp: *const *const c_char,
    ref_len: usize,
    ref_a1: *const i32,
    ref_a1_len: usize,
    ref_a2: *const i32,
    ref_a2_len: usize,
    col_indices: *const i32,
    col_indices_nrow: usize,
    col_indices_ncol: usize,
    provided_n: *const f64,
    provided_n_len: usize,
    n_multipliers: *const f64,
    n_multipliers_len: usize,
    info_filter: f64,
    maf_filter: f64,
    n_threads: usize,
    out_counts: *mut i32,
    out_counts_len: usize,
    out_rows: *mut usize,
    out_rows_len: usize,
    unsupported: *mut i32,
    unsupported_len: usize,
) -> i32 {
    let result = (|| {
        let filenames = checked_cstrs(filenames, file_count)?;
        let output_paths = checked_cstrs(output_paths, file_count)?;
        let ref_snp = checked_cstrs(ref_snp, ref_len)?;
        let ref_a1 = checked_slice(ref_a1, ref_a1_len)?;
        let ref_a2 = checked_slice(ref_a2, ref_a2_len)?;
        let col_indices = checked_slice(col_indices, col_indices_nrow * col_indices_ncol)?;
        let provided_n = checked_slice(provided_n, provided_n_len)?;
        let n_multipliers = checked_slice(n_multipliers, n_multipliers_len)?;
        let out_counts = checked_slice_mut(out_counts, out_counts_len)?;
        let out_rows = checked_slice_mut(out_rows, out_rows_len)?;
        let unsupported = checked_slice_mut(unsupported, unsupported_len)?;

        if output_paths.len() != file_count
            || provided_n.len() != file_count
            || n_multipliers.len() != file_count
            || col_indices_nrow != file_count
            || col_indices_ncol < genomicssem_core::MUNGE_FUSED_COLS
            || out_counts.len() < file_count * genomicssem_core::MUNGE_QC_COUNT_LEN
            || out_rows.len() < file_count * 3
            || unsupported.len() < file_count
        {
            return Err(KernelError::BadDimensions);
        }

        let mut per_file_col_indices = Vec::with_capacity(file_count);
        for file_i in 0..file_count {
            let mut file_cols = Vec::with_capacity(col_indices_ncol);
            for col_i in 0..col_indices_ncol {
                file_cols.push(col_indices[file_i + col_i * col_indices_nrow]);
            }
            per_file_col_indices.push(file_cols);
        }

        let inputs: Vec<_> = (0..file_count)
            .map(|file_i| MungeFusedInput {
                filename: filenames[file_i],
                output_path: output_paths[file_i],
                col_indices: per_file_col_indices[file_i].as_slice(),
                provided_n: provided_n[file_i].is_finite().then_some(provided_n[file_i]),
                n_multiplier: n_multipliers[file_i],
            })
            .collect();

        let reports = munge_fused_batch(
            &inputs,
            &ref_snp,
            ref_a1,
            ref_a2,
            info_filter,
            maf_filter,
            n_threads,
        )?;

        for (file_i, report) in reports.iter().enumerate() {
            for count_i in 0..genomicssem_core::MUNGE_QC_COUNT_LEN {
                out_counts[file_i + count_i * file_count] = report.counts[count_i];
            }
            out_rows[file_i] = report.rows_total;
            out_rows[file_i + file_count] = report.rows_joined;
            out_rows[file_i + 2 * file_count] = report.rows_written;
            unsupported[file_i] = i32::from(report.unsupported);
        }
        Ok(())
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_sumstats_qc(
    a1_ref: *const i32,
    a2_ref: *const i32,
    a1_file: *const i32,
    a2_file: *const i32,
    effect: *const f64,
    se: *const f64,
    p: *const f64,
    n_values: *const f64,
    maf_ref: *const f64,
    n: usize,
    maf_file: *const f64,
    maf_file_len: usize,
    info: *const f64,
    info_len: usize,
    info_filter: f64,
    ols: i32,
    beta_is_character: i32,
    linprob: i32,
    se_logit: i32,
    keep: *mut i32,
    beta: *mut f64,
    se_out: *mut f64,
    out_len: usize,
    counts: *mut i32,
    counts_len: usize,
    out_count: *mut usize,
) -> i32 {
    let result = (|| {
        let a1_ref = checked_slice(a1_ref, n)?;
        let a2_ref = checked_slice(a2_ref, n)?;
        let a1_file = checked_slice(a1_file, n)?;
        let a2_file = checked_slice(a2_file, n)?;
        let effect = checked_slice(effect, n)?;
        let se = checked_slice(se, n)?;
        let p = checked_slice(p, n)?;
        let n_values = checked_slice(n_values, n)?;
        let maf_ref = checked_slice(maf_ref, n)?;
        let maf_file = if maf_file_len == 0 {
            None
        } else {
            Some(checked_slice(maf_file, maf_file_len)?)
        };
        let info = if info_len == 0 {
            None
        } else {
            Some(checked_slice(info, info_len)?)
        };
        let keep = checked_slice_mut(keep, out_len)?;
        let beta = checked_slice_mut(beta, out_len)?;
        let se_out = checked_slice_mut(se_out, out_len)?;
        let counts = checked_slice_mut(counts, counts_len)?;
        if out_count.is_null() {
            return Err(KernelError::NullPointer);
        }

        let mut out = SumstatsQcOutput {
            keep,
            beta,
            se_out,
            counts,
        };
        *out_count = sumstats_qc(
            a1_ref,
            a2_ref,
            a1_file,
            a2_file,
            effect,
            se,
            p,
            n_values,
            maf_ref,
            maf_file,
            info,
            info_filter,
            ols != 0,
            beta_is_character != 0,
            linprob != 0,
            se_logit != 0,
            &mut out,
        )?;
        Ok(())
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_sumstats_fused(
    filename: *const c_char,
    ref_snp: *const *const c_char,
    ref_len: usize,
    ref_a1: *const i32,
    ref_a1_len: usize,
    ref_a2: *const i32,
    ref_a2_len: usize,
    ref_maf: *const f64,
    ref_maf_len: usize,
    col_indices: *const i32,
    col_indices_len: usize,
    provided_n: f64,
    info_filter: f64,
    ols: i32,
    beta_is_character: i32,
    linprob: i32,
    se_logit: i32,
    keep: *mut i32,
    beta: *mut f64,
    se: *mut f64,
    out_len: usize,
    out_counts: *mut i32,
    out_counts_len: usize,
    rows_total: *mut usize,
    rows_duplicate_removed: *mut usize,
    rows_joined: *mut usize,
    rows_written: *mut usize,
    unsupported: *mut i32,
) -> i32 {
    let result = (|| {
        if rows_total.is_null()
            || rows_duplicate_removed.is_null()
            || rows_joined.is_null()
            || rows_written.is_null()
            || unsupported.is_null()
        {
            return Err(KernelError::NullPointer);
        }
        let filename = checked_cstr(filename)?;
        let ref_snp = checked_cstrs(ref_snp, ref_len)?;
        let ref_a1 = checked_slice(ref_a1, ref_a1_len)?;
        let ref_a2 = checked_slice(ref_a2, ref_a2_len)?;
        let ref_maf = checked_slice(ref_maf, ref_maf_len)?;
        let col_indices = checked_slice(col_indices, col_indices_len)?;
        let keep = checked_slice_mut(keep, out_len)?;
        let beta = checked_slice_mut(beta, out_len)?;
        let se = checked_slice_mut(se, out_len)?;
        let out_counts = checked_slice_mut(out_counts, out_counts_len)?;
        if out_counts_len < genomicssem_core::SUMSTATS_QC_COUNT_LEN {
            return Err(KernelError::BadDimensions);
        }

        let report = sumstats_fused(
            filename,
            &ref_snp,
            ref_a1,
            ref_a2,
            ref_maf,
            col_indices,
            provided_n.is_finite().then_some(provided_n),
            info_filter,
            ols != 0,
            beta_is_character != 0,
            linprob != 0,
            se_logit != 0,
            &mut SumstatsFusedOutput { keep, beta, se },
        )?;

        out_counts[..genomicssem_core::SUMSTATS_QC_COUNT_LEN].copy_from_slice(&report.counts);
        *rows_total = report.rows_total;
        *rows_duplicate_removed = report.rows_duplicate_removed;
        *rows_joined = report.rows_joined;
        *rows_written = report.rows_written;
        *unsupported = i32::from(report.unsupported);
        Ok(())
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_sumstats_fused_batch(
    filenames: *const *const c_char,
    file_count: usize,
    ref_snp: *const *const c_char,
    ref_len: usize,
    ref_a1: *const i32,
    ref_a1_len: usize,
    ref_a2: *const i32,
    ref_a2_len: usize,
    ref_maf: *const f64,
    ref_maf_len: usize,
    col_indices: *const i32,
    col_indices_nrow: usize,
    col_indices_ncol: usize,
    provided_n: *const f64,
    provided_n_len: usize,
    ols: *const i32,
    ols_len: usize,
    beta_is_character: *const i32,
    beta_is_character_len: usize,
    linprob: *const i32,
    linprob_len: usize,
    se_logit: *const i32,
    se_logit_len: usize,
    info_filter: f64,
    n_threads: usize,
    keep: *mut i32,
    keep_len: usize,
    beta: *mut f64,
    beta_len: usize,
    se: *mut f64,
    se_len: usize,
    out_counts: *mut i32,
    out_counts_len: usize,
    out_rows: *mut usize,
    out_rows_len: usize,
    out_mean_abs_z: *mut f64,
    out_mean_abs_z_len: usize,
    out_n: *mut usize,
    unsupported: *mut i32,
    unsupported_len: usize,
) -> i32 {
    let result = (|| {
        if out_n.is_null() {
            return Err(KernelError::NullPointer);
        }
        let filenames = checked_cstrs(filenames, file_count)?;
        let ref_snp = checked_cstrs(ref_snp, ref_len)?;
        let ref_a1 = checked_slice(ref_a1, ref_a1_len)?;
        let ref_a2 = checked_slice(ref_a2, ref_a2_len)?;
        let ref_maf = checked_slice(ref_maf, ref_maf_len)?;
        let col_indices = checked_slice(col_indices, col_indices_nrow * col_indices_ncol)?;
        let provided_n = checked_slice(provided_n, provided_n_len)?;
        let ols = checked_slice(ols, ols_len)?;
        let beta_is_character = checked_slice(beta_is_character, beta_is_character_len)?;
        let linprob = checked_slice(linprob, linprob_len)?;
        let se_logit = checked_slice(se_logit, se_logit_len)?;
        let keep = checked_slice_mut(keep, keep_len)?;
        let beta = checked_slice_mut(beta, beta_len)?;
        let se = checked_slice_mut(se, se_len)?;
        let out_counts = checked_slice_mut(out_counts, out_counts_len)?;
        let out_rows = checked_slice_mut(out_rows, out_rows_len)?;
        let out_mean_abs_z = checked_slice_mut(out_mean_abs_z, out_mean_abs_z_len)?;
        let unsupported = checked_slice_mut(unsupported, unsupported_len)?;

        if ref_a1.len() != ref_len
            || ref_a2.len() != ref_len
            || ref_maf.len() != ref_len
            || provided_n.len() != file_count
            || ols.len() != file_count
            || beta_is_character.len() != file_count
            || linprob.len() != file_count
            || se_logit.len() != file_count
            || col_indices_nrow != file_count
            || col_indices_ncol < genomicssem_core::SUMSTATS_FUSED_COLS
            || keep.len() < ref_len
            || beta.len() < ref_len * file_count
            || se.len() < ref_len * file_count
            || out_counts.len() < file_count * genomicssem_core::SUMSTATS_QC_COUNT_LEN
            || out_rows.len() < file_count * 4
            || out_mean_abs_z.len() < file_count
            || unsupported.len() < file_count
        {
            return Err(KernelError::BadDimensions);
        }

        let mut per_file_col_indices = Vec::with_capacity(file_count);
        for file_i in 0..file_count {
            let mut file_cols = Vec::with_capacity(col_indices_ncol);
            for col_i in 0..col_indices_ncol {
                file_cols.push(col_indices[file_i + col_i * col_indices_nrow]);
            }
            per_file_col_indices.push(file_cols);
        }

        let inputs: Vec<_> = (0..file_count)
            .map(|file_i| SumstatsFusedInput {
                filename: filenames[file_i],
                col_indices: per_file_col_indices[file_i].as_slice(),
                provided_n: provided_n[file_i].is_finite().then_some(provided_n[file_i]),
                ols: ols[file_i] != 0,
                beta_is_character: beta_is_character[file_i] != 0,
                linprob: linprob[file_i] != 0,
                se_logit: se_logit[file_i] != 0,
            })
            .collect();

        let output = sumstats_fused_batch(
            &inputs,
            &ref_snp,
            ref_a1,
            ref_a2,
            ref_maf,
            info_filter,
            n_threads,
        )?;

        for (file_i, report) in output.reports.iter().enumerate() {
            for count_i in 0..genomicssem_core::SUMSTATS_QC_COUNT_LEN {
                out_counts[file_i + count_i * file_count] = report.counts[count_i];
            }
            out_rows[file_i] = report.rows_total;
            out_rows[file_i + file_count] = report.rows_duplicate_removed;
            out_rows[file_i + 2 * file_count] = report.rows_joined;
            out_rows[file_i + 3 * file_count] = report.rows_written;
            out_mean_abs_z[file_i] = report.mean_abs_z;
            unsupported[file_i] = i32::from(report.unsupported);
        }

        *out_n = output.keep.len();
        keep[..output.keep.len()].copy_from_slice(&output.keep);
        beta[..output.beta.len()].copy_from_slice(&output.beta);
        se[..output.se.len()].copy_from_slice(&output.se);
        Ok(())
    })();

    result.map(|()| OK).unwrap_or_else(code)
}

#[no_mangle]
pub unsafe extern "C" fn genomicssem_ldsc_block_products(
    weighted_ld: *const f64,
    n_snps: usize,
    n_annot: usize,
    weighted_chi: *const f64,
    chi_len: usize,
    n_blocks: usize,
    xty_block: *mut f64,
    xty_block_len: usize,
    xtx_block: *mut f64,
    xtx_block_len: usize,
    xty: *mut f64,
    xty_len: usize,
    xtx: *mut f64,
    xtx_len: usize,
    n_threads: usize,
) -> i32 {
    let result = (|| {
        let weighted_ld = checked_slice(weighted_ld, n_snps * n_annot)?;
        let weighted_chi = checked_slice(weighted_chi, chi_len)?;
        let xty_block = checked_slice_mut(xty_block, xty_block_len)?;
        let xtx_block = checked_slice_mut(xtx_block, xtx_block_len)?;
        let xty = checked_slice_mut(xty, xty_len)?;
        let xtx = checked_slice_mut(xtx, xtx_len)?;

        ldsc_block_products(
            weighted_ld,
            n_snps,
            n_annot,
            weighted_chi,
            n_blocks,
            xty_block,
            xtx_block,
            xty,
            xtx,
            n_threads,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}
