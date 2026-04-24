#![allow(unsafe_op_in_unsafe_fn)]

use genomicssem_core::{
    fill_s_full, fill_v_full, fill_v_snp, fill_v_snp_batch, fill_z_pre, fit_commonfactor_batch,
    fit_commonfactor_main, fit_commonfactor_q, fit_generic_sem, GenomicControl, KernelError,
};
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
            max_iter,
            tol,
            out,
        )
    })();

    result.map(|()| OK).unwrap_or_else(code)
}
