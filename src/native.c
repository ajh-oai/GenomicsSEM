#include <R.h>
#include <Rinternals.h>
#include <R_ext/Error.h>

#include <stddef.h>
#include <string.h>

#include "genomicssem_native.h"

static void matrix_dims(SEXP x, const char *name, size_t *nrow, size_t *ncol) {
  if (!Rf_isMatrix(x)) {
    Rf_error("'%s' must be a matrix", name);
  }

  SEXP dims = Rf_getAttrib(x, R_DimSymbol);
  if (dims == R_NilValue || TYPEOF(dims) != INTSXP || Rf_length(dims) != 2) {
    Rf_error("'%s' must have two integer dimensions", name);
  }

  int nr = INTEGER(dims)[0];
  int nc = INTEGER(dims)[1];
  if (nr < 0 || nc < 0) {
    Rf_error("'%s' dimensions must be non-negative", name);
  }

  *nrow = (size_t)nr;
  *ncol = (size_t)nc;
}

static SEXP protect_real_matrix(SEXP x, const char *name, int *nprotect) {
  size_t nrow, ncol;
  matrix_dims(x, name, &nrow, &ncol);
  if (TYPEOF(x) == REALSXP) {
    return x;
  }
  if (!Rf_isNumeric(x)) {
    Rf_error("'%s' must be numeric", name);
  }
  ++(*nprotect);
  SEXP out = PROTECT(Rf_coerceVector(x, REALSXP));
  if (!Rf_isMatrix(out)) {
    SEXP dims = PROTECT(Rf_allocVector(INTSXP, 2));
    ++(*nprotect);
    INTEGER(dims)[0] = (int)nrow;
    INTEGER(dims)[1] = (int)ncol;
    Rf_setAttrib(out, R_DimSymbol, dims);
  }
  return out;
}

static SEXP protect_real_vector(SEXP x, const char *name, int *nprotect) {
  if (!Rf_isNumeric(x) || Rf_isMatrix(x)) {
    Rf_error("'%s' must be a numeric vector", name);
  }
  if (TYPEOF(x) == REALSXP) {
    return x;
  }
  ++(*nprotect);
  return PROTECT(Rf_coerceVector(x, REALSXP));
}

static SEXP protect_int_matrix(SEXP x, const char *name, int *nprotect) {
  size_t nrow, ncol;
  matrix_dims(x, name, &nrow, &ncol);
  if (TYPEOF(x) == INTSXP) {
    return x;
  }
  if (!Rf_isNumeric(x)) {
    Rf_error("'%s' must be numeric", name);
  }
  ++(*nprotect);
  SEXP out = PROTECT(Rf_coerceVector(x, INTSXP));
  if (!Rf_isMatrix(out)) {
    SEXP dims = PROTECT(Rf_allocVector(INTSXP, 2));
    ++(*nprotect);
    INTEGER(dims)[0] = (int)nrow;
    INTEGER(dims)[1] = (int)ncol;
    Rf_setAttrib(out, R_DimSymbol, dims);
  }
  return out;
}

static int scalar_int(SEXP x, const char *name) {
  int value = Rf_asInteger(x);
  if (value == NA_INTEGER) {
    Rf_error("'%s' must be an integer scalar", name);
  }
  return value;
}

static double scalar_real(SEXP x, const char *name) {
  double value = Rf_asReal(x);
  if (!R_FINITE(value)) {
    Rf_error("'%s' must be a finite numeric scalar", name);
  }
  return value;
}

static int gc_code(SEXP x) {
  if (TYPEOF(x) != STRSXP || XLENGTH(x) != 1 || STRING_ELT(x, 0) == NA_STRING) {
    Rf_error("'GC' must be one of 'conserv', 'standard', or 'none'");
  }

  const char *value = CHAR(STRING_ELT(x, 0));
  if (strcmp(value, "conserv") == 0) return 0;
  if (strcmp(value, "standard") == 0) return 1;
  if (strcmp(value, "none") == 0) return 2;

  Rf_error("'GC' must be one of 'conserv', 'standard', or 'none'");
  return -1;
}

static void check_status(int status, const char *kernel) {
  switch (status) {
    case 0:
      return;
    case 1:
      Rf_error("%s failed: incompatible input dimensions", kernel);
    case 2:
      Rf_error("%s failed: index out of range", kernel);
    case 3:
      Rf_error("%s failed: unknown genomic-control mode", kernel);
    case 4:
      Rf_error("%s failed: internal null pointer", kernel);
    case 5:
      Rf_error("%s failed: could not build native thread pool", kernel);
    case 6:
      Rf_error("%s failed: singular normal equations", kernel);
    default:
      Rf_error("%s failed: unknown native error %d", kernel, status);
  }
}

SEXP genomicssem_get_v_snp_call(
    SEXP se_snp_,
    SEXP i_ld_,
    SEXP var_snp_,
    SEXP gc_,
    SEXP coords_,
    SEXP k_,
    SEXP i_) {
  int nprotect = 0;

  SEXP se_snp = protect_real_matrix(se_snp_, "SE_SNP", &nprotect);
  SEXP i_ld = protect_real_matrix(i_ld_, "I_LD", &nprotect);
  SEXP var_snp = protect_real_vector(var_snp_, "varSNP", &nprotect);
  SEXP coords = protect_int_matrix(coords_, "coords", &nprotect);

  size_t se_nrow, se_ncol, i_ld_nrow, i_ld_ncol, coords_nrow, coords_ncol;
  matrix_dims(se_snp, "SE_SNP", &se_nrow, &se_ncol);
  matrix_dims(i_ld, "I_LD", &i_ld_nrow, &i_ld_ncol);
  matrix_dims(coords, "coords", &coords_nrow, &coords_ncol);

  int k_int = scalar_int(k_, "k");
  int i_int = scalar_int(i_, "i");
  if (k_int < 0 || i_int <= 0) {
    Rf_error("'k' must be non-negative and 'i' must be one-based");
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, k_int, k_int));
  ++nprotect;

  int status = genomicssem_fill_v_snp(
      REAL(se_snp),
      se_nrow,
      se_ncol,
      (size_t)(i_int - 1),
      REAL(i_ld),
      i_ld_nrow,
      i_ld_ncol,
      REAL(var_snp),
      (size_t)XLENGTH(var_snp),
      INTEGER(coords),
      coords_nrow,
      coords_ncol,
      (size_t)k_int,
      gc_code(gc_),
      REAL(out),
      (size_t)XLENGTH(out));

  check_status(status, "genomicssem_fill_v_snp");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_get_v_full_call(
    SEXP k_,
    SEXP v_ld_,
    SEXP var_snp_se2_,
    SEXP v_snp_) {
  int nprotect = 0;
  SEXP v_ld = protect_real_matrix(v_ld_, "V_LD", &nprotect);
  SEXP v_snp = protect_real_matrix(v_snp_, "V_SNP", &nprotect);

  size_t v_ld_nrow, v_ld_ncol, v_snp_nrow, v_snp_ncol;
  matrix_dims(v_ld, "V_LD", &v_ld_nrow, &v_ld_ncol);
  matrix_dims(v_snp, "V_SNP", &v_snp_nrow, &v_snp_ncol);

  int k_int = scalar_int(k_, "k");
  if (k_int < 0) {
    Rf_error("'k' must be non-negative");
  }

  size_t n = ((size_t)k_int + 1) * ((size_t)k_int + 2) / 2;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int)n, (int)n));
  ++nprotect;

  int status = genomicssem_fill_v_full(
      (size_t)k_int,
      REAL(v_ld),
      v_ld_nrow,
      v_ld_ncol,
      scalar_real(var_snp_se2_, "varSNPSE2"),
      REAL(v_snp),
      v_snp_nrow,
      v_snp_ncol,
      REAL(out),
      (size_t)XLENGTH(out));

  check_status(status, "genomicssem_fill_v_full");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_get_v_snp_batch_call(
    SEXP se_snp_,
    SEXP i_ld_,
    SEXP var_snp_,
    SEXP gc_,
    SEXP coords_,
    SEXP k_,
    SEXP n_threads_) {
  int nprotect = 0;

  SEXP se_snp = protect_real_matrix(se_snp_, "SE_SNP", &nprotect);
  SEXP i_ld = protect_real_matrix(i_ld_, "I_LD", &nprotect);
  SEXP var_snp = protect_real_vector(var_snp_, "varSNP", &nprotect);
  SEXP coords = protect_int_matrix(coords_, "coords", &nprotect);

  size_t se_nrow, se_ncol, i_ld_nrow, i_ld_ncol, coords_nrow, coords_ncol;
  matrix_dims(se_snp, "SE_SNP", &se_nrow, &se_ncol);
  matrix_dims(i_ld, "I_LD", &i_ld_nrow, &i_ld_ncol);
  matrix_dims(coords, "coords", &coords_nrow, &coords_ncol);

  int k_int = scalar_int(k_, "k");
  int n_threads_int = scalar_int(n_threads_, "n_threads");
  if (k_int < 0 || n_threads_int <= 0) {
    Rf_error("'k' must be non-negative and 'n_threads' must be positive");
  }

  SEXP dims = PROTECT(Rf_allocVector(INTSXP, 3));
  ++nprotect;
  INTEGER(dims)[0] = k_int;
  INTEGER(dims)[1] = k_int;
  INTEGER(dims)[2] = (int)se_nrow;

  SEXP out = PROTECT(Rf_allocArray(REALSXP, dims));
  ++nprotect;

  int status = genomicssem_fill_v_snp_batch(
      REAL(se_snp),
      se_nrow,
      se_ncol,
      REAL(i_ld),
      i_ld_nrow,
      i_ld_ncol,
      REAL(var_snp),
      (size_t)XLENGTH(var_snp),
      INTEGER(coords),
      coords_nrow,
      coords_ncol,
      (size_t)k_int,
      gc_code(gc_),
      (size_t)n_threads_int,
      REAL(out),
      (size_t)XLENGTH(out));

  check_status(status, "genomicssem_fill_v_snp_batch");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_get_s_full_call(
    SEXP k_,
    SEXP s_ld_,
    SEXP var_snp_,
    SEXP beta_snp_,
    SEXP i_) {
  int nprotect = 0;
  SEXP s_ld = protect_real_matrix(s_ld_, "S_LD", &nprotect);
  SEXP var_snp = protect_real_vector(var_snp_, "varSNP", &nprotect);
  SEXP beta_snp = protect_real_matrix(beta_snp_, "beta_SNP", &nprotect);

  size_t s_ld_nrow, s_ld_ncol, beta_nrow, beta_ncol;
  matrix_dims(s_ld, "S_LD", &s_ld_nrow, &s_ld_ncol);
  matrix_dims(beta_snp, "beta_SNP", &beta_nrow, &beta_ncol);

  int k_int = scalar_int(k_, "k");
  int i_int = scalar_int(i_, "i");
  if (k_int < 0 || i_int <= 0) {
    Rf_error("'k' must be non-negative and 'i' must be one-based");
  }

  int n = k_int + 1;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  ++nprotect;

  int status = genomicssem_fill_s_full(
      (size_t)k_int,
      REAL(s_ld),
      s_ld_nrow,
      s_ld_ncol,
      REAL(var_snp),
      (size_t)XLENGTH(var_snp),
      REAL(beta_snp),
      beta_nrow,
      beta_ncol,
      (size_t)(i_int - 1),
      REAL(out),
      (size_t)XLENGTH(out));

  check_status(status, "genomicssem_fill_s_full");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_get_z_pre_call(
    SEXP i_,
    SEXP beta_snp_,
    SEXP se_snp_,
    SEXP i_ld_,
    SEXP gc_) {
  int nprotect = 0;
  SEXP beta_snp = protect_real_matrix(beta_snp_, "beta_SNP", &nprotect);
  SEXP se_snp = protect_real_matrix(se_snp_, "SE_SNP", &nprotect);
  SEXP i_ld = protect_real_matrix(i_ld_, "I_LD", &nprotect);

  size_t beta_nrow, beta_ncol, se_nrow, se_ncol, i_ld_nrow, i_ld_ncol;
  matrix_dims(beta_snp, "beta_SNP", &beta_nrow, &beta_ncol);
  matrix_dims(se_snp, "SE_SNP", &se_nrow, &se_ncol);
  matrix_dims(i_ld, "I_LD", &i_ld_nrow, &i_ld_ncol);

  int i_int = scalar_int(i_, "i");
  if (i_int <= 0) {
    Rf_error("'i' must be one-based");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)beta_ncol));
  ++nprotect;

  int status = genomicssem_fill_z_pre(
      REAL(beta_snp),
      beta_nrow,
      beta_ncol,
      REAL(se_snp),
      se_nrow,
      se_ncol,
      REAL(i_ld),
      i_ld_nrow,
      i_ld_ncol,
      (size_t)(i_int - 1),
      gc_code(gc_),
      REAL(out),
      (size_t)XLENGTH(out));

  check_status(status, "genomicssem_fill_z_pre");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_fit_commonfactor_main_call(
    SEXP k_,
    SEXP s_full_,
    SEXP v_full_reorder_,
    SEXP w_diag_,
    SEXP start_,
    SEXP max_iter_,
    SEXP tol_) {
  int nprotect = 0;
  SEXP s_full = protect_real_matrix(s_full_, "S_Fullrun", &nprotect);
  SEXP v_full_reorder = protect_real_matrix(v_full_reorder_, "V_Full_Reorder", &nprotect);
  SEXP w_diag = protect_real_vector(w_diag_, "W_diag", &nprotect);
  SEXP start = protect_real_vector(start_, "start", &nprotect);

  size_t s_nrow, s_ncol, v_nrow, v_ncol;
  matrix_dims(s_full, "S_Fullrun", &s_nrow, &s_ncol);
  matrix_dims(v_full_reorder, "V_Full_Reorder", &v_nrow, &v_ncol);

  int k_int = scalar_int(k_, "k");
  int max_iter_int = scalar_int(max_iter_, "max_iter");
  if (k_int <= 0 || max_iter_int <= 0) {
    Rf_error("'k' and 'max_iter' must be positive");
  }

  size_t q = 2 * (size_t)k_int + 2;
  size_t out_len = 2 * q + 3;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)out_len));
  ++nprotect;

  int status = genomicssem_fit_commonfactor_main(
      (size_t)k_int,
      REAL(s_full),
      s_nrow,
      s_ncol,
      REAL(v_full_reorder),
      v_nrow,
      v_ncol,
      REAL(w_diag),
      (size_t)XLENGTH(w_diag),
      REAL(start),
      (size_t)XLENGTH(start),
      (size_t)max_iter_int,
      scalar_real(tol_, "tol"),
      REAL(out),
      out_len);

  check_status(status, "genomicssem_fit_commonfactor_main");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_fit_commonfactor_q_call(
    SEXP k_,
    SEXP s_full_,
    SEXP v_full_reorder_,
    SEXP w_diag_,
    SEXP fixed_,
    SEXP start_,
    SEXP max_iter_,
    SEXP tol_) {
  int nprotect = 0;
  SEXP s_full = protect_real_matrix(s_full_, "S_Fullrun", &nprotect);
  SEXP v_full_reorder = protect_real_matrix(v_full_reorder_, "V_Full_Reorder", &nprotect);
  SEXP w_diag = protect_real_vector(w_diag_, "W_diag", &nprotect);
  SEXP fixed = protect_real_vector(fixed_, "fixed", &nprotect);
  SEXP start = protect_real_vector(start_, "start", &nprotect);

  size_t s_nrow, s_ncol, v_nrow, v_ncol;
  matrix_dims(s_full, "S_Fullrun", &s_nrow, &s_ncol);
  matrix_dims(v_full_reorder, "V_Full_Reorder", &v_nrow, &v_ncol);

  int k_int = scalar_int(k_, "k");
  int max_iter_int = scalar_int(max_iter_, "max_iter");
  if (k_int <= 0 || max_iter_int <= 0) {
    Rf_error("'k' and 'max_iter' must be positive");
  }

  size_t q = 2 * (size_t)k_int;
  size_t out_len = q + (size_t)k_int * (size_t)k_int + 3;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)out_len));
  ++nprotect;

  int status = genomicssem_fit_commonfactor_q(
      (size_t)k_int,
      REAL(s_full),
      s_nrow,
      s_ncol,
      REAL(v_full_reorder),
      v_nrow,
      v_ncol,
      REAL(w_diag),
      (size_t)XLENGTH(w_diag),
      REAL(fixed),
      (size_t)XLENGTH(fixed),
      REAL(start),
      (size_t)XLENGTH(start),
      (size_t)max_iter_int,
      scalar_real(tol_, "tol"),
      REAL(out),
      out_len);

  check_status(status, "genomicssem_fit_commonfactor_q");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_fit_commonfactor_batch_call(
    SEXP k_,
    SEXP s_ld_,
    SEXP v_ld_,
    SEXP i_ld_,
    SEXP beta_snp_,
    SEXP se_snp_,
    SEXP var_snp_,
    SEXP gc_,
    SEXP coords_,
    SEXP var_snp_se2_,
    SEXP start_,
    SEXP max_iter_main_,
    SEXP max_iter_q_,
    SEXP tol_,
    SEXP n_threads_) {
  int nprotect = 0;
  SEXP s_ld = protect_real_matrix(s_ld_, "S_LD", &nprotect);
  SEXP v_ld = protect_real_matrix(v_ld_, "V_LD", &nprotect);
  SEXP i_ld = protect_real_matrix(i_ld_, "I_LD", &nprotect);
  SEXP beta_snp = protect_real_matrix(beta_snp_, "beta_SNP", &nprotect);
  SEXP se_snp = protect_real_matrix(se_snp_, "SE_SNP", &nprotect);
  SEXP var_snp = protect_real_vector(var_snp_, "varSNP", &nprotect);
  SEXP coords = protect_int_matrix(coords_, "coords", &nprotect);
  SEXP start = protect_real_vector(start_, "start", &nprotect);

  size_t s_ld_nrow, s_ld_ncol, v_ld_nrow, v_ld_ncol, i_ld_nrow, i_ld_ncol;
  size_t beta_nrow, beta_ncol, se_nrow, se_ncol, coords_nrow, coords_ncol;
  matrix_dims(s_ld, "S_LD", &s_ld_nrow, &s_ld_ncol);
  matrix_dims(v_ld, "V_LD", &v_ld_nrow, &v_ld_ncol);
  matrix_dims(i_ld, "I_LD", &i_ld_nrow, &i_ld_ncol);
  matrix_dims(beta_snp, "beta_SNP", &beta_nrow, &beta_ncol);
  matrix_dims(se_snp, "SE_SNP", &se_nrow, &se_ncol);
  matrix_dims(coords, "coords", &coords_nrow, &coords_ncol);

  int k_int = scalar_int(k_, "k");
  int max_iter_main_int = scalar_int(max_iter_main_, "max_iter_main");
  int max_iter_q_int = scalar_int(max_iter_q_, "max_iter_q");
  int n_threads_int = scalar_int(n_threads_, "n_threads");
  if (k_int <= 0 || max_iter_main_int <= 0 || max_iter_q_int <= 0 || n_threads_int <= 0) {
    Rf_error("'k', iteration counts, and 'n_threads' must be positive");
  }

  const size_t out_cols = 7;
  size_t out_len = beta_nrow * out_cols;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)out_len));
  ++nprotect;

  int status = genomicssem_fit_commonfactor_batch(
      (size_t)k_int,
      REAL(s_ld),
      s_ld_nrow,
      s_ld_ncol,
      REAL(v_ld),
      v_ld_nrow,
      v_ld_ncol,
      REAL(i_ld),
      i_ld_nrow,
      i_ld_ncol,
      REAL(beta_snp),
      beta_nrow,
      beta_ncol,
      REAL(se_snp),
      se_nrow,
      se_ncol,
      REAL(var_snp),
      (size_t)XLENGTH(var_snp),
      INTEGER(coords),
      coords_nrow,
      coords_ncol,
      scalar_real(var_snp_se2_, "varSNPSE2"),
      gc_code(gc_),
      REAL(start),
      (size_t)XLENGTH(start),
      (size_t)max_iter_main_int,
      (size_t)max_iter_q_int,
      scalar_real(tol_, "tol"),
      (size_t)n_threads_int,
      REAL(out),
      out_len);

  check_status(status, "genomicssem_fit_commonfactor_batch");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_fit_generic_sem_call(
    SEXP obs_n_,
    SEXP total_n_,
    SEXP s_full_,
    SEXP v_full_reorder_,
    SEXP w_diag_,
    SEXP b_fixed_,
    SEXP psi_fixed_,
    SEXP b_free_,
    SEXP psi_free_,
    SEXP start_,
    SEXP max_iter_,
    SEXP tol_) {
  int nprotect = 0;
  SEXP s_full = protect_real_matrix(s_full_, "S_Fullrun", &nprotect);
  SEXP v_full_reorder = protect_real_matrix(v_full_reorder_, "V_Full_Reorder", &nprotect);
  SEXP w_diag = protect_real_vector(w_diag_, "W_diag", &nprotect);
  SEXP b_fixed = protect_real_matrix(b_fixed_, "B_fixed", &nprotect);
  SEXP psi_fixed = protect_real_matrix(psi_fixed_, "Psi_fixed", &nprotect);
  SEXP b_free = protect_int_matrix(b_free_, "B_free", &nprotect);
  SEXP psi_free = protect_int_matrix(psi_free_, "Psi_free", &nprotect);
  SEXP start = protect_real_vector(start_, "start", &nprotect);

  size_t s_nrow, s_ncol, v_nrow, v_ncol, b_nrow, b_ncol, psi_nrow, psi_ncol;
  size_t b_free_nrow, b_free_ncol, psi_free_nrow, psi_free_ncol;
  matrix_dims(s_full, "S_Fullrun", &s_nrow, &s_ncol);
  matrix_dims(v_full_reorder, "V_Full_Reorder", &v_nrow, &v_ncol);
  matrix_dims(b_fixed, "B_fixed", &b_nrow, &b_ncol);
  matrix_dims(psi_fixed, "Psi_fixed", &psi_nrow, &psi_ncol);
  matrix_dims(b_free, "B_free", &b_free_nrow, &b_free_ncol);
  matrix_dims(psi_free, "Psi_free", &psi_free_nrow, &psi_free_ncol);

  int obs_n_int = scalar_int(obs_n_, "obs_n");
  int total_n_int = scalar_int(total_n_, "total_n");
  int max_iter_int = scalar_int(max_iter_, "max_iter");
  if (obs_n_int <= 0 || total_n_int <= 0 || max_iter_int <= 0) {
    Rf_error("'obs_n', 'total_n', and 'max_iter' must be positive");
  }

  size_t q = (size_t)XLENGTH(start);
  size_t out_len = 2 * q + (size_t)obs_n_int * (size_t)obs_n_int + 3;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)out_len));
  ++nprotect;

  int status = genomicssem_fit_generic_sem(
      (size_t)obs_n_int,
      (size_t)total_n_int,
      REAL(s_full),
      s_nrow,
      s_ncol,
      REAL(v_full_reorder),
      v_nrow,
      v_ncol,
      REAL(w_diag),
      (size_t)XLENGTH(w_diag),
      REAL(b_fixed),
      (size_t)XLENGTH(b_fixed),
      REAL(psi_fixed),
      (size_t)XLENGTH(psi_fixed),
      INTEGER(b_free),
      (size_t)XLENGTH(b_free),
      INTEGER(psi_free),
      (size_t)XLENGTH(psi_free),
      REAL(start),
      q,
      (size_t)max_iter_int,
      scalar_real(tol_, "tol"),
      REAL(out),
      out_len);

  check_status(status, "genomicssem_fit_generic_sem");
  UNPROTECT(nprotect);
  return out;
}
