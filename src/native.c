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

static SEXP protect_int_vector(SEXP x, const char *name, int *nprotect) {
  if (!Rf_isInteger(x) || Rf_isMatrix(x)) {
    if (!Rf_isNumeric(x) || Rf_isMatrix(x)) {
      Rf_error("'%s' must be an integer vector", name);
    }
    ++(*nprotect);
    return PROTECT(Rf_coerceVector(x, INTSXP));
  }
  return x;
}

static SEXP protect_string_vector(SEXP x, const char *name, int *nprotect) {
  if (TYPEOF(x) == STRSXP && !Rf_isMatrix(x)) {
    return x;
  }
  ++(*nprotect);
  SEXP out = PROTECT(Rf_coerceVector(x, STRSXP));
  if (Rf_isMatrix(out)) {
    Rf_error("'%s' must be a character vector", name);
  }
  return out;
}

static const char *scalar_string(SEXP x, const char *name) {
  if (TYPEOF(x) != STRSXP || XLENGTH(x) != 1 || STRING_ELT(x, 0) == NA_STRING) {
    Rf_error("'%s' must be a non-missing character scalar", name);
  }
  return CHAR(STRING_ELT(x, 0));
}

static const char **string_ptrs(SEXP x) {
  R_xlen_t n = XLENGTH(x);
  const char **out = (const char **)R_alloc((size_t)n, sizeof(const char *));
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP value = STRING_ELT(x, i);
    out[i] = value == NA_STRING ? "" : CHAR(value);
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

static int scalar_logical(SEXP x, const char *name) {
  int value = Rf_asLogical(x);
  if (value == NA_LOGICAL) {
    Rf_error("'%s' must be TRUE or FALSE", name);
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

SEXP genomicssem_fit_generic_sem_batch_call(
    SEXP obs_n_,
    SEXP total_n_,
    SEXP s_ld_,
    SEXP v_ld_,
    SEXP i_ld_,
    SEXP beta_snp_,
    SEXP se_snp_,
    SEXP var_snp_,
    SEXP gc_,
    SEXP coords_,
    SEXP var_snp_se2_,
    SEXP order_,
    SEXP spec_to_original_,
    SEXP b_fixed_,
    SEXP psi_fixed_,
    SEXP b_free_,
    SEXP psi_free_,
    SEXP start_,
    SEXP q_snp_indices_,
    SEXP q_snp_lengths_,
    SEXP max_iter_,
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
  SEXP order = protect_int_matrix(order_, "order", &nprotect);
  SEXP spec_to_original = protect_int_matrix(spec_to_original_, "spec_to_original", &nprotect);
  SEXP b_fixed = protect_real_matrix(b_fixed_, "B_fixed", &nprotect);
  SEXP psi_fixed = protect_real_matrix(psi_fixed_, "Psi_fixed", &nprotect);
  SEXP b_free = protect_int_matrix(b_free_, "B_free", &nprotect);
  SEXP psi_free = protect_int_matrix(psi_free_, "Psi_free", &nprotect);
  SEXP start = protect_real_vector(start_, "start", &nprotect);
  SEXP q_snp_indices = protect_int_matrix(q_snp_indices_, "q_snp_indices", &nprotect);
  SEXP q_snp_lengths = protect_int_matrix(q_snp_lengths_, "q_snp_lengths", &nprotect);

  size_t s_ld_nrow, s_ld_ncol, v_ld_nrow, v_ld_ncol, i_ld_nrow, i_ld_ncol;
  size_t beta_nrow, beta_ncol, se_nrow, se_ncol, coords_nrow, coords_ncol;
  size_t order_nrow, order_ncol, spec_nrow, spec_ncol, b_nrow, b_ncol, psi_nrow, psi_ncol;
  size_t b_free_nrow, b_free_ncol, psi_free_nrow, psi_free_ncol, q_idx_nrow, q_idx_ncol;
  size_t q_len_nrow, q_len_ncol;
  matrix_dims(s_ld, "S_LD", &s_ld_nrow, &s_ld_ncol);
  matrix_dims(v_ld, "V_LD", &v_ld_nrow, &v_ld_ncol);
  matrix_dims(i_ld, "I_LD", &i_ld_nrow, &i_ld_ncol);
  matrix_dims(beta_snp, "beta_SNP", &beta_nrow, &beta_ncol);
  matrix_dims(se_snp, "SE_SNP", &se_nrow, &se_ncol);
  matrix_dims(coords, "coords", &coords_nrow, &coords_ncol);
  matrix_dims(order, "order", &order_nrow, &order_ncol);
  matrix_dims(spec_to_original, "spec_to_original", &spec_nrow, &spec_ncol);
  matrix_dims(b_fixed, "B_fixed", &b_nrow, &b_ncol);
  matrix_dims(psi_fixed, "Psi_fixed", &psi_nrow, &psi_ncol);
  matrix_dims(b_free, "B_free", &b_free_nrow, &b_free_ncol);
  matrix_dims(psi_free, "Psi_free", &psi_free_nrow, &psi_free_ncol);
  matrix_dims(q_snp_indices, "q_snp_indices", &q_idx_nrow, &q_idx_ncol);
  matrix_dims(q_snp_lengths, "q_snp_lengths", &q_len_nrow, &q_len_ncol);

  int obs_n_int = scalar_int(obs_n_, "obs_n");
  int total_n_int = scalar_int(total_n_, "total_n");
  int max_iter_int = scalar_int(max_iter_, "max_iter");
  int n_threads_int = scalar_int(n_threads_, "n_threads");
  if (obs_n_int <= 0 || total_n_int <= 0 || max_iter_int <= 0 || n_threads_int <= 0) {
    Rf_error("'obs_n', 'total_n', 'max_iter', and 'n_threads' must be positive");
  }

  size_t q = (size_t)XLENGTH(start);
  size_t out_cols = 2 * q + 1 + q_idx_ncol + 2;
  size_t out_len = beta_nrow * out_cols;
  SEXP out = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)out_len));
  ++nprotect;

  int status = genomicssem_fit_generic_sem_batch(
      (size_t)obs_n_int,
      (size_t)total_n_int,
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
      INTEGER(order),
      (size_t)XLENGTH(order),
      INTEGER(spec_to_original),
      (size_t)XLENGTH(spec_to_original),
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
      INTEGER(q_snp_indices),
      q_idx_nrow,
      q_idx_ncol,
      INTEGER(q_snp_lengths),
      (size_t)XLENGTH(q_snp_lengths),
      (size_t)max_iter_int,
      scalar_real(tol_, "tol"),
      (size_t)n_threads_int,
      REAL(out),
      out_len);

  check_status(status, "genomicssem_fit_generic_sem_batch");
  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_munge_qc_call(
    SEXP a1_ref_,
    SEXP a2_ref_,
    SEXP a1_file_,
    SEXP a2_file_,
    SEXP effect_,
    SEXP p_,
    SEXP info_,
    SEXP maf_,
    SEXP info_filter_,
    SEXP maf_filter_) {
  int nprotect = 0;
  SEXP a1_ref = protect_int_vector(a1_ref_, "A1.x", &nprotect);
  SEXP a2_ref = protect_int_vector(a2_ref_, "A2.x", &nprotect);
  SEXP a1_file = protect_int_vector(a1_file_, "A1.y", &nprotect);
  SEXP a2_file = protect_int_vector(a2_file_, "A2.y", &nprotect);
  SEXP effect = protect_real_vector(effect_, "effect", &nprotect);
  SEXP p = protect_real_vector(p_, "P", &nprotect);
  SEXP info = protect_real_vector(info_, "INFO", &nprotect);
  SEXP maf = protect_real_vector(maf_, "MAF", &nprotect);

  R_xlen_t n_xlen = XLENGTH(effect);
  if (XLENGTH(a1_ref) != n_xlen || XLENGTH(a2_ref) != n_xlen ||
      XLENGTH(a1_file) != n_xlen || XLENGTH(a2_file) != n_xlen || XLENGTH(p) != n_xlen) {
    Rf_error("munge QC inputs must have the same length");
  }
  if (XLENGTH(info) != 0 && XLENGTH(info) != n_xlen) {
    Rf_error("'INFO' must be empty or match 'effect' length");
  }
  if (XLENGTH(maf) != 0 && XLENGTH(maf) != n_xlen) {
    Rf_error("'MAF' must be empty or match 'effect' length");
  }

  SEXP keep = PROTECT(Rf_allocVector(INTSXP, n_xlen));
  SEXP z = PROTECT(Rf_allocVector(REALSXP, n_xlen));
  SEXP counts = PROTECT(Rf_allocVector(INTSXP, 8));
  nprotect += 3;

  size_t out_count = 0;
  int status = genomicssem_munge_qc(
      INTEGER(a1_ref),
      INTEGER(a2_ref),
      INTEGER(a1_file),
      INTEGER(a2_file),
      REAL(effect),
      REAL(p),
      (size_t)n_xlen,
      XLENGTH(info) == 0 ? NULL : REAL(info),
      (size_t)XLENGTH(info),
      XLENGTH(maf) == 0 ? NULL : REAL(maf),
      (size_t)XLENGTH(maf),
      scalar_real(info_filter_, "info.filter"),
      scalar_real(maf_filter_, "maf.filter"),
      INTEGER(keep),
      REAL(z),
      (size_t)n_xlen,
      INTEGER(counts),
      8,
      &out_count);

  check_status(status, "genomicssem_munge_qc");

  SEXP out_count_sexp = PROTECT(Rf_ScalarInteger((int)out_count));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  nprotect += 3;
  SET_VECTOR_ELT(out, 0, keep);
  SET_VECTOR_ELT(out, 1, z);
  SET_VECTOR_ELT(out, 2, counts);
  SET_VECTOR_ELT(out, 3, out_count_sexp);
  SET_STRING_ELT(names, 0, Rf_mkChar("keep"));
  SET_STRING_ELT(names, 1, Rf_mkChar("z"));
  SET_STRING_ELT(names, 2, Rf_mkChar("counts"));
  SET_STRING_ELT(names, 3, Rf_mkChar("n"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_munge_fused_call(
    SEXP filename_,
    SEXP output_path_,
    SEXP ref_snp_,
    SEXP ref_a1_,
    SEXP ref_a2_,
    SEXP col_indices_,
    SEXP provided_n_,
    SEXP n_multiplier_,
    SEXP info_filter_,
    SEXP maf_filter_) {
  int nprotect = 0;
  SEXP ref_snp = protect_string_vector(ref_snp_, "ref$SNP", &nprotect);
  SEXP ref_a1 = protect_int_vector(ref_a1_, "ref$A1", &nprotect);
  SEXP ref_a2 = protect_int_vector(ref_a2_, "ref$A2", &nprotect);
  SEXP col_indices = protect_int_vector(col_indices_, "col_indices", &nprotect);

  R_xlen_t ref_len = XLENGTH(ref_snp);
  if (XLENGTH(ref_a1) != ref_len || XLENGTH(ref_a2) != ref_len) {
    Rf_error("reference SNP and allele vectors must have the same length");
  }

  SEXP counts = PROTECT(Rf_allocVector(INTSXP, 8));
  ++nprotect;

  size_t rows_total = 0;
  size_t rows_joined = 0;
  size_t rows_written = 0;
  int unsupported = 0;
  int status = genomicssem_munge_fused(
      scalar_string(filename_, "filename"),
      scalar_string(output_path_, "output_path"),
      string_ptrs(ref_snp),
      (size_t)ref_len,
      INTEGER(ref_a1),
      (size_t)XLENGTH(ref_a1),
      INTEGER(ref_a2),
      (size_t)XLENGTH(ref_a2),
      INTEGER(col_indices),
      (size_t)XLENGTH(col_indices),
      Rf_asReal(provided_n_),
      scalar_real(n_multiplier_, "n_multiplier"),
      scalar_real(info_filter_, "info.filter"),
      scalar_real(maf_filter_, "maf.filter"),
      INTEGER(counts),
      8,
      &rows_total,
      &rows_joined,
      &rows_written,
      &unsupported);

  check_status(status, "genomicssem_munge_fused");
  if (unsupported) {
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  SEXP rows = PROTECT(Rf_allocVector(REALSXP, 3));
  nprotect += 3;

  REAL(rows)[0] = (double)rows_total;
  REAL(rows)[1] = (double)rows_joined;
  REAL(rows)[2] = (double)rows_written;
  SET_VECTOR_ELT(out, 0, counts);
  SET_VECTOR_ELT(out, 1, rows);
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger((int)rows_written));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(1));
  SET_STRING_ELT(names, 0, Rf_mkChar("counts"));
  SET_STRING_ELT(names, 1, Rf_mkChar("rows"));
  SET_STRING_ELT(names, 2, Rf_mkChar("n"));
  SET_STRING_ELT(names, 3, Rf_mkChar("used"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_sumstats_qc_call(
    SEXP a1_ref_,
    SEXP a2_ref_,
    SEXP a1_file_,
    SEXP a2_file_,
    SEXP effect_,
    SEXP se_,
    SEXP p_,
    SEXP n_values_,
    SEXP maf_ref_,
    SEXP maf_file_,
    SEXP info_,
    SEXP info_filter_,
    SEXP ols_,
    SEXP beta_is_character_,
    SEXP linprob_,
    SEXP se_logit_) {
  int nprotect = 0;
  SEXP a1_ref = protect_int_vector(a1_ref_, "A1.x", &nprotect);
  SEXP a2_ref = protect_int_vector(a2_ref_, "A2.x", &nprotect);
  SEXP a1_file = protect_int_vector(a1_file_, "A1.y", &nprotect);
  SEXP a2_file = protect_int_vector(a2_file_, "A2.y", &nprotect);
  SEXP effect = protect_real_vector(effect_, "effect", &nprotect);
  SEXP se = protect_real_vector(se_, "SE", &nprotect);
  SEXP p = protect_real_vector(p_, "P", &nprotect);
  SEXP n_values = protect_real_vector(n_values_, "N", &nprotect);
  SEXP maf_ref = protect_real_vector(maf_ref_, "MAF", &nprotect);
  SEXP maf_file = protect_real_vector(maf_file_, "MAF.y", &nprotect);
  SEXP info = protect_real_vector(info_, "INFO", &nprotect);

  R_xlen_t n_xlen = XLENGTH(effect);
  if (XLENGTH(a1_ref) != n_xlen || XLENGTH(a2_ref) != n_xlen ||
      XLENGTH(a1_file) != n_xlen || XLENGTH(a2_file) != n_xlen ||
      XLENGTH(se) != n_xlen || XLENGTH(p) != n_xlen ||
      XLENGTH(n_values) != n_xlen || XLENGTH(maf_ref) != n_xlen) {
    Rf_error("sumstats QC inputs must have the same length");
  }
  if (XLENGTH(maf_file) != 0 && XLENGTH(maf_file) != n_xlen) {
    Rf_error("'MAF.y' must be empty or match 'effect' length");
  }
  if (XLENGTH(info) != 0 && XLENGTH(info) != n_xlen) {
    Rf_error("'INFO' must be empty or match 'effect' length");
  }

  SEXP keep = PROTECT(Rf_allocVector(INTSXP, n_xlen));
  SEXP beta = PROTECT(Rf_allocVector(REALSXP, n_xlen));
  SEXP se_out = PROTECT(Rf_allocVector(REALSXP, n_xlen));
  SEXP counts = PROTECT(Rf_allocVector(INTSXP, 10));
  nprotect += 4;

  size_t out_count = 0;
  int status = genomicssem_sumstats_qc(
      INTEGER(a1_ref),
      INTEGER(a2_ref),
      INTEGER(a1_file),
      INTEGER(a2_file),
      REAL(effect),
      REAL(se),
      REAL(p),
      REAL(n_values),
      REAL(maf_ref),
      (size_t)n_xlen,
      XLENGTH(maf_file) == 0 ? NULL : REAL(maf_file),
      (size_t)XLENGTH(maf_file),
      XLENGTH(info) == 0 ? NULL : REAL(info),
      (size_t)XLENGTH(info),
      scalar_real(info_filter_, "info.filter"),
      scalar_logical(ols_, "OLS"),
      scalar_logical(beta_is_character_, "beta_is_character"),
      scalar_logical(linprob_, "linprob"),
      scalar_logical(se_logit_, "se.logit"),
      INTEGER(keep),
      REAL(beta),
      REAL(se_out),
      (size_t)n_xlen,
      INTEGER(counts),
      10,
      &out_count);

  check_status(status, "genomicssem_sumstats_qc");

  SEXP out_count_sexp = PROTECT(Rf_ScalarInteger((int)out_count));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
  nprotect += 3;
  SET_VECTOR_ELT(out, 0, keep);
  SET_VECTOR_ELT(out, 1, beta);
  SET_VECTOR_ELT(out, 2, se_out);
  SET_VECTOR_ELT(out, 3, counts);
  SET_VECTOR_ELT(out, 4, out_count_sexp);
  SET_STRING_ELT(names, 0, Rf_mkChar("keep"));
  SET_STRING_ELT(names, 1, Rf_mkChar("beta"));
  SET_STRING_ELT(names, 2, Rf_mkChar("se"));
  SET_STRING_ELT(names, 3, Rf_mkChar("counts"));
  SET_STRING_ELT(names, 4, Rf_mkChar("n"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_sumstats_fused_call(
    SEXP filename_,
    SEXP ref_snp_,
    SEXP ref_a1_,
    SEXP ref_a2_,
    SEXP ref_maf_,
    SEXP col_indices_,
    SEXP provided_n_,
    SEXP info_filter_,
    SEXP ols_,
    SEXP beta_is_character_,
    SEXP linprob_,
    SEXP se_logit_) {
  int nprotect = 0;
  SEXP ref_snp = protect_string_vector(ref_snp_, "ref$SNP", &nprotect);
  SEXP ref_a1 = protect_int_vector(ref_a1_, "ref$A1", &nprotect);
  SEXP ref_a2 = protect_int_vector(ref_a2_, "ref$A2", &nprotect);
  SEXP ref_maf = protect_real_vector(ref_maf_, "ref$MAF", &nprotect);
  SEXP col_indices = protect_int_vector(col_indices_, "col_indices", &nprotect);

  R_xlen_t ref_len = XLENGTH(ref_snp);
  if (XLENGTH(ref_a1) != ref_len || XLENGTH(ref_a2) != ref_len || XLENGTH(ref_maf) != ref_len) {
    Rf_error("reference SNP, allele, and MAF vectors must have the same length");
  }

  SEXP keep = PROTECT(Rf_allocVector(INTSXP, ref_len));
  SEXP beta = PROTECT(Rf_allocVector(REALSXP, ref_len));
  SEXP se = PROTECT(Rf_allocVector(REALSXP, ref_len));
  SEXP counts = PROTECT(Rf_allocVector(INTSXP, 10));
  nprotect += 4;

  size_t rows_total = 0;
  size_t rows_duplicate_removed = 0;
  size_t rows_joined = 0;
  size_t rows_written = 0;
  int unsupported = 0;
  int status = genomicssem_sumstats_fused(
      scalar_string(filename_, "filename"),
      string_ptrs(ref_snp),
      (size_t)ref_len,
      INTEGER(ref_a1),
      (size_t)XLENGTH(ref_a1),
      INTEGER(ref_a2),
      (size_t)XLENGTH(ref_a2),
      REAL(ref_maf),
      (size_t)XLENGTH(ref_maf),
      INTEGER(col_indices),
      (size_t)XLENGTH(col_indices),
      Rf_asReal(provided_n_),
      scalar_real(info_filter_, "info.filter"),
      scalar_logical(ols_, "OLS"),
      scalar_logical(beta_is_character_, "beta_is_character"),
      scalar_logical(linprob_, "linprob"),
      scalar_logical(se_logit_, "se.logit"),
      INTEGER(keep),
      REAL(beta),
      REAL(se),
      (size_t)ref_len,
      INTEGER(counts),
      10,
      &rows_total,
      &rows_duplicate_removed,
      &rows_joined,
      &rows_written,
      &unsupported);

  check_status(status, "genomicssem_sumstats_fused");
  if (unsupported) {
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  SEXP out_count = PROTECT(Rf_ScalarInteger((int)rows_written));
  SEXP rows = PROTECT(Rf_allocVector(REALSXP, 4));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));
  nprotect += 4;

  REAL(rows)[0] = (double)rows_total;
  REAL(rows)[1] = (double)rows_duplicate_removed;
  REAL(rows)[2] = (double)rows_joined;
  REAL(rows)[3] = (double)rows_written;
  SET_VECTOR_ELT(out, 0, keep);
  SET_VECTOR_ELT(out, 1, beta);
  SET_VECTOR_ELT(out, 2, se);
  SET_VECTOR_ELT(out, 3, counts);
  SET_VECTOR_ELT(out, 4, out_count);
  SET_VECTOR_ELT(out, 5, rows);
  SET_VECTOR_ELT(out, 6, Rf_ScalarLogical(1));
  SET_STRING_ELT(names, 0, Rf_mkChar("keep"));
  SET_STRING_ELT(names, 1, Rf_mkChar("beta"));
  SET_STRING_ELT(names, 2, Rf_mkChar("se"));
  SET_STRING_ELT(names, 3, Rf_mkChar("counts"));
  SET_STRING_ELT(names, 4, Rf_mkChar("n"));
  SET_STRING_ELT(names, 5, Rf_mkChar("rows"));
  SET_STRING_ELT(names, 6, Rf_mkChar("used"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}

SEXP genomicssem_ldsc_block_products_call(
    SEXP weighted_ld_,
    SEXP weighted_chi_,
    SEXP n_blocks_,
    SEXP n_threads_) {
  int nprotect = 0;
  SEXP weighted_ld = protect_real_matrix(weighted_ld_, "weighted.LD", &nprotect);
  SEXP weighted_chi = protect_real_vector(weighted_chi_, "weighted.chi", &nprotect);

  size_t n_snps, n_annot;
  matrix_dims(weighted_ld, "weighted.LD", &n_snps, &n_annot);

  int n_blocks_int = scalar_int(n_blocks_, "n.blocks");
  if (n_blocks_int <= 0) {
    Rf_error("'n.blocks' must be positive");
  }
  if ((size_t)XLENGTH(weighted_chi) != n_snps) {
    Rf_error("'weighted.chi' must have one value per row of 'weighted.LD'");
  }
  int n_threads_int = scalar_int(n_threads_, "n.threads");
  if (n_threads_int <= 0) {
    Rf_error("'n.threads' must be positive");
  }

  SEXP xty_block = PROTECT(Rf_allocMatrix(REALSXP, n_blocks_int, (int)n_annot));
  SEXP xtx_block = PROTECT(Rf_allocMatrix(REALSXP, n_blocks_int * (int)n_annot, (int)n_annot));
  SEXP xty = PROTECT(Rf_allocMatrix(REALSXP, (int)n_annot, 1));
  SEXP xtx = PROTECT(Rf_allocMatrix(REALSXP, (int)n_annot, (int)n_annot));
  nprotect += 4;

  int status = genomicssem_ldsc_block_products(
      REAL(weighted_ld),
      n_snps,
      n_annot,
      REAL(weighted_chi),
      (size_t)XLENGTH(weighted_chi),
      (size_t)n_blocks_int,
      REAL(xty_block),
      (size_t)XLENGTH(xty_block),
      REAL(xtx_block),
      (size_t)XLENGTH(xtx_block),
      REAL(xty),
      (size_t)XLENGTH(xty),
      REAL(xtx),
      (size_t)XLENGTH(xtx),
      (size_t)n_threads_int);

  check_status(status, "genomicssem_ldsc_block_products");

  SEXP delete_from = PROTECT(Rf_allocVector(REALSXP, n_blocks_int));
  SEXP delete_to = PROTECT(Rf_allocVector(INTSXP, n_blocks_int));
  nprotect += 2;
  for (int i = 0; i < n_blocks_int; ++i) {
    REAL(delete_from)[i] = (double)(i * (int)n_annot + 1);
    INTEGER(delete_to)[i] = (i + 1) * (int)n_annot;
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
  nprotect += 2;
  SET_VECTOR_ELT(out, 0, xty_block);
  SET_VECTOR_ELT(out, 1, xtx_block);
  SET_VECTOR_ELT(out, 2, xty);
  SET_VECTOR_ELT(out, 3, xtx);
  SET_VECTOR_ELT(out, 4, delete_from);
  SET_VECTOR_ELT(out, 5, delete_to);
  SET_STRING_ELT(names, 0, Rf_mkChar("xty.block.values"));
  SET_STRING_ELT(names, 1, Rf_mkChar("xtx.block.values"));
  SET_STRING_ELT(names, 2, Rf_mkChar("xty"));
  SET_STRING_ELT(names, 3, Rf_mkChar("xtx"));
  SET_STRING_ELT(names, 4, Rf_mkChar("delete.from"));
  SET_STRING_ELT(names, 5, Rf_mkChar("delete.to"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}
