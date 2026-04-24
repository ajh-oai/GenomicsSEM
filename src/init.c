#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP genomicssem_get_v_snp_call(
    SEXP se_snp_,
    SEXP i_ld_,
    SEXP var_snp_,
    SEXP gc_,
    SEXP coords_,
    SEXP k_,
    SEXP i_);
extern SEXP genomicssem_get_v_full_call(
    SEXP k_,
    SEXP v_ld_,
    SEXP var_snp_se2_,
    SEXP v_snp_);
extern SEXP genomicssem_get_v_snp_batch_call(
    SEXP se_snp_,
    SEXP i_ld_,
    SEXP var_snp_,
    SEXP gc_,
    SEXP coords_,
    SEXP k_,
    SEXP n_threads_);
extern SEXP genomicssem_get_s_full_call(
    SEXP k_,
    SEXP s_ld_,
    SEXP var_snp_,
    SEXP beta_snp_,
    SEXP i_);
extern SEXP genomicssem_get_z_pre_call(
    SEXP i_,
    SEXP beta_snp_,
    SEXP se_snp_,
    SEXP i_ld_,
    SEXP gc_);
extern SEXP genomicssem_fit_commonfactor_main_call(
    SEXP k_,
    SEXP s_full_,
    SEXP v_full_reorder_,
    SEXP w_diag_,
    SEXP start_,
    SEXP max_iter_,
    SEXP tol_);
extern SEXP genomicssem_fit_commonfactor_q_call(
    SEXP k_,
    SEXP s_full_,
    SEXP v_full_reorder_,
    SEXP w_diag_,
    SEXP fixed_,
    SEXP start_,
    SEXP max_iter_,
    SEXP tol_);
extern SEXP genomicssem_fit_commonfactor_batch_call(
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
    SEXP n_threads_);
extern SEXP genomicssem_fit_generic_sem_call(
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
    SEXP tol_);

static const R_CallMethodDef CallEntries[] = {
    {"genomicssem_get_v_snp_call", (DL_FUNC)&genomicssem_get_v_snp_call, 7},
    {"genomicssem_get_v_full_call", (DL_FUNC)&genomicssem_get_v_full_call, 4},
    {"genomicssem_get_v_snp_batch_call", (DL_FUNC)&genomicssem_get_v_snp_batch_call, 7},
    {"genomicssem_get_s_full_call", (DL_FUNC)&genomicssem_get_s_full_call, 5},
    {"genomicssem_get_z_pre_call", (DL_FUNC)&genomicssem_get_z_pre_call, 5},
    {"genomicssem_fit_commonfactor_main_call", (DL_FUNC)&genomicssem_fit_commonfactor_main_call, 7},
    {"genomicssem_fit_commonfactor_q_call", (DL_FUNC)&genomicssem_fit_commonfactor_q_call, 8},
    {"genomicssem_fit_commonfactor_batch_call", (DL_FUNC)&genomicssem_fit_commonfactor_batch_call, 15},
    {"genomicssem_fit_generic_sem_call", (DL_FUNC)&genomicssem_fit_generic_sem_call, 12},
    {NULL, NULL, 0}
};

void R_init_GenomicSEM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
