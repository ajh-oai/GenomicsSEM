library(GenomicSEM)

genomicssem_ns <- asNamespace("GenomicSEM")
.get_V_SNP <- get(".get_V_SNP", envir = genomicssem_ns)
.get_V_SNP_r <- get(".get_V_SNP_r", envir = genomicssem_ns)
.get_V_SNP_batch <- get(".get_V_SNP_batch", envir = genomicssem_ns)
.get_V_full <- get(".get_V_full", envir = genomicssem_ns)
.get_V_full_r <- get(".get_V_full_r", envir = genomicssem_ns)
.get_S_Full <- get(".get_S_Full", envir = genomicssem_ns)
.get_S_Full_r <- get(".get_S_Full_r", envir = genomicssem_ns)
.get_Z_pre <- get(".get_Z_pre", envir = genomicssem_ns)
.get_Z_pre_r <- get(".get_Z_pre_r", envir = genomicssem_ns)

stopifnot_rust <- function(x) {
  if (!isTRUE(x)) {
    stop("assertion failed", call. = FALSE)
  }
}

stopifnot_rust(is.loaded("genomicssem_get_v_snp_call"))

withr_seed <- function(seed, code) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  set.seed(seed)
  on.exit({
    if (is.null(old)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  })
  force(code)
}

run_case <- function(seed, k, n_snp, gc_mode) {
  withr_seed(seed, {
    beta_SNP <- matrix(rnorm(n_snp * k), n_snp, k)
    SE_SNP <- matrix(runif(n_snp * k, 0.01, 0.5), n_snp, k)

    I_LD <- matrix(runif(k * k, 0.05, 0.4), k, k)
    I_LD <- 0.5 * (I_LD + t(I_LD))
    diag(I_LD) <- runif(k, 1.0, 1.8)

    S_LD <- matrix(runif(k * k, -0.2, 0.6), k, k)
    S_LD <- 0.5 * (S_LD + t(S_LD))
    diag(S_LD) <- runif(k, 0.8, 1.4)
    colnames(S_LD) <- paste0("trait", seq_len(k))

    v_dim <- k * (k + 1) / 2
    V_LD <- matrix(runif(v_dim * v_dim, 0.01, 0.5), v_dim, v_dim)
    V_LD <- 0.5 * (V_LD + t(V_LD))

    varSNP <- runif(n_snp, 0.05, 0.95)
    varSNPSE2 <- runif(1, 0.001, 0.1)
    coords <- which(!is.na(I_LD), arr.ind = TRUE)

    rust_batch <- .get_V_SNP_batch(SE_SNP, I_LD, varSNP, gc_mode, coords, k, n_threads = 2L)

    for (i in seq_len(n_snp)) {
      rust_v_snp <- .get_V_SNP(SE_SNP, I_LD, varSNP, gc_mode, coords, k, i)
      ref_v_snp <- .get_V_SNP_r(SE_SNP, I_LD, varSNP, gc_mode, coords, k, i)
      stopifnot_rust(identical(rust_v_snp, ref_v_snp))
      stopifnot_rust(identical(rust_batch[, , i], ref_v_snp))

      rust_v_full <- .get_V_full(k, V_LD, varSNPSE2, rust_v_snp)
      ref_v_full <- .get_V_full_r(k, V_LD, varSNPSE2, ref_v_snp)
      stopifnot_rust(identical(rust_v_full, ref_v_full))

      rust_s_full <- .get_S_Full(k, S_LD, varSNP, beta_SNP, FALSE, i)
      ref_s_full <- .get_S_Full_r(k, S_LD, varSNP, beta_SNP, FALSE, i)
      stopifnot_rust(identical(rust_s_full, ref_s_full))

      rust_z_pre <- .get_Z_pre(i, beta_SNP, SE_SNP, I_LD, gc_mode)
      ref_z_pre <- .get_Z_pre_r(i, beta_SNP, SE_SNP, I_LD, gc_mode)
      stopifnot_rust(identical(rust_z_pre, ref_z_pre))
    }
  })
}

for (gc_mode in c("conserv", "standard", "none")) {
  for (seed in 1:50) {
    run_case(seed = seed, k = 2 + seed %% 5, n_snp = 1 + seed %% 7, gc_mode = gc_mode)
  }
}

cat("Rust kernel parity tests passed\n")
