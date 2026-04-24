library(GenomicSEM)

make_inputs <- function(n_snp = 12L, k = 5L) {
  set.seed(29)
  traits <- paste0("T", seq_len(k))

  a <- matrix(rnorm(k * k), k, k)
  S_LD <- cov2cor(crossprod(a) + diag(k))
  dimnames(S_LD) <- list(traits, traits)

  V_LD <- diag(runif(k * (k + 1) / 2, 1e-4, 5e-4))

  I_LD <- matrix(0.02, k, k)
  diag(I_LD) <- runif(k, 1.0, 1.2)
  dimnames(I_LD) <- list(traits, traits)

  SNPs <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    CHR = rep(1L, n_snp),
    BP = seq_len(n_snp),
    MAF = runif(n_snp, 0.05, 0.5),
    A1 = rep("A", n_snp),
    A2 = rep("G", n_snp),
    check.names = FALSE
  )

  beta <- as.data.frame(matrix(rnorm(n_snp * k, sd = 0.03), n_snp, k), check.names = FALSE)
  names(beta) <- paste0("beta.", traits)

  se <- as.data.frame(matrix(runif(n_snp * k, 0.03, 0.08), n_snp, k), check.names = FALSE)
  names(se) <- paste0("se.", traits)

  list(
    covstruc = list(V = V_LD, S = S_LD, I = I_LD),
    SNPs = cbind(SNPs, beta, se)
  )
}

inputs <- make_inputs()

run_commonfactor <- function(fast) {
  options(GenomicSEM.use_rust = TRUE)
  options(GenomicSEM.fast_commonfactor_fit = fast)
  suppressWarnings(commonfactorGWAS(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    parallel = FALSE,
    GC = "standard"
  ))
}

slow <- NULL
fast <- NULL
invisible(capture.output({
  slow <- run_commonfactor(FALSE)
  fast <- run_commonfactor(TRUE)
}))

stopifnot(identical(names(slow), names(fast)))
stopifnot(max(abs(slow$est - fast$est), na.rm = TRUE) < 1e-5)
stopifnot(max(abs(slow$se_c - fast$se_c), na.rm = TRUE) < 1e-6)
stopifnot(max(abs(slow$Q - fast$Q), na.rm = TRUE) < 1e-2)
stopifnot(max(abs(slow$Q_pval - fast$Q_pval), na.rm = TRUE) < 1e-4)
stopifnot(all(fast$warning == 0))
