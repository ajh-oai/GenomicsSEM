library(GenomicSEM)

make_inputs <- function(n_snp = 8L, k = 4L) {
  set.seed(31)
  traits <- paste0("T", seq_len(k))

  S_LD <- matrix(0.35, k, k)
  diag(S_LD) <- 1
  dimnames(S_LD) <- list(traits, traits)

  V_LD <- diag(1e-4, k * (k + 1) / 2)

  I_LD <- matrix(0.02, k, k)
  diag(I_LD) <- 1.05
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

  beta <- as.data.frame(matrix(rnorm(n_snp * k, sd = 0.015), n_snp, k), check.names = FALSE)
  names(beta) <- paste0("beta.", traits)

  se <- as.data.frame(matrix(runif(n_snp * k, min = 0.03, max = 0.08), n_snp, k), check.names = FALSE)
  names(se) <- paste0("se.", traits)

  model <- paste(
    paste0("F1 =~ ", paste(traits, collapse = " + ")),
    "F1 ~ SNP",
    paste0(traits, " ~ 0*SNP", collapse = "\n"),
    sep = "\n"
  )

  list(
    covstruc = list(V = V_LD, S = S_LD, I = I_LD),
    SNPs = cbind(SNPs, beta, se),
    model = model
  )
}

inputs <- make_inputs()

run_usergwas <- function(fast, q_snp = FALSE) {
  options(GenomicSEM.use_rust = TRUE)
  options(GenomicSEM.fast_usergwas_fit = fast)
  suppressWarnings(userGWAS(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    model = inputs$model,
    estimation = "DWLS",
    parallel = FALSE,
    GC = "standard",
    fix_measurement = TRUE,
    Q_SNP = q_snp,
    printwarn = TRUE
  ))
}

compare_runs <- function(q_snp) {
  slow <- NULL
  fast <- NULL
  invisible(capture.output({
    slow <- run_usergwas(FALSE, q_snp = q_snp)
    fast <- run_usergwas(TRUE, q_snp = q_snp)
  }))

  slow <- do.call(rbind, slow)
  fast <- do.call(rbind, fast)

  stopifnot(identical(names(slow), names(fast)))
  stopifnot(max(abs(as.numeric(slow$est) - as.numeric(fast$est)), na.rm = TRUE) < 1e-5)
  stopifnot(max(abs(as.numeric(slow$SE) - as.numeric(fast$SE)), na.rm = TRUE) < 1e-6)
  stopifnot(max(abs(as.numeric(slow$chisq) - as.numeric(fast$chisq)), na.rm = TRUE) < 1e-5)
  stopifnot(all(fast$error == 0))
  stopifnot(all(fast$warning == 0))

  if (q_snp) {
    q_rows <- !is.na(slow$Q_SNP)
    stopifnot(any(q_rows))
    stopifnot(max(abs(as.numeric(slow$Q_SNP[q_rows]) - as.numeric(fast$Q_SNP[q_rows])), na.rm = TRUE) < 1e-5)
  }
}

compare_runs(FALSE)
compare_runs(TRUE)
