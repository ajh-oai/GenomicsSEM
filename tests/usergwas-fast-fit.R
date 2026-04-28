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
constrained_inputs <- make_inputs()
constraint_traits <- colnames(constrained_inputs$covstruc$S)
constraint_labels <- paste0("theta", seq_along(constraint_traits))
constrained_inputs$model <- paste(
  constrained_inputs$model,
  paste0(constraint_traits, " ~~ ", constraint_labels, "*", constraint_traits, collapse = "\n"),
  paste0(constraint_labels, " > 0.001", collapse = "\n"),
  sep = "\n"
)

bound_spec <- getFromNamespace(".sem_fast_compile", "GenomicSEM")(
  lavaan::lavaanify(
    "F1 =~ T1 + T2 + T3\nT1 ~~ theta*T1\ntheta > 0.001",
    auto.var = TRUE,
    auto.fix.first = TRUE
  ),
  paste0("T", 1:3)
)
stopifnot(isTRUE(bound_spec$supported))
stopifnot(any(abs(bound_spec$lower - 0.001) < .Machine$double.eps^0.5))

run_usergwas <- function(fast, q_snp = FALSE, input_data = inputs) {
  options(GenomicSEM.use_rust = TRUE)
  options(GenomicSEM.fast_usergwas_fit = fast)
  suppressWarnings(userGWAS(
    covstruc = input_data$covstruc,
    SNPs = input_data$SNPs,
    model = input_data$model,
    estimation = "DWLS",
    parallel = FALSE,
    GC = "standard",
    fix_measurement = TRUE,
    Q_SNP = q_snp,
    printwarn = TRUE
  ))
}

run_usergwas_sub <- function(fast, input_data = inputs) {
  options(GenomicSEM.use_rust = TRUE)
  options(GenomicSEM.fast_usergwas_fit = fast)
  suppressWarnings(userGWAS(
    covstruc = input_data$covstruc,
    SNPs = input_data$SNPs,
    model = input_data$model,
    sub = "F1~SNP",
    estimation = "DWLS",
    parallel = FALSE,
    GC = "standard",
    fix_measurement = TRUE,
    Q_SNP = TRUE,
    printwarn = TRUE
  ))
}

compare_runs <- function(q_snp, input_data = inputs) {
  slow <- NULL
  fast <- NULL
  invisible(capture.output({
    slow <- run_usergwas(FALSE, q_snp = q_snp, input_data = input_data)
    fast <- run_usergwas(TRUE, q_snp = q_snp, input_data = input_data)
  }))

  stopifnot(identical(attr(fast, "GenomicSEM.fast_path"), "rust_usergwas_batch"))
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

compare_sub_runs <- function() {
  slow <- NULL
  fast <- NULL
  invisible(capture.output({
    slow <- run_usergwas_sub(FALSE)
    fast <- run_usergwas_sub(TRUE)
  }))

  stopifnot(length(slow) == 1L)
  stopifnot(length(fast) == 1L)
  stopifnot(identical(names(slow[[1]]), names(fast[[1]])))
  stopifnot(identical(slow[[1]]$SNP, fast[[1]]$SNP))
  stopifnot(max(abs(as.numeric(slow[[1]]$est) - as.numeric(fast[[1]]$est)), na.rm = TRUE) < 1e-5)
  stopifnot(max(abs(as.numeric(slow[[1]]$SE) - as.numeric(fast[[1]]$SE)), na.rm = TRUE) < 1e-6)
  stopifnot(max(abs(as.numeric(slow[[1]]$Q_SNP) - as.numeric(fast[[1]]$Q_SNP)), na.rm = TRUE) < 1e-5)
  stopifnot(all(fast[[1]]$error == 0))
  stopifnot(all(fast[[1]]$warning == 0))
}

compare_runs(FALSE)
compare_runs(TRUE)
compare_runs(FALSE, constrained_inputs)
compare_runs(TRUE, constrained_inputs)
compare_sub_runs()
