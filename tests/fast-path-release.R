library(GenomicSEM)

make_release_inputs <- function(n_snp = 4L, k = 3L) {
  set.seed(41)
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

old_options <- options(
  GenomicSEM.use_rust = TRUE,
  GenomicSEM.fast_commonfactor_fit = TRUE,
  GenomicSEM.fast_usergwas_fit = TRUE,
  GenomicSEM.fast_diagnostics = FALSE,
  GenomicSEM.fast_strict = TRUE
)
on.exit(options(old_options), add = TRUE)

inputs <- make_release_inputs()

commonfactor_out <- NULL
invisible(capture.output({
  commonfactor_out <- suppressWarnings(commonfactorGWAS(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    parallel = TRUE,
    cores = 2,
    GC = "standard"
  ))
}))

stopifnot(identical(attr(commonfactor_out, "GenomicSEM.fast_path"), "rust_commonfactor_batch"))
stopifnot(identical(as.integer(attr(commonfactor_out, "GenomicSEM.fast_threads")), 2L))
stopifnot(all(is.finite(commonfactor_out$est)))
stopifnot(all(is.finite(commonfactor_out$Q)))

usergwas_out <- NULL
invisible(capture.output({
  usergwas_out <- suppressWarnings(userGWAS(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    model = inputs$model,
    estimation = "DWLS",
    parallel = TRUE,
    cores = 2,
    GC = "standard",
    fix_measurement = TRUE,
    Q_SNP = TRUE,
    printwarn = TRUE
  ))
}))

stopifnot(identical(attr(usergwas_out, "GenomicSEM.fast_path"), "rust_usergwas_batch"))
stopifnot(identical(as.integer(attr(usergwas_out, "GenomicSEM.fast_threads")), 2L))
usergwas_df <- do.call(rbind, usergwas_out)
stopifnot(all(usergwas_df$error == 0))
stopifnot(all(usergwas_df$warning == 0))

strict_error <- tryCatch({
  invisible(capture.output({
    suppressWarnings(commonfactorGWAS(
      covstruc = inputs$covstruc,
      SNPs = inputs$SNPs,
      parallel = FALSE,
      smooth_check = TRUE,
      GC = "standard"
    ))
  }))
  NULL
}, error = function(e) conditionMessage(e))

stopifnot(is.character(strict_error))
stopifnot(grepl("smooth_check=TRUE", strict_error, fixed = TRUE))
