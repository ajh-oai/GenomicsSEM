library(GenomicSEM)

make_inputs <- function(n_snp = 5L, k = 3L) {
  set.seed(11)
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

  beta <- matrix(rnorm(n_snp * k, sd = 0.015), n_snp, k)
  beta <- as.data.frame(beta, check.names = FALSE)
  names(beta) <- paste0("beta.", traits)

  se <- matrix(runif(n_snp * k, min = 0.03, max = 0.08), n_snp, k)
  se <- as.data.frame(se, check.names = FALSE)
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
result <- NULL
invisible(capture.output({
  result <- suppressWarnings(userGWAS(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    model = inputs$model,
    parallel = FALSE,
    GC = "standard",
    fix_measurement = TRUE,
    Q_SNP = FALSE,
    printwarn = FALSE
  ))
}))

stopifnot(length(result) == nrow(inputs$SNPs))
stopifnot(all(vapply(result, function(x) !any(c("error", "warning") %in% names(x)), logical(1))))
