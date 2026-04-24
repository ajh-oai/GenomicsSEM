make_synthetic_genomicsem_inputs <- function(n_snp = 1000L, k = 6L, seed = 1L) {
  set.seed(seed)

  traits <- paste0("T", seq_len(k))

  S_LD <- matrix(0.35, k, k)
  diag(S_LD) <- 1
  dimnames(S_LD) <- list(traits, traits)

  vech_dim <- k * (k + 1) / 2
  V_LD <- diag(1e-4, vech_dim)

  I_LD <- matrix(0.02, k, k)
  diag(I_LD) <- 1.05
  dimnames(I_LD) <- list(traits, traits)

  maf <- runif(n_snp, 0.05, 0.5)
  beta <- matrix(rnorm(n_snp * k, mean = 0, sd = 0.015), n_snp, k)
  se <- matrix(runif(n_snp * k, min = 0.03, max = 0.08), n_snp, k)

  snps <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    CHR = rep(1L, n_snp),
    BP = seq_len(n_snp),
    MAF = maf,
    A1 = rep("A", n_snp),
    A2 = rep("G", n_snp),
    check.names = FALSE
  )

  beta_df <- as.data.frame(beta, check.names = FALSE)
  names(beta_df) <- paste0("beta.", traits)
  se_df <- as.data.frame(se, check.names = FALSE)
  names(se_df) <- paste0("se.", traits)
  snps <- cbind(snps, beta_df, se_df)

  model <- paste(
    paste0("F1 =~ ", paste(traits, collapse = " + ")),
    "F1 ~ SNP",
    paste0(traits, " ~ 0*SNP", collapse = "\n"),
    sep = "\n"
  )

  list(
    covstruc = list(V = V_LD, S = S_LD, I = I_LD),
    SNPs = snps,
    model = model
  )
}

checksum_result <- function(x) {
  values <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  sum(values[is.finite(values)], na.rm = TRUE)
}

count_result_rows <- function(x) {
  if (is.data.frame(x)) {
    return(nrow(x))
  }

  if (is.list(x)) {
    return(sum(vapply(x, count_result_rows, integer(1))))
  }

  length(x)
}

count_result_cols <- function(x) {
  if (is.data.frame(x)) {
    return(ncol(x))
  }

  if (is.list(x)) {
    cols <- vapply(x, count_result_cols, integer(1))
    return(if (length(cols) == 0L) 0L else max(cols))
  }

  1L
}
