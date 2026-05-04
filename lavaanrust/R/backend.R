.lav_matrix_vech_rust <- function(matrix) {
  matrix[lower.tri(matrix, diag = TRUE)]
}

.is_one_factor_dwls_model <- function(model) {
  if (!is.character(model) || length(model) != 1L) {
    return(FALSE)
  }

  lines <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]

  if (length(lines) != 2L) {
    return(FALSE)
  }

  loading_line <- lines[grepl("=~", lines, fixed = TRUE)]
  variance_line <- lines[grepl("~~", lines, fixed = TRUE)]

  if (length(loading_line) != 1L || length(variance_line) != 1L) {
    return(FALSE)
  }

  loading_parts <- strsplit(loading_line, "=~", fixed = TRUE)[[1L]]
  if (length(loading_parts) != 2L) {
    return(FALSE)
  }

  latent <- trimws(loading_parts[[1L]])
  variance_no_space <- gsub("[[:space:]]+", "", variance_line)

  identical(variance_no_space, paste0(latent, "~~1*", latent))
}

.one_factor_latent_name <- function(model) {
  loading_line <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  loading_line <- loading_line[grepl("=~", loading_line, fixed = TRUE)][[1L]]
  trimws(strsplit(loading_line, "=~", fixed = TRUE)[[1L]][[1L]])
}

.one_factor_par_table <- function(observed_names, latent_name, fit) {
  k <- length(observed_names)
  ids <- seq_len(2L * k + 1L)
  free <- c(seq_len(k), 0L, k + seq_len(k))
  est <- c(fit$loadings, 1, fit$residuals)
  se <- c(fit$naive_se[seq_len(k)], 0, fit$naive_se[k + seq_len(k)])
  starts <- c(fit$loadings, 1, fit$residuals)

  data.frame(
    id = ids,
    lhs = c(rep(latent_name, k), latent_name, observed_names),
    op = c(rep("=~", k), "~~", rep("~~", k)),
    rhs = c(observed_names, latent_name, observed_names),
    user = c(rep(1L, k + 1L), rep(0L, k)),
    block = rep(1L, 2L * k + 1L),
    group = rep(1L, 2L * k + 1L),
    free = free,
    ustart = c(rep(NA_real_, k), 1, rep(NA_real_, k)),
    exo = rep(0L, 2L * k + 1L),
    label = rep("", 2L * k + 1L),
    plabel = paste0(".p", ids, "."),
    start = starts,
    est = est,
    se = se,
    stringsAsFactors = FALSE
  )
}

.new_one_factor_fit <- function(model, sample.cov, WLS.V, fit) {
  observed_names <- colnames(sample.cov)
  if (is.null(observed_names)) {
    observed_names <- rownames(sample.cov)
  }
  if (is.null(observed_names)) {
    observed_names <- paste0("V", seq_len(ncol(sample.cov)))
  }

  dimnames(sample.cov) <- list(observed_names, observed_names)
  latent_name <- .one_factor_latent_name(model)
  k <- ncol(sample.cov)
  n_stats <- k * (k + 1L) / 2L
  delta <- matrix(fit$delta, nrow = n_stats, ncol = 2L * k)
  implied <- matrix(fit$implied, nrow = k, ncol = k)
  dimnames(implied) <- dimnames(sample.cov)
  stat_names <- character(n_stats)
  idx <- 1L
  for (col in seq_len(k)) {
    for (row in col:k) {
      stat_names[[idx]] <- paste0(observed_names[[col]], "~~", observed_names[[row]])
      idx <- idx + 1L
    }
  }

  dimnames(delta) <- list(
    stat_names,
    c(
      paste0(latent_name, "=~", observed_names),
      paste0(observed_names, "~~", observed_names)
    )
  )
  cor.lv <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(latent_name, latent_name))
  par_table <- .one_factor_par_table(observed_names, latent_name, fit)
  npar <- 2L * k
  df <- n_stats - npar
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = df,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = cor.lv,
    se = list(theta = implied),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = latent_name,
    model_kind = "one_factor_dwls"
  )
}

#' Rust-backed experimental replacement for the supported `lavaan::sem()` slice.
#'
#' Unsupported syntax deliberately falls back to upstream lavaan so callers can
#' adopt the backend incrementally while keeping identical control flow.
#' @export
sem_rust <- function(model, sample.cov, estimator = "ML", WLS.V = NULL, ...) {
  if (
    identical(estimator, "DWLS") &&
      .is_one_factor_dwls_model(model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    fit <- fit_one_factor_dwls(
      as.double(sample.cov),
      as.double(WLS.V),
      as.integer(nrow(sample.cov)),
      200L,
      1e-12
    )
    return(.new_one_factor_fit(model, sample.cov, WLS.V, fit))
  }

  lavaan::sem(
    model = model,
    sample.cov = sample.cov,
    estimator = estimator,
    WLS.V = WLS.V,
    ...
  )
}

#' Rust-backed experimental replacement for `lavaan::lavaan()`.
#'
#' Model reuse is not implemented yet, so this currently delegates to lavaan.
#' @export
lavaan_rust <- function(...) {
  lavaan::lavaan(...)
}

#' Rust-backed experimental replacement for `lavaan::lavInspect()`.
#' @export
lavInspect_rust <- function(object, what, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    return(lavaan::lavInspect(object, what, ...))
  }

  switch(
    what,
    "delta" = object@delta,
    "WLS.V" = object@WLS.V,
    "converged" = object@converged,
    "cor.lv" = object@cor.lv,
    "fit" = object@fit,
    stop("Unsupported lavaan_rust inspection key: ", what, call. = FALSE)
  )
}

#' Rust-backed experimental replacement for `lavaan::inspect()`.
#' @export
inspect_rust <- function(object, what, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    return(lavaan::inspect(object, what, ...))
  }

  if (missing(what)) {
    return(list(object@implied$cov))
  }

  switch(
    what,
    "list" = object@ParTable,
    "se" = object@se,
    stop("Unsupported lavaan_rust inspection key: ", what, call. = FALSE)
  )
}

#' Rust-backed experimental replacement for `lavaan::parTable()`.
#' @export
parTable_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    return(lavaan::parTable(object, ...))
  }

  object@ParTable
}

#' Rust-backed experimental replacement for `stats::fitted()`.
#' @export
fitted_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    return(stats::fitted(object, ...))
  }

  object@implied
}

#' Rust-backed experimental replacement for `stats::resid()`.
#' @export
resid_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    return(stats::resid(object, ...))
  }

  stop("resid_rust() is not implemented for lavaan_rust_fit objects yet", call. = FALSE)
}

#' Rust-backed experimental replacement for lavaan's parameter extractor.
#' @export
lav_model_get_parameters_rust <- function(...) {
  lavaan::lav_model_get_parameters(...)
}

#' Rust-backed experimental replacement for lavaan's complex-step Jacobian.
#' @export
lav_func_jacobian_complex_rust <- function(...) {
  lavaan::lav_func_jacobian_complex(...)
}
