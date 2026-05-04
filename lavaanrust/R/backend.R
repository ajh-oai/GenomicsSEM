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

.stat_names <- function(observed_names) {
  k <- length(observed_names)
  n_stats <- k * (k + 1L) / 2L
  stat_names <- character(n_stats)
  idx <- 1L

  for (col in seq_len(k)) {
    for (row in col:k) {
      stat_names[[idx]] <- paste0(observed_names[[col]], "~~", observed_names[[row]])
      idx <- idx + 1L
    }
  }

  stat_names
}

.vech_rust <- function(matrix) {
  matrix[lower.tri(matrix, diag = TRUE)]
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
  dimnames(delta) <- list(
    .stat_names(observed_names),
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
    observed = sample.cov,
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

.is_commonfactor_null_model <- function(model) {
  is.character(model) &&
    length(model) == 1L &&
    grepl("VF1 =~ 1*", model, fixed = TRUE) &&
    grepl("~~ 0*", model, fixed = TRUE) &&
    !grepl(":=", model, fixed = TRUE)
}

.commonfactor_null_par_table <- function(observed_names) {
  k <- length(observed_names)
  latent_names <- paste0("VF", seq_len(k))
  rows <- vector("list", 0L)

  add_row <- function(lhs, op, rhs, free, ustart, est) {
    rows[[length(rows) + 1L]] <<- list(
      lhs = lhs,
      op = op,
      rhs = rhs,
      free = free,
      ustart = ustart,
      est = est
    )
  }

  for (idx in seq_len(k)) {
    add_row(observed_names[[idx]], "~~", observed_names[[idx]], idx, NA_real_, NA_real_)
  }

  for (idx in seq_len(k)) {
    add_row(latent_names[[idx]], "=~", observed_names[[idx]], 0L, 1, 1)
  }

  if (k >= 2L) {
    for (lhs in seq_len(k - 1L)) {
      for (rhs in seq.int(lhs + 1L, k)) {
        add_row(latent_names[[lhs]], "~~", latent_names[[rhs]], 0L, 0, 0)
      }
    }

    for (lhs in seq_len(k - 1L)) {
      for (rhs in seq.int(lhs + 1L, k)) {
        add_row(observed_names[[lhs]], "~~", observed_names[[rhs]], 0L, 0, 0)
      }
    }
  }

  for (idx in seq_len(k)) {
    add_row(latent_names[[idx]], "~~", latent_names[[idx]], 0L, 0, 0)
  }

  raw <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  n <- nrow(raw)
  data.frame(
    id = seq_len(n),
    lhs = raw$lhs,
    op = raw$op,
    rhs = raw$rhs,
    user = rep(1L, n),
    block = rep(1L, n),
    group = rep(1L, n),
    free = as.integer(raw$free),
    ustart = as.numeric(raw$ustart),
    exo = rep(0L, n),
    label = rep("", n),
    plabel = paste0(".p", seq_len(n), "."),
    start = as.numeric(raw$est),
    est = as.numeric(raw$est),
    se = rep(0, n),
    stringsAsFactors = FALSE
  )
}

.observed_covariance_spec <- function(par_table, sample.cov) {
  observed_names <- colnames(sample.cov)
  if (is.null(observed_names)) {
    observed_names <- rownames(sample.cov)
  }
  if (is.null(observed_names)) {
    observed_names <- paste0("V", seq_len(ncol(sample.cov)))
  }

  dimnames(sample.cov) <- list(observed_names, observed_names)
  stat_names <- .stat_names(observed_names)
  free_mask <- integer(length(stat_names))
  fixed_values <- numeric(length(stat_names))
  row_lookup <- integer(length(stat_names))
  loading_rows <- par_table[par_table$op == "=~", , drop = FALSE]
  latent_for_observed <- stats::setNames(loading_rows$lhs, loading_rows$rhs)

  par_value <- function(row_idx) {
    if (!length(row_idx)) {
      return(0)
    }

    if (!is.na(par_table$ustart[[row_idx]])) {
      return(par_table$ustart[[row_idx]])
    }

    if (!is.na(par_table$est[[row_idx]])) {
      return(par_table$est[[row_idx]])
    }

    0
  }

  for (idx in seq_along(stat_names)) {
    parts <- strsplit(stat_names[[idx]], "~~", fixed = TRUE)[[1L]]
    direct_row <- which(
      par_table$op == "~~" &
        par_table$lhs == parts[[1L]] &
        par_table$rhs == parts[[2L]]
    )

    if (!length(direct_row)) {
      direct_row <- which(
        par_table$op == "~~" &
          par_table$lhs == parts[[2L]] &
          par_table$rhs == parts[[1L]]
      )
    }

    if (!length(direct_row)) {
      stop("Observed covariance row missing from parameter table: ", stat_names[[idx]], call. = FALSE)
    }

    direct_row <- direct_row[[1L]]

    if (identical(parts[[1L]], parts[[2L]]) && parts[[1L]] %in% names(latent_for_observed)) {
      latent <- latent_for_observed[[parts[[1L]]]]
      latent_var_row <- which(
        par_table$op == "~~" &
          par_table$lhs == latent &
          par_table$rhs == latent
      )
      latent_var_row <- if (length(latent_var_row)) latent_var_row[[1L]] else integer()

      direct_free <- par_table$free[[direct_row]] > 0L
      latent_free <- length(latent_var_row) && par_table$free[[latent_var_row]] > 0L

      if (direct_free && latent_free) {
        stop("Unsupported model: both observed and latent variance rows are free", call. = FALSE)
      }

      if (direct_free) {
        free_mask[[idx]] <- 1L
        row_lookup[[idx]] <- direct_row
        fixed_values[[idx]] <- par_value(latent_var_row)
      } else if (latent_free) {
        free_mask[[idx]] <- 1L
        row_lookup[[idx]] <- latent_var_row
        fixed_values[[idx]] <- par_value(direct_row)
      } else {
        fixed_values[[idx]] <- par_value(direct_row) + par_value(latent_var_row)
      }
    } else {
      row_lookup[[idx]] <- direct_row
      free_mask[[idx]] <- as.integer(par_table$free[[direct_row]] > 0L)

      if (!free_mask[[idx]]) {
        fixed_values[[idx]] <- par_value(direct_row)
      }
    }
  }

  list(
    observed_names = observed_names,
    sample_cov = sample.cov,
    free_mask = free_mask,
    fixed_values = fixed_values,
    row_lookup = row_lookup,
    stat_names = stat_names
  )
}

.new_observed_covariance_fit <- function(par_table, sample.cov, WLS.V, fit, model_kind) {
  spec <- .observed_covariance_spec(par_table, sample.cov)
  implied <- matrix(fit$implied, nrow = nrow(sample.cov), ncol = ncol(sample.cov))
  dimnames(implied) <- list(spec$observed_names, spec$observed_names)
  delta <- matrix(fit$delta, nrow = length(spec$stat_names))
  dimnames(delta) <- list(spec$stat_names, spec$stat_names[spec$free_mask > 0L])
  par_table <- as.data.frame(par_table, stringsAsFactors = FALSE)
  free_rows <- spec$row_lookup[spec$free_mask > 0L]
  par_table$est[free_rows] <- fit$estimates
  par_table$se[] <- 0
  par_table$se[free_rows] <- fit$naive_se
  npar <- sum(spec$free_mask)
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = length(spec$stat_names) - npar,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = spec$sample_cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(numeric(), nrow = 0L, ncol = 0L),
    se = list(theta = matrix(0, nrow = nrow(sample.cov), ncol = ncol(sample.cov))),
    converged = isTRUE(fit$converged),
    observed_names = spec$observed_names,
    latent_name = character(),
    model_kind = model_kind
  )
}

.fit_observed_covariance_model <- function(par_table, sample.cov, WLS.V, model_kind) {
  spec <- .observed_covariance_spec(par_table, sample.cov)
  fit <- fit_observed_covariance_dwls(
    as.double(spec$sample_cov),
    as.double(WLS.V),
    as.integer(spec$free_mask),
    as.double(spec$fixed_values),
    as.integer(nrow(spec$sample_cov))
  )

  .new_observed_covariance_fit(par_table, spec$sample_cov, WLS.V, fit, model_kind)
}

#' Rust-backed experimental replacement for the supported `lavaan::sem()` slice.
#'
#' @param model A supported lavaan-style model string or parameter table.
#' @param sample.cov Observed covariance matrix.
#' @param estimator Estimator name. Only `"DWLS"` is currently supported.
#' @param WLS.V DWLS weight matrix.
#' @param ... Additional lavaan-style arguments. Present for signature
#'   compatibility; unsupported paths still error.
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

  if (
    identical(estimator, "DWLS") &&
      .is_commonfactor_null_model(model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    observed_names <- colnames(sample.cov)
    if (is.null(observed_names)) {
      observed_names <- rownames(sample.cov)
    }
    if (is.null(observed_names)) {
      observed_names <- paste0("V", seq_len(ncol(sample.cov)))
    }

    return(.fit_observed_covariance_model(
      .commonfactor_null_par_table(observed_names),
      sample.cov,
      WLS.V,
      "commonfactor_null_dwls"
    ))
  }

  if (
    identical(estimator, "DWLS") &&
      is.data.frame(model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    return(.fit_observed_covariance_model(
      model,
      sample.cov,
      WLS.V,
      "observed_covariance_par_table_dwls"
    ))
  }

  stop(
    "Unsupported sem_rust() model path. The rust wrapper does not fall back to lavaan.",
    call. = FALSE
  )
}

#' Rust-backed experimental replacement for `lavaan::lavaan()`.
#'
#' @param ... Reserved for future lavaan-compatible model-reuse arguments.
#' @export
lavaan_rust <- function(...) {
  stop(
    "lavaan_rust() model-reuse refits are not implemented yet.",
    call. = FALSE
  )
}

#' Rust-backed experimental replacement for `lavaan::lavInspect()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param what Inspection key.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
lavInspect_rust <- function(object, what, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("lavInspect_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
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
#'
#' @param object A `lavaan_rust_fit` object.
#' @param what Inspection key.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
inspect_rust <- function(object, what, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("inspect_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
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
#'
#' @param object A `lavaan_rust_fit` object.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
parTable_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("parTable_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  object@ParTable
}

#' Rust-backed experimental replacement for `stats::fitted()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
fitted_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("fitted_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  object@implied
}

#' Rust-backed experimental replacement for `stats::resid()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
resid_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("resid_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  list(cov = object@observed - object@implied$cov)
}

#' Rust-backed experimental replacement for lavaan's parameter extractor.
#'
#' @param ... Reserved for a future compatible implementation.
#' @export
lav_model_get_parameters_rust <- function(...) {
  stop("lav_model_get_parameters_rust() is not implemented yet.", call. = FALSE)
}

#' Rust-backed experimental replacement for lavaan's complex-step Jacobian.
#'
#' @param ... Reserved for a future compatible implementation.
#' @export
lav_func_jacobian_complex_rust <- function(...) {
  stop("lav_func_jacobian_complex_rust() is not implemented yet.", call. = FALSE)
}
