.lavaan_fast_required_columns <- c("lhs", "op", "rhs", "free")

.lavaan_fast_value <- function(values) {
  finite_values <- values[is.finite(values)]
  if (length(finite_values)) {
    return(finite_values[[1L]])
  }

  0
}

.lavaan_fast_row_value <- function(par_table, row_idx) {
  .lavaan_fast_value(c(
    par_table$est[[row_idx]],
    par_table$ustart[[row_idx]],
    par_table$start[[row_idx]]
  ))
}

.lavaan_fast_compile_par_table <- function(par_table, observed_names) {
  if (!is.data.frame(par_table)) {
    stop("lavaan_fast compiler expects a data.frame parameter table.", call. = FALSE)
  }

  missing_columns <- setdiff(.lavaan_fast_required_columns, names(par_table))
  if (length(missing_columns)) {
    stop(
      sprintf(
        "lavaan_fast compiler parameter table is missing columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.character(observed_names) || !length(observed_names) || any(!nzchar(observed_names))) {
    stop("lavaan_fast compiler requires non-empty observed variable names.", call. = FALSE)
  }

  supported_rows <- par_table$op %in% c("=~", "~", "~~")
  if (!all(supported_rows)) {
    stop(
      sprintf(
        "lavaan_fast compiler does not yet support operators: %s",
        paste(sort(unique(par_table$op[!supported_rows])), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  latent_names <- unique(par_table$lhs[par_table$op == "=~"])
  variable_names <- unique(c(observed_names, latent_names, par_table$lhs, par_table$rhs))

  unknown_observed <- setdiff(observed_names, variable_names)
  if (length(unknown_observed)) {
    stop("lavaan_fast compiler observed variables must appear in the parameter table.", call. = FALSE)
  }

  free_ids <- sort(unique(as.integer(par_table$free[par_table$free > 0L])))
  free_index <- match(as.integer(par_table$free), free_ids)
  free_index[par_table$free <= 0L] <- NA_integer_
  default_free_values <- vapply(
    free_ids,
    function(free_id) {
      rows <- which(par_table$free == free_id)
      .lavaan_fast_value(c(
        par_table$est[rows],
        par_table$ustart[rows],
        par_table$start[rows]
      ))
    },
    numeric(1L)
  )

  lhs_index <- match(par_table$lhs, variable_names)
  rhs_index <- match(par_table$rhs, variable_names)
  fixed_values <- vapply(seq_len(nrow(par_table)), function(row_idx) {
    if (!is.na(free_index[[row_idx]])) {
      return(NA_real_)
    }

    .lavaan_fast_row_value(par_table, row_idx)
  }, numeric(1L))

  row_kind <- ifelse(par_table$op %in% c("=~", "~"), "A", "S")
  op_code <- match(par_table$op, c("=~", "~", "~~"))
  edge_key <- paste(row_kind, lhs_index, rhs_index, sep = ":")
  duplicated_edges <- duplicated(edge_key) & row_kind == "A"
  if (any(duplicated_edges)) {
    stop("lavaan_fast compiler does not yet support duplicate directed edges.", call. = FALSE)
  }

  selector <- matrix(
    0,
    nrow = length(observed_names),
    ncol = length(variable_names),
    dimnames = list(observed_names, variable_names)
  )
  selector[cbind(seq_along(observed_names), match(observed_names, variable_names))] <- 1
  free_row_groups <- lapply(seq_along(free_ids), function(free_position) {
    which(free_index == free_position)
  })
  free_row_offsets <- c(0L, cumsum(lengths(free_row_groups)))
  free_row_indices <- as.integer(unlist(free_row_groups, use.names = FALSE))
  free_labels <- vapply(seq_along(free_ids), function(free_position) {
    row_idx <- free_row_groups[[free_position]][[1L]]
    paste0(par_table$lhs[[row_idx]], par_table$op[[row_idx]], par_table$rhs[[row_idx]])
  }, character(1L))
  stat_names <- .stat_names(observed_names)

  structure(
    list(
      observed_names = observed_names,
      latent_names = latent_names,
      variable_names = variable_names,
      par_table = as.data.frame(par_table, stringsAsFactors = FALSE),
      row_kind = row_kind,
      op_code = op_code,
      lhs_index = lhs_index,
      rhs_index = rhs_index,
      free_ids = free_ids,
      free_index = free_index,
      native_free_index = as.integer(ifelse(is.na(free_index), 0L, free_index)),
      fixed_values = fixed_values,
      native_fixed_values = as.double(ifelse(is.na(fixed_values), 0, fixed_values)),
      default_free_values = default_free_values,
      selector = selector,
      observed_index = as.integer(match(observed_names, variable_names)),
      free_row_offsets = as.integer(free_row_offsets),
      free_row_indices = free_row_indices,
      free_labels = free_labels,
      stat_names = stat_names,
      n_variables = length(variable_names),
      n_observed = length(observed_names),
      n_stats = length(stat_names),
      n_free = length(free_ids)
    ),
    class = "lavaan_fast_compiled"
  )
}

.lavaan_fast_ram_matrices <- function(compiled, free_values = compiled$default_free_values) {
  if (!inherits(compiled, "lavaan_fast_compiled")) {
    stop("lavaan_fast RAM builder expects a compiled model.", call. = FALSE)
  }

  if (length(free_values) != length(compiled$free_ids)) {
    stop("lavaan_fast free-value vector has the wrong length.", call. = FALSE)
  }

  n_variables <- length(compiled$variable_names)
  directed <- matrix(
    0,
    nrow = n_variables,
    ncol = n_variables,
    dimnames = list(compiled$variable_names, compiled$variable_names)
  )
  covariance <- directed

  row_values <- compiled$fixed_values
  free_rows <- which(!is.na(compiled$free_index))
  row_values[free_rows] <- free_values[compiled$free_index[free_rows]]

  for (row_idx in seq_along(compiled$row_kind)) {
    lhs <- compiled$lhs_index[[row_idx]]
    rhs <- compiled$rhs_index[[row_idx]]
    value <- row_values[[row_idx]]

    if (compiled$row_kind[[row_idx]] == "A") {
      if (compiled$par_table$op[[row_idx]] == "=~") {
        directed[rhs, lhs] <- value
      } else {
        directed[lhs, rhs] <- value
      }
    } else {
      covariance[lhs, rhs] <- value
      covariance[rhs, lhs] <- value
    }
  }

  list(A = directed, S = covariance, F = compiled$selector)
}

.lavaan_fast_free_labels <- function(compiled) {
  compiled$free_labels
}

.lavaan_fast_derivative_matrices <- function(compiled, free_position) {
  n_variables <- length(compiled$variable_names)
  d_directed <- matrix(0, nrow = n_variables, ncol = n_variables)
  d_covariance <- d_directed

  rows <- which(compiled$free_index == free_position)
  for (row_idx in rows) {
    lhs <- compiled$lhs_index[[row_idx]]
    rhs <- compiled$rhs_index[[row_idx]]
    op <- compiled$par_table$op[[row_idx]]

    if (op == "=~") {
      d_directed[rhs, lhs] <- d_directed[rhs, lhs] + 1
    } else if (op == "~") {
      d_directed[lhs, rhs] <- d_directed[lhs, rhs] + 1
    } else {
      d_covariance[lhs, rhs] <- d_covariance[lhs, rhs] + 1
      if (lhs != rhs) {
        d_covariance[rhs, lhs] <- d_covariance[rhs, lhs] + 1
      }
    }
  }

  list(A = d_directed, S = d_covariance)
}

.lavaan_fast_implied_covariance <- function(compiled, free_values = compiled$default_free_values) {
  matrices <- .lavaan_fast_ram_matrices(compiled, free_values)
  inverse <- solve(diag(nrow(matrices$A)) - matrices$A)
  implied <- matrices$F %*% inverse %*% matrices$S %*% t(inverse) %*% t(matrices$F)
  dimnames(implied) <- list(compiled$observed_names, compiled$observed_names)
  implied
}

.lavaan_fast_implied_jacobian <- function(compiled, free_values = compiled$default_free_values) {
  matrices <- .lavaan_fast_ram_matrices(compiled, free_values)
  inverse <- solve(diag(nrow(matrices$A)) - matrices$A)
  n_observed <- length(compiled$observed_names)
  jacobian <- matrix(
    0,
    nrow = n_observed * (n_observed + 1L) / 2L,
    ncol = length(compiled$free_ids),
    dimnames = list(compiled$stat_names, compiled$free_labels)
  )

  for (free_position in seq_along(compiled$free_ids)) {
    derivative <- .lavaan_fast_derivative_matrices(compiled, free_position)
    d_inverse <- inverse %*% derivative$A %*% inverse
    d_implied <- matrices$F %*% (
      d_inverse %*% matrices$S %*% t(inverse) +
        inverse %*% derivative$S %*% t(inverse) +
        inverse %*% matrices$S %*% t(d_inverse)
    ) %*% t(matrices$F)

    jacobian[, free_position] <- d_implied[lower.tri(d_implied, diag = TRUE)]
  }

  jacobian
}

.lavaan_fast_implied_surfaces_rust_flat <- function(compiled, free_values = compiled$default_free_values) {
  evaluate_ram_surfaces(
    as.integer(compiled$lhs_index),
    as.integer(compiled$rhs_index),
    as.integer(compiled$op_code),
    compiled$native_free_index,
    compiled$native_fixed_values,
    as.double(free_values),
    compiled$observed_index,
    compiled$free_row_offsets,
    compiled$free_row_indices,
    as.integer(compiled$n_variables)
  )
}

.lavaan_fast_implied_surfaces_rust <- function(compiled, free_values = compiled$default_free_values) {
  surfaces <- .lavaan_fast_implied_surfaces_rust_flat(compiled, free_values)

  implied <- matrix(
    surfaces$implied,
    nrow = compiled$n_observed,
    ncol = compiled$n_observed,
    dimnames = list(compiled$observed_names, compiled$observed_names)
  )
  jacobian <- matrix(
    surfaces$delta,
    nrow = compiled$n_stats,
    ncol = compiled$n_free,
    dimnames = list(compiled$stat_names, compiled$free_labels)
  )

  list(implied = implied, delta = jacobian)
}

.lavaan_fast_fit_dwls_rust <- function(compiled, sample_cov, wls_v, max_iter = 400L, tol = 1e-12) {
  fit <- fit_ram_dwls(
    as.double(sample_cov),
    as.double(wls_v),
    as.integer(compiled$lhs_index),
    as.integer(compiled$rhs_index),
    as.integer(compiled$op_code),
    compiled$native_free_index,
    compiled$native_fixed_values,
    as.double(compiled$default_free_values),
    compiled$observed_index,
    compiled$free_row_offsets,
    compiled$free_row_indices,
    as.integer(compiled$n_variables),
    as.integer(max_iter),
    tol
  )

  fit$implied <- matrix(
    fit$implied,
    nrow = length(compiled$observed_names),
    ncol = length(compiled$observed_names),
    dimnames = list(compiled$observed_names, compiled$observed_names)
  )
  fit$delta <- matrix(
    fit$delta,
    nrow = compiled$n_stats,
    ncol = compiled$n_free,
    dimnames = list(compiled$stat_names, compiled$free_labels)
  )

  fit
}
