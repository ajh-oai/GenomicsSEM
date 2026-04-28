#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicSEM))

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

parse_int_list <- function(x) {
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_args <- function(args) {
  out <- list(
    data_dir = file.path("repro", "data", "pfactor_practical"),
    sumstats = NULL,
    out_dir = file.path("repro", "results"),
    sizes = c(1000L, 5000L, 10000L, 25000L, 50000L),
    cores = 16L,
    q_snp = FALSE,
    new_full_sec = NA_real_
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else NA_character_

    if (identical(key, "--data-dir")) {
      out$data_dir <- value
      i <- i + 2L
    } else if (identical(key, "--sumstats")) {
      out$sumstats <- value
      i <- i + 2L
    } else if (identical(key, "--out-dir")) {
      out$out_dir <- value
      i <- i + 2L
    } else if (identical(key, "--sizes")) {
      out$sizes <- parse_int_list(value)
      i <- i + 2L
    } else if (identical(key, "--cores")) {
      out$cores <- as.integer(value)
      i <- i + 2L
    } else if (identical(key, "--q-snp")) {
      out$q_snp <- TRUE
      i <- i + 1L
    } else if (identical(key, "--new-full-sec")) {
      out$new_full_sec <- as.numeric(value)
      i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }

  out$data_dir <- normalizePath(out$data_dir, mustWork = FALSE)
  out$out_dir <- normalizePath(out$out_dir, mustWork = FALSE)
  if (is.null(out$sumstats)) {
    out$sumstats <- file.path(out$data_dir, "sumstats_new_1000000.rds")
  }
  out$sumstats <- normalizePath(out$sumstats, mustWork = FALSE)
  out$sizes <- sort(unique(out$sizes[is.finite(out$sizes) & out$sizes > 0L]))
  out$cores <- max(1L, out$cores)

  if (length(out$sizes) == 0L) {
    stop("--sizes must contain at least one positive integer", call. = FALSE)
  }

  out
}

old_options <- function() {
  list(
    GenomicSEM.use_rust = FALSE,
    GenomicSEM.fast_table_read = FALSE,
    GenomicSEM.fast_snp_join = FALSE,
    GenomicSEM.fast_munge_qc = FALSE,
    GenomicSEM.fast_sumstats_qc = FALSE,
    GenomicSEM.fast_munge_engine = FALSE,
    GenomicSEM.fast_sumstats_engine = FALSE,
    GenomicSEM.fast_ldsc_read = FALSE,
    GenomicSEM.fast_ldsc_blocks = FALSE,
    GenomicSEM.fast_commonfactor_fit = FALSE,
    GenomicSEM.fast_usergwas_fit = FALSE,
    GenomicSEM.fast_diagnostics = FALSE,
    GenomicSEM.fast_strict = FALSE
  )
}

with_options <- function(opts, code) {
  old <- options(opts)
  on.exit(options(old), add = TRUE)
  force(code)
}

time_expr <- function(expr) {
  gc()
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = proc.time()[["elapsed"]] - start)
}

numeric_checksum <- function(x) {
  values <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  sum(values[is.finite(values)], na.rm = TRUE)
}

count_rows <- function(x) {
  if (is.data.frame(x)) {
    return(nrow(x))
  }
  if (is.list(x)) {
    return(sum(vapply(x, count_rows, integer(1))))
  }
  length(x)
}

count_cols <- function(x) {
  if (is.data.frame(x)) {
    return(ncol(x))
  }
  if (is.list(x)) {
    cols <- vapply(x, count_cols, integer(1))
    return(if (length(cols) == 0L) 0L else max(cols))
  }
  1L
}

load_practical_covstruc <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env$PSYCH_COV
}

pfactor_model <- function() {
  paste(
    "F1=~SCZ+BIP+MDD",
    "F1~SNP",
    "SCZ~~a*SCZ",
    "BIP~~b*BIP",
    "MDD~~c*MDD",
    "a > .001",
    "b > .001",
    "c > .001",
    sep = "\n"
  )
}

run_old_lavaan <- function(covstruc, snps, cores, q_snp) {
  with_options(old_options(), {
    timed <- time_expr({
      capture.output({
        result <- suppressWarnings(userGWAS(
          covstruc = covstruc,
          SNPs = snps,
          model = pfactor_model(),
          estimation = "DWLS",
          sub = "F1~SNP",
          parallel = cores > 1L,
          cores = if (cores > 1L) cores else NULL,
          Q_SNP = q_snp,
          fix_measurement = TRUE,
          GC = "standard",
          printwarn = TRUE
        ))
      })
      result
    })
  })

  list(
    result = timed$value,
    elapsed = timed$elapsed,
    fast_path = attr(timed$value, "GenomicSEM.fast_path") %||% NA_character_
  )
}

fit_projection <- function(results, new_full_sec) {
  fit <- lm(elapsed_sec ~ n_snp, data = results)
  coef <- coefficients(fit)
  projected <- as.numeric(coef[["(Intercept)"]] + coef[["n_snp"]] * 1000000)

  data.frame(
    model = "elapsed_sec ~ intercept + slope*n_snp",
    intercept_sec = unname(coef[["(Intercept)"]]),
    slope_sec_per_snp = unname(coef[["n_snp"]]),
    projected_1m_sec = projected,
    projected_1m_hours = projected / 3600,
    r_squared = summary(fit)$r.squared,
    new_full_1m_sec = new_full_sec,
    projected_speedup_vs_new = if (is.finite(new_full_sec)) projected / new_full_sec else NA_real_,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

  cov_path <- file.path(args$data_dir, "GenomicSEMPractical.RData")
  if (!file.exists(cov_path)) {
    stop("Missing covariance object: ", cov_path, call. = FALSE)
  }
  if (!file.exists(args$sumstats)) {
    stop("Missing sumstats RDS: ", args$sumstats, call. = FALSE)
  }

  covstruc <- load_practical_covstruc(cov_path)
  all_snps <- readRDS(args$sumstats)
  if (is.list(all_snps) && is.data.frame(all_snps$value)) {
    all_snps <- all_snps$value
  }
  if (!is.data.frame(all_snps)) {
    stop("Expected a data.frame or a list with data.frame element `value` in ", args$sumstats, call. = FALSE)
  }
  max_size <- max(args$sizes)
  if (nrow(all_snps) < max_size) {
    stop("sumstats has only ", nrow(all_snps), " rows; largest requested size is ", max_size, call. = FALSE)
  }

  results <- vector("list", length(args$sizes))
  for (idx in seq_along(args$sizes)) {
    n <- args$sizes[[idx]]
    message(sprintf("[%s] old lavaan userGWAS: n_snp=%d cores=%d", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n, args$cores))
    run <- run_old_lavaan(covstruc, all_snps[seq_len(n), , drop = FALSE], args$cores, args$q_snp)
    result <- run$result

    results[[idx]] <- data.frame(
      stage = "public_pfactor_old_lavaan_scaling",
      backend = "old_r_lavaan",
      cores = args$cores,
      n_snp = n,
      elapsed_sec = run$elapsed,
      per_snp_ms = 1000 * run$elapsed / n,
      rows = count_rows(result),
      cols = count_cols(result),
      checksum = numeric_checksum(result),
      fast_path = run$fast_path,
      q_snp = args$q_snp,
      stringsAsFactors = FALSE
    )
    print(results[[idx]], row.names = FALSE)
    rm(result, run)
    gc()
  }

  results <- do.call(rbind, results)
  projection <- fit_projection(results, args$new_full_sec)

  result_path <- file.path(args$out_dir, paste0("pfactor_old_lavaan_scaling_", timestamp, ".csv"))
  projection_path <- file.path(args$out_dir, paste0("pfactor_old_lavaan_scaling_projection_", timestamp, ".csv"))
  write.csv(results, result_path, row.names = FALSE)
  write.csv(projection, projection_path, row.names = FALSE)

  cat("Wrote ", result_path, "\n", sep = "")
  cat("Wrote ", projection_path, "\n", sep = "")
  print(results, row.names = FALSE)
  print(projection, row.names = FALSE)
}

main()
