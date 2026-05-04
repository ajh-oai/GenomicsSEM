#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
})

parse_args <- function(args) {
  values <- list(
    n_snp = 1000L,
    k = 5L,
    cores = c(1L, 2L, 4L, 8L, 16L),
    repeats = 1L,
    seed = 20260504L
  )

  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      stop("Arguments must have the form key=value.", call. = FALSE)
    }

    key <- parts[[1]]
    value <- parts[[2]]

    if (key == "n_snp") {
      values$n_snp <- as.integer(value)
    } else if (key == "k") {
      values$k <- as.integer(value)
    } else if (key == "cores") {
      values$cores <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
    } else if (key == "repeats") {
      values$repeats <- as.integer(value)
    } else if (key == "seed") {
      values$seed <- as.integer(value)
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }

  values$cores <- sort(unique(values$cores[values$cores > 0L]))
  values
}

make_inputs <- function(n_snp, k, seed) {
  set.seed(seed)
  traits <- paste0("T", seq_len(k))

  s_ld <- matrix(0.35, k, k)
  diag(s_ld) <- 1
  dimnames(s_ld) <- list(traits, traits)

  v_ld <- diag(1e-4, k * (k + 1L) / 2L)

  i_ld <- matrix(0.02, k, k)
  diag(i_ld) <- 1.05
  dimnames(i_ld) <- list(traits, traits)

  snps <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    CHR = rep(1L, n_snp),
    BP = seq_len(n_snp),
    MAF = runif(n_snp, 0.05, 0.5),
    A1 = rep("A", n_snp),
    A2 = rep("G", n_snp),
    check.names = FALSE
  )

  beta <- as.data.frame(
    matrix(rnorm(n_snp * k, sd = 0.02), n_snp, k),
    check.names = FALSE
  )
  names(beta) <- paste0("beta.", traits)

  se <- as.data.frame(
    matrix(runif(n_snp * k, min = 0.03, max = 0.08), n_snp, k),
    check.names = FALSE
  )
  names(se) <- paste0("se.", traits)

  model <- paste(
    paste0("F1 =~ ", paste(traits, collapse = " + ")),
    "F1 ~ SNP",
    sep = "\n"
  )

  list(
    covstruc = list(V = v_ld, S = s_ld, I = i_ld),
    SNPs = cbind(snps, beta, se),
    model = model
  )
}

bind_usergwas <- function(result) {
  do.call(rbind, result)
}

numeric_max_abs_diff <- function(left, right, columns) {
  max(vapply(
    columns,
    function(column) {
      max(abs(as.numeric(left[[column]]) - as.numeric(right[[column]])), na.rm = TRUE)
    },
    numeric(1)
  ))
}

run_commonfactor <- function(inputs, backend, parallel, cores) {
  fun <- if (backend == "lavaan") commonfactorGWAS else commonfactorGWAS_rust

  suppressWarnings(fun(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    estimation = "DWLS",
    parallel = parallel,
    cores = if (parallel) cores else NULL,
    GC = "standard"
  ))
}

run_usergwas <- function(inputs, backend, parallel, cores) {
  fun <- if (backend == "lavaan") userGWAS else userGWAS_rust

  suppressWarnings(fun(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    model = inputs$model,
    estimation = "DWLS",
    parallel = parallel,
    cores = if (parallel) cores else NULL,
    GC = "standard",
    fix_measurement = TRUE,
    Q_SNP = TRUE,
    printwarn = TRUE
  ))
}

benchmark_once <- function(workflow, backend, parallel, cores, inputs) {
  gc()
  elapsed <- system.time({
    result <- if (workflow == "commonfactorGWAS") {
      run_commonfactor(inputs, backend, parallel, cores)
    } else {
      run_usergwas(inputs, backend, parallel, cores)
    }
  })[["elapsed"]]

  data.frame(
    workflow = workflow,
    backend = backend,
    mode = if (parallel) "parallel" else "sequential",
    cores = if (parallel) cores else 1L,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

check_parallel_equivalence <- function(inputs) {
  commonfactor_seq <- run_commonfactor(inputs, "lavaanrust", FALSE, 1L)
  commonfactor_par <- run_commonfactor(inputs, "lavaanrust", TRUE, 2L)
  usergwas_seq <- bind_usergwas(run_usergwas(inputs, "lavaanrust", FALSE, 1L))
  usergwas_par <- bind_usergwas(run_usergwas(inputs, "lavaanrust", TRUE, 2L))

  data.frame(
    workflow = c("commonfactorGWAS", "userGWAS"),
    max_abs_diff = c(
      numeric_max_abs_diff(
        commonfactor_seq,
        commonfactor_par,
        c("est", "se_c", "Q", "Z_Estimate", "Pval_Estimate")
      ),
      numeric_max_abs_diff(
        usergwas_seq,
        usergwas_par,
        c("est", "SE", "chisq", "Q_SNP", "Q_SNP_df", "Q_SNP_pval")
      )
    ),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
inputs <- make_inputs(args$n_snp, args$k, args$seed)

cat("equivalence_check\n")
print(check_parallel_equivalence(inputs), row.names = FALSE)

timings <- do.call(rbind, lapply(seq_len(args$repeats), function(repeat_id) {
  do.call(rbind, lapply(c("commonfactorGWAS", "userGWAS"), function(workflow) {
    do.call(rbind, lapply(c("lavaan", "lavaanrust"), function(backend) {
      rows <- list(
        benchmark_once(workflow, backend, FALSE, 1L, inputs)
      )
      rows <- c(rows, lapply(args$cores, function(core_count) {
        benchmark_once(workflow, backend, TRUE, core_count, inputs)
      }))

      out <- do.call(rbind, rows)
      out$repeat_id <- repeat_id
      out
    }))
  }))
}))

timings$n_snp <- args$n_snp
timings$k <- args$k
timings$seed <- args$seed

cat("timings\n")
print(timings, row.names = FALSE)
