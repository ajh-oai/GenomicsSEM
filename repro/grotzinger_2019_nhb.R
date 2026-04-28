#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicSEM))

parse_args <- function(args) {
  out <- list(
    model_snps = 1000L,
    prep_snps = 50000L,
    cores = c(1L, 4L),
    out_dir = file.path("repro", "results"),
    skip_prep = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else NA_character_

    if (identical(key, "--model-snps")) {
      out$model_snps <- as.integer(value)
      i <- i + 2L
    } else if (identical(key, "--prep-snps")) {
      out$prep_snps <- as.integer(value)
      i <- i + 2L
    } else if (identical(key, "--cores")) {
      out$cores <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
      i <- i + 2L
    } else if (identical(key, "--out-dir")) {
      out$out_dir <- value
      i <- i + 2L
    } else if (identical(key, "--skip-prep")) {
      out$skip_prep <- TRUE
      i <- i + 1L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }

  out$cores <- sort(unique(out$cores[is.finite(out$cores) & out$cores > 0L]))
  if (length(out$cores) == 0L) {
    out$cores <- 1L
  }

  out
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

as_compare_frame <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x) && all(vapply(x, is.data.frame, logical(1)))) {
    out <- do.call(rbind, Map(function(item, name) {
      item$.result_component <- name
      item
    }, x, names(x) %||% as.character(seq_along(x))))
    rownames(out) <- NULL
    return(out)
  }
  as.data.frame(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

compare_outputs <- function(old, new, tolerance = 1e-6) {
  old_df <- as_compare_frame(old)
  new_df <- as_compare_frame(new)
  common <- intersect(names(old_df), names(new_df))
  numeric_cols <- common[
    vapply(old_df[common], is.numeric, logical(1)) &
      vapply(new_df[common], is.numeric, logical(1))
  ]
  character_cols <- common[
    vapply(old_df[common], is.character, logical(1)) &
      vapply(new_df[common], is.character, logical(1))
  ]

  max_abs_diff <- 0
  if (length(numeric_cols) > 0L && nrow(old_df) == nrow(new_df)) {
    diffs <- mapply(function(a, b) {
      suppressWarnings(max(abs(a - b), na.rm = TRUE))
    }, old_df[numeric_cols], new_df[numeric_cols])
    diffs[!is.finite(diffs)] <- 0
    max_abs_diff <- max(diffs, 0)
  } else if (nrow(old_df) != nrow(new_df)) {
    max_abs_diff <- Inf
  }

  character_equal <- TRUE
  if (length(character_cols) > 0L && nrow(old_df) == nrow(new_df)) {
    character_equal <- all(vapply(character_cols, function(col) {
      identical(old_df[[col]], new_df[[col]])
    }, logical(1)))
  } else if (nrow(old_df) != nrow(new_df)) {
    character_equal <- FALSE
  }

  list(
    max_abs_diff = max_abs_diff,
    equivalent = is.finite(max_abs_diff) && max_abs_diff <= tolerance && character_equal,
    numeric_cols = paste(numeric_cols, collapse = ","),
    character_equal = character_equal
  )
}

with_options <- function(opts, code) {
  old <- options(opts)
  on.exit(options(old), add = TRUE)
  force(code)
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

new_options <- function() {
  list(
    GenomicSEM.use_rust = TRUE,
    GenomicSEM.fast_table_read = TRUE,
    GenomicSEM.fast_snp_join = FALSE,
    GenomicSEM.fast_munge_qc = TRUE,
    GenomicSEM.fast_sumstats_qc = TRUE,
    GenomicSEM.fast_munge_engine = TRUE,
    GenomicSEM.fast_sumstats_engine = TRUE,
    GenomicSEM.fast_ldsc_read = TRUE,
    GenomicSEM.fast_ldsc_blocks = TRUE,
    GenomicSEM.fast_commonfactor_fit = TRUE,
    GenomicSEM.fast_usergwas_fit = TRUE,
    GenomicSEM.fast_diagnostics = FALSE,
    GenomicSEM.fast_strict = TRUE
  )
}

make_grotzinger_ldsc <- function() {
  traits <- c("ALCH", "PTSD", "MDD", "ANX")

  S <- matrix(c(
    0.13938790, 0.05977947, 0.05943021, 0.08500236,
    0.05977947, 0.23937808, 0.05799439, 0.11428679,
    0.05943021, 0.05799439, 0.08503281, 0.12667327,
    0.08500236, 0.11428679, 0.12667327, 0.23329360
  ), nrow = 4L, byrow = TRUE, dimnames = list(traits, traits))

  I <- matrix(c(
    1.01069811, 0.023791915, 0.015906770, 0.017004740,
    0.02379192, 0.993973587, 0.004565391, 0.003975452,
    0.01590677, 0.004565391, 0.994969603, 0.262655462,
    0.01700474, 0.003975452, 0.262655462, 1.007313416
  ), nrow = 4L, byrow = TRUE, dimnames = list(traits, traits))

  V <- matrix(c(
    6.041603e-04, 4.630259e-05, 1.329209e-05, 2.110818e-05, 1.699375e-04, 1.928842e-05, 5.939947e-06, 3.078558e-06, 6.866515e-06, 1.292371e-05,
    4.630259e-05, 1.683117e-03, 5.491624e-06, 1.012613e-04, 7.541207e-05, 1.794540e-05, 9.955430e-05, 2.224531e-07, 3.119102e-06, 1.406092e-05,
    1.329209e-05, 5.491624e-06, 3.692467e-05, 4.865869e-05, 8.966010e-05, 1.773068e-06, 2.035833e-05, 6.028298e-06, 9.777927e-06, 7.764292e-06,
    2.110818e-05, 1.012613e-04, 4.865869e-05, 2.627033e-04, 1.244898e-05, -7.693756e-08, 1.510451e-05, 6.091000e-06, 1.858721e-05, 3.524720e-05,
    1.699375e-04, 7.541207e-05, 8.966010e-05, 1.244898e-05, 1.135989e-02, 4.536925e-05, 3.350365e-04, -1.606011e-05, 5.419991e-05, 1.212087e-06,
    1.928842e-05, 1.794540e-05, 1.773068e-06, -7.693756e-08, 4.536925e-05, 1.336556e-04, 1.255636e-04, 7.752931e-06, 1.314640e-05, 3.972090e-05,
    5.939947e-06, 9.955430e-05, 2.035833e-05, 1.510451e-05, 3.350365e-04, 1.255636e-04, 8.250458e-04, -4.762316e-07, 5.242814e-06, 3.447767e-05,
    3.078558e-06, 2.224531e-07, 6.028298e-06, 6.091000e-06, -1.606011e-05, 7.752931e-06, -4.762316e-07, 1.219319e-05, 1.921343e-05, 1.590190e-05,
    6.866515e-06, 3.119102e-06, 9.777927e-06, 1.858721e-05, 5.419991e-05, 1.314640e-05, 5.242814e-06, 1.921343e-05, 5.109221e-05, 5.444954e-05,
    1.292371e-05, 1.406092e-05, 7.764292e-06, 3.524720e-05, 1.212087e-06, 3.972090e-05, 3.447767e-05, 1.590190e-05, 5.444954e-05, 2.843342e-04
  ), nrow = 10L, byrow = TRUE)

  list(
    V = V,
    S = S,
    I = I,
    N = matrix(c(21647.63, 11235.43, 96488.81, 32830.69, 5831.346, 49943.61, 16951.93, 427751, 145314.9, 49279.86), nrow = 1L),
    m = 1173569
  )
}

reproduce_common_factor <- function(covstruc) {
  expected_fit <- c(chisq = 1.283452, df = 2, p_chisq = 0.526383, AIC = 17.28345, CFI = 1, SRMR = 0.03621695)
  expected_loadings <- data.frame(
    rhs = c("MDD", "PTSD", "ALCH", "ANX"),
    Unstand_Est = c(0.283806747, 0.221278068, 0.205225639, 0.445784745),
    STD_Genotype = c(0.97326125, 0.45226835, 0.54969157, 0.92294074),
    stringsAsFactors = FALSE
  )

  timed <- time_expr({
    capture.output({
      result <- suppressWarnings(usermodel(
        covstruc = covstruc,
        model = "F1=~MDD+PTSD+ALCH+ANX",
        std.lv = TRUE
      ))
    })
    result
  })

  fit <- unlist(timed$value$modelfit[1, names(expected_fit)])
  result_loadings <- timed$value$results[timed$value$results$op == "=~", ]
  result_loadings <- result_loadings[match(expected_loadings$rhs, result_loadings$rhs), ]

  data.frame(
    stage = "published_common_factor",
    backend = "GenomicSEM_usermodel",
    cores = 1L,
    elapsed_sec = timed$elapsed,
    n_snp = NA_integer_,
    rows = nrow(timed$value$results),
    cols = ncol(timed$value$results),
    checksum = numeric_checksum(timed$value),
    max_abs_diff_vs_old = NA_real_,
    equivalent_to_old = NA,
    max_abs_diff_vs_published = max(
      abs(as.numeric(fit) - expected_fit),
      abs(result_loadings$Unstand_Est - expected_loadings$Unstand_Est),
      abs(result_loadings$STD_Genotype - expected_loadings$STD_Genotype),
      na.rm = TRUE
    ),
    equivalent_to_published = max(
      abs(as.numeric(fit) - expected_fit),
      abs(result_loadings$Unstand_Est - expected_loadings$Unstand_Est),
      abs(result_loadings$STD_Genotype - expected_loadings$STD_Genotype),
      na.rm = TRUE
    ) < 5e-6,
    fast_path = NA_character_,
    fast_threads = NA_integer_,
    compared_numeric_cols = NA_character_,
    stringsAsFactors = FALSE
  )
}

make_snp_fixture <- function(n_snp, covstruc, seed = 2019L) {
  set.seed(seed)
  traits <- colnames(covstruc$S)
  maf <- runif(n_snp, 0.05, 0.49)
  common_effect <- rnorm(n_snp, sd = 0.002)
  loadings <- c(ALCH = 0.205225639, PTSD = 0.221278068, MDD = 0.283806747, ANX = 0.445784745)

  snps <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    CHR = rep(4L, n_snp),
    BP = 68000L + seq_len(n_snp),
    MAF = maf,
    A1 = rep("A", n_snp),
    A2 = rep("G", n_snp),
    check.names = FALSE
  )
  if (n_snp >= 5L) {
    snps$SNP[1:5] <- c("rs10030871", "rs6599368", "rs7678633", "rs13130581", "rs13125929")
    snps$BP[1:5] <- c(68786L, 69567L, 69713L, 70392L, 71566L)
    snps$MAF[1:5] <- c(0.0765408, 0.0755467, 0.0755467, 0.0725646, 0.0725646)
    snps$A1[1:5] <- c("T", "A", "G", "A", "T")
    snps$A2[1:5] <- c("C", "T", "A", "G", "C")
  }

  beta <- sapply(traits, function(trait) {
    loadings[[trait]] * common_effect + rnorm(n_snp, sd = 0.001)
  })
  se <- sapply(traits, function(trait) {
    runif(n_snp, min = 0.006, max = 0.012)
  })

  beta_df <- as.data.frame(beta, check.names = FALSE)
  names(beta_df) <- paste0("beta.", traits)
  se_df <- as.data.frame(se, check.names = FALSE)
  names(se_df) <- paste0("se.", traits)

  cbind(snps, beta_df, se_df)
}

run_usergwas <- function(covstruc, SNPs, cores, use_new) {
  opts <- if (use_new) new_options() else old_options()
  with_options(opts, {
    timed <- time_expr({
      capture.output({
        result <- suppressWarnings(userGWAS(
          covstruc = covstruc,
          SNPs = SNPs,
          model = "F1=~MDD+PTSD+ALCH+ANX\nF1~SNP",
          sub = "F1~SNP",
          parallel = cores > 1L,
          cores = if (cores > 1L) cores else NULL,
          Q_SNP = TRUE,
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
    fast_path = attr(timed$value, "GenomicSEM.fast_path") %||% NA_character_,
    fast_threads = attr(timed$value, "GenomicSEM.fast_threads") %||% NA_integer_
  )
}

run_commonfactor_gwas <- function(covstruc, SNPs, cores, use_new) {
  opts <- if (use_new) new_options() else old_options()
  with_options(opts, {
    timed <- time_expr({
      capture.output({
        result <- suppressWarnings(commonfactorGWAS(
          covstruc = covstruc,
          SNPs = SNPs,
          parallel = cores > 1L,
          cores = if (cores > 1L) cores else NULL,
          GC = "standard"
        ))
      })
      result
    })
  })

  list(
    result = timed$value,
    elapsed = timed$elapsed,
    fast_path = attr(timed$value, "GenomicSEM.fast_path") %||% NA_character_,
    fast_threads = attr(timed$value, "GenomicSEM.fast_threads") %||% NA_integer_
  )
}

benchmark_model_workflow <- function(stage, run_fun, covstruc, SNPs, cores) {
  old <- run_fun(covstruc, SNPs, cores, FALSE)
  new <- run_fun(covstruc, SNPs, cores, TRUE)
  comparison <- compare_outputs(old$result, new$result, tolerance = 1e-5)

  rbind(
    data.frame(
      stage = stage,
      backend = "old_r_lavaan",
      cores = cores,
      elapsed_sec = old$elapsed,
      n_snp = nrow(SNPs),
      rows = count_rows(old$result),
      cols = count_cols(old$result),
      checksum = numeric_checksum(old$result),
      max_abs_diff_vs_old = 0,
      equivalent_to_old = TRUE,
      max_abs_diff_vs_published = NA_real_,
      equivalent_to_published = NA,
      fast_path = old$fast_path,
      fast_threads = old$fast_threads,
      compared_numeric_cols = comparison$numeric_cols,
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = stage,
      backend = "new_rust_binding",
      cores = cores,
      elapsed_sec = new$elapsed,
      n_snp = nrow(SNPs),
      rows = count_rows(new$result),
      cols = count_cols(new$result),
      checksum = numeric_checksum(new$result),
      max_abs_diff_vs_old = comparison$max_abs_diff,
      equivalent_to_old = comparison$equivalent,
      max_abs_diff_vs_published = NA_real_,
      equivalent_to_published = NA,
      fast_path = new$fast_path,
      fast_threads = new$fast_threads,
      compared_numeric_cols = comparison$numeric_cols,
      stringsAsFactors = FALSE
    )
  )
}

make_prep_fixture <- function(n_snp, dir, seed = 2020L) {
  set.seed(seed)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  alleles <- c("A", "C", "G", "T")
  a1 <- sample(alleles, n_snp, replace = TRUE)
  a2 <- sample(alleles, n_snp, replace = TRUE)
  same <- a1 == a2
  a2[same] <- alleles[(match(a1[same], alleles) %% length(alleles)) + 1L]

  ref <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    A1 = a1,
    A2 = a2,
    MAF = runif(n_snp, 0.02, 0.49),
    stringsAsFactors = FALSE
  )
  write.table(ref, file.path(dir, "reference.1000G.paper_shape.txt"), row.names = FALSE, quote = FALSE)
  write.table(ref[c("SNP", "A1", "A2")], file.path(dir, "w_hm3.paper_shape.snplist"), row.names = FALSE, quote = FALSE)

  files <- file.path(dir, paste0(c("ALCH", "PTSD", "MDD", "ANX"), "_paper_shape.txt"))
  for (i in seq_along(files)) {
    flip <- (seq_len(n_snp) + i) %% 5L == 0L
    trait_a1 <- ref$A1
    trait_a2 <- ref$A2
    trait_a1[flip] <- ref$A2[flip]
    trait_a2[flip] <- ref$A1[flip]

    z <- rnorm(n_snp)
    trait <- data.frame(
      SNP = ref$SNP,
      A1 = trait_a1,
      A2 = trait_a2,
      BETA = z * runif(n_snp, 0.006, 0.012),
      SE = runif(n_snp, 0.006, 0.012),
      P = format.pval(2 * pnorm(-abs(z)), digits = 8),
      N = sample(20000:200000, n_snp, replace = TRUE),
      INFO = runif(n_snp, 0.91, 1.0),
      MAF = ref$MAF,
      stringsAsFactors = FALSE
    )
    write.table(trait, files[[i]], row.names = FALSE, quote = FALSE)
  }

  list(
    files = files,
    trait_names = c("ALCH", "PTSD", "MDD", "ANX"),
    ref = file.path(dir, "reference.1000G.paper_shape.txt"),
    hm3 = file.path(dir, "w_hm3.paper_shape.snplist")
  )
}

benchmark_sumstats <- function(fixture, use_new) {
  opts <- if (use_new) new_options() else old_options()
  with_options(opts, {
    time_expr({
      capture.output({
        result <- suppressWarnings(sumstats(
          files = fixture$files,
          ref = fixture$ref,
          trait.names = fixture$trait_names,
          se.logit = rep(TRUE, length(fixture$files)),
          parallel = FALSE
        ))
      })
      result
    })
  })
}

benchmark_munge <- function(fixture, use_new) {
  opts <- if (use_new) new_options() else old_options()
  trait_names <- paste0(fixture$trait_names, if (use_new) "_new" else "_old")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(dirname(fixture$hm3))

  local_files <- basename(fixture$files)
  with_options(opts, {
    time_expr({
      capture.output({
        suppressWarnings(munge(
          files = local_files,
          hm3 = basename(fixture$hm3),
          trait.names = trait_names,
          parallel = FALSE,
          overwrite = TRUE,
          column.names = list(effect = "BETA")
        ))
      })
      lapply(trait_names, function(name) {
        read.table(gzfile(paste0(name, ".sumstats.gz")), header = TRUE)
      })
    })
  })
}

benchmark_prep_workflow <- function(stage, old_timed, new_timed, n_snp) {
  comparison <- compare_outputs(old_timed$value, new_timed$value, tolerance = 1e-5)
  rbind(
    data.frame(
      stage = stage,
      backend = "old_r_prep",
      cores = 1L,
      elapsed_sec = old_timed$elapsed,
      n_snp = n_snp,
      rows = count_rows(old_timed$value),
      cols = count_cols(old_timed$value),
      checksum = numeric_checksum(old_timed$value),
      max_abs_diff_vs_old = 0,
      equivalent_to_old = TRUE,
      max_abs_diff_vs_published = NA_real_,
      equivalent_to_published = NA,
      fast_path = NA_character_,
      fast_threads = NA_integer_,
      compared_numeric_cols = comparison$numeric_cols,
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = stage,
      backend = "new_rust_binding",
      cores = 1L,
      elapsed_sec = new_timed$elapsed,
      n_snp = n_snp,
      rows = count_rows(new_timed$value),
      cols = count_cols(new_timed$value),
      checksum = numeric_checksum(new_timed$value),
      max_abs_diff_vs_old = comparison$max_abs_diff,
      equivalent_to_old = comparison$equivalent,
      max_abs_diff_vs_published = NA_real_,
      equivalent_to_published = NA,
      fast_path = NA_character_,
      fast_threads = NA_integer_,
      compared_numeric_cols = comparison$numeric_cols,
      stringsAsFactors = FALSE
    )
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

  covstruc <- make_grotzinger_ldsc()
  SNPs <- make_snp_fixture(args$model_snps, covstruc)

  rows <- list(reproduce_common_factor(covstruc))
  for (core in args$cores) {
    rows[[length(rows) + 1L]] <- benchmark_model_workflow(
      "paper_shaped_userGWAS_Q_SNP",
      run_usergwas,
      covstruc,
      SNPs,
      core
    )
    rows[[length(rows) + 1L]] <- benchmark_model_workflow(
      "paper_shaped_commonfactorGWAS",
      run_commonfactor_gwas,
      covstruc,
      SNPs,
      core
    )
  }

  if (!args$skip_prep) {
    prep_dir <- tempfile("grotzinger-prep-")
    fixture <- make_prep_fixture(args$prep_snps, prep_dir)
    rows[[length(rows) + 1L]] <- benchmark_prep_workflow(
      "paper_shaped_sumstats",
      benchmark_sumstats(fixture, FALSE),
      benchmark_sumstats(fixture, TRUE),
      args$prep_snps
    )
    rows[[length(rows) + 1L]] <- benchmark_prep_workflow(
      "paper_shaped_munge",
      benchmark_munge(fixture, FALSE),
      benchmark_munge(fixture, TRUE),
      args$prep_snps
    )
  }

  results <- do.call(rbind, rows)
  csv_path <- file.path(args$out_dir, paste0("grotzinger_2019_nhb_", timestamp, ".csv"))
  write.csv(results, csv_path, row.names = FALSE)

  cat("Wrote ", csv_path, "\n", sep = "")
  print(results, row.names = FALSE)
}

main()
