#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicSEM))

parse_args <- function(args) {
  out <- list(
    model_snps = 100L,
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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
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

compare_outputs <- function(old, new, tolerance = 1e-5) {
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

make_grotzinger_2025_ldsc <- function() {
  traits <- c("AN", "OCD", "TS", "SCZ", "BIP", "ASD", "ADHD", "PTSD", "MD", "ANX", "CUD", "AUD", "NIC", "OUD")

  S_reported <- matrix(c(
    0.155, 0.09, 0.018, 0.043, 0.033, 0.014, 0.013, 0.028, 0.026, 0.042, 0.003, 0.024, 0.012, 0.032,
    0.09, 0.161, 0.077, 0.065, 0.057, 0.045, 0.038, 0.042, 0.047, 0.086, 0.014, 0.029, 0.01, 0.021,
    0.018, 0.077, 0.224, 0.02, 0.018, 0.021, 0.045, 0.024, 0.025, 0.044, 0.011, 0.006, -0.006, -0.004,
    0.043, 0.065, 0.02, 0.223, 0.138, 0.041, 0.039, 0.042, 0.042, 0.064, 0.071, 0.062, 0.027, 0.04,
    0.033, 0.057, 0.018, 0.138, 0.187, 0.031, 0.041, 0.047, 0.051, 0.063, 0.042, 0.049, 0.033, 0.028,
    0.014, 0.045, 0.021, 0.041, 0.031, 0.119, 0.062, 0.025, 0.028, 0.033, 0.009, 0.0005, 0.025, -0.006,
    0.013, 0.038, 0.045, 0.039, 0.041, 0.062, 0.18, 0.073, 0.063, 0.071, 0.107, 0.059, 0.061, 0.045,
    0.028, 0.042, 0.024, 0.042, 0.047, 0.025, 0.073, 0.053, 0.057, 0.066, 0.054, 0.035, 0.032, 0.029,
    0.026, 0.047, 0.025, 0.042, 0.051, 0.028, 0.063, 0.057, 0.062, 0.073, 0.045, 0.038, 0.033, 0.031,
    0.042, 0.086, 0.044, 0.064, 0.063, 0.033, 0.071, 0.066, 0.073, 0.108, 0.049, 0.044, 0.036, 0.041,
    0.003, 0.014, 0.011, 0.071, 0.042, 0.009, 0.107, 0.054, 0.045, 0.049, 0.19, 0.092, 0.056, 0.085,
    0.024, 0.029, 0.006, 0.062, 0.049, 0.0005, 0.059, 0.035, 0.038, 0.044, 0.092, 0.115, 0.041, 0.069,
    0.012, 0.01, -0.006, 0.027, 0.033, 0.025, 0.061, 0.032, 0.033, 0.036, 0.056, 0.041, 0.085, 0.038,
    0.032, 0.021, -0.004, 0.04, 0.028, -0.006, 0.045, 0.029, 0.031, 0.041, 0.085, 0.069, 0.038, 0.064
  ), nrow = 14L, byrow = TRUE, dimnames = list(traits, traits))

  ldsc_se <- matrix(c(
    0.01, 0.008, 0.012, 0.006, 0.006, 0.008, 0.006, 0.004, 0.003, 0.004, 0.009, 0.006, 0.007, 0.006,
    0.008, 0.011, 0.012, 0.007, 0.006, 0.007, 0.006, 0.003, 0.003, 0.004, 0.01, 0.006, 0.007, 0.008,
    0.012, 0.012, 0.026, 0.009, 0.009, 0.01, 0.009, 0.005, 0.004, 0.006, 0.014, 0.008, 0.012, 0.012,
    0.006, 0.007, 0.009, 0.008, 0.006, 0.006, 0.005, 0.003, 0.003, 0.004, 0.009, 0.005, 0.006, 0.006,
    0.006, 0.006, 0.009, 0.006, 0.008, 0.006, 0.006, 0.003, 0.003, 0.004, 0.008, 0.005, 0.006, 0.005,
    0.008, 0.007, 0.01, 0.006, 0.006, 0.01, 0.007, 0.003, 0.003, 0.004, 0.009, 0.006, 0.008, 0.007,
    0.006, 0.006, 0.009, 0.005, 0.006, 0.007, 0.008, 0.004, 0.003, 0.004, 0.008, 0.005, 0.007, 0.006,
    0.004, 0.003, 0.005, 0.003, 0.003, 0.003, 0.004, 0.002, 0.002, 0.003, 0.004, 0.003, 0.004, 0.003,
    0.003, 0.003, 0.004, 0.003, 0.003, 0.003, 0.003, 0.002, 0.002, 0.003, 0.004, 0.002, 0.003, 0.003,
    0.004, 0.004, 0.006, 0.004, 0.004, 0.004, 0.004, 0.003, 0.003, 0.004, 0.006, 0.003, 0.004, 0.004,
    0.009, 0.01, 0.014, 0.009, 0.008, 0.009, 0.008, 0.004, 0.004, 0.006, 0.018, 0.008, 0.009, 0.009,
    0.006, 0.006, 0.008, 0.005, 0.005, 0.006, 0.005, 0.003, 0.002, 0.003, 0.008, 0.007, 0.006, 0.005,
    0.007, 0.007, 0.012, 0.006, 0.006, 0.008, 0.007, 0.004, 0.003, 0.004, 0.009, 0.006, 0.012, 0.008,
    0.006, 0.008, 0.012, 0.006, 0.005, 0.007, 0.006, 0.003, 0.003, 0.004, 0.009, 0.005, 0.008, 0.01
  ), nrow = 14L, byrow = TRUE, dimnames = list(traits, traits))

  I <- matrix(c(
    1.02004753951182, 0.101580091111402, 0.00867879117916542, 0.0604592419998623, 0.0717558038104117, 0.12590151118142, 0.0738465627381517, 0.0316543382228407, 0.0902575961032426, 0.0587501481129721, 0.0128020244583516, 0.00365141665784929, -0.0155671845468513, -0.00791605733351862,
    0.101580091111402, 0.997150773604804, 0.0287880491235204, 0.0297181578277577, 0.0605160980441217, 0.0765826072166044, 0.059755087862351, 0.038612098876486, 0.0613311833066134, 0.0636889418869513, 0.00524046477601784, -0.00355524451103263, -0.00729012474606044, 0.0104627183266632,
    0.00867879117916542, 0.0287880491235204, 1.01381954926641, 0.0164350301149663, 0.0371701620588568, 0.00566112835885535, 0.0100516142261298, 0.000219926338451098, 0.0192528968264795, 0.00147409429930324, 0.0101454317512308, 0.0000358383973074061, 0.0104003413752521, 0.00494443990553433,
    0.0604592419998623, 0.0297181578277577, 0.0164350301149663, 1.07802450842275, 0.211707622884005, 0.0240593417446992, 0.036159781891134, 0.0221015581293727, 0.0655378877483725, 0.0351622988914932, 0.0145095634258482, 0.019673776614872, 0.00785055379131484, 0.00996193534311699,
    0.0717558038104117, 0.0605160980441217, 0.0371701620588568, 0.211707622884005, 1.02964300094291, 0.0390901507397886, 0.0629874390114312, 0.025245279746279, 0.0783480862138948, 0.0463883820978562, 0.030760164133941, 0.0216100405840245, 0.00590293806296044, 0.00603421242067271,
    0.12590151118142, 0.0765826072166044, 0.00566112835885535, 0.0240593417446992, 0.0390901507397886, 1.00801961848236, 0.267388886645242, 0.0762162564947683, 0.0782940404757744, 0.049570609596332, 0.0163150127976701, 0.00826771115651044, -0.00803467671394816, 0.00567286440016974,
    0.0738465627381517, 0.059755087862351, 0.0100516142261298, 0.036159781891134, 0.0629874390114312, 0.267388886645242, 1.02529187827736, 0.104033584429931, 0.123856592477088, 0.037198327424167, 0.0705617007565768, 0.00504457407814051, 0.019661904807875, 0.00527052695818655,
    0.0316543382228407, 0.038612098876486, 0.000219926338451098, 0.0221015581293727, 0.025245279746279, 0.0762162564947683, 0.104033584429931, 1.01814564355308, 0.200715672196561, 0.158354051456104, 0.0128999304792121, 0.062949163900661, -0.00201829126129822, 0.0298435457851606,
    0.0902575961032426, 0.0613311833066134, 0.0192528968264795, 0.0655378877483725, 0.0783480862138948, 0.0782940404757744, 0.123856592477088, 0.200715672196561, 1.05978363626061, 0.259877032381568, 0.0499868240451633, 0.081757023344427, 0.0111856453710174, 0.0354604049009183,
    0.0587501481129721, 0.0636889418869513, 0.00147409429930324, 0.0351622988914932, 0.0463883820978562, 0.049570609596332, 0.037198327424167, 0.158354051456104, 0.259877032381568, 1.03452106190489, 0.0130935177030827, 0.0516032543280675, 0.00476958310192178, 0.0237876050975194,
    0.0128020244583516, 0.00524046477601784, 0.0101454317512308, 0.0145095634258482, 0.030760164133941, 0.0163150127976701, 0.0705617007565768, 0.0128999304792121, 0.0499868240451633, 0.0130935177030827, 0.99153848322887, 0.0222928249075566, 0.054059639075802, -0.0102151447605616,
    0.00365141665784929, -0.00355524451103263, 0.0000358383973074061, 0.019673776614872, 0.0216100405840245, 0.00826771115651044, 0.00504457407814051, 0.062949163900661, 0.081757023344427, 0.0516032543280675, 0.0222928249075566, 1.01042340194426, 0.00760803953499768, 0.185604560675876,
    -0.0155671845468513, -0.00729012474606044, 0.0104003413752521, 0.00785055379131484, 0.00590293806296044, -0.00803467671394816, 0.019661904807875, -0.00201829126129822, 0.0111856453710174, 0.00476958310192178, 0.054059639075802, 0.00760803953499768, 0.989498280718446, 0.00452250880589401,
    -0.00791605733351862, 0.0104627183266632, 0.00494443990553433, 0.00996193534311699, 0.00603421242067271, 0.00567286440016974, 0.00527052695818655, 0.0298435457851606, 0.0354604049009183, 0.0237876050975194, -0.0102151447605616, 0.185604560675876, 0.00452250880589401, 1.03329693158974
  ), nrow = 14L, byrow = TRUE, dimnames = list(traits, traits))

  S <- S_reported
  min_eigen_reported <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eigen_reported <= 0) {
    S <- as.matrix(Matrix::nearPD(S, corr = FALSE, keepDiag = TRUE)$mat)
    dimnames(S) <- list(traits, traits)
  }

  V <- diag(ldsc_se[lower.tri(ldsc_se, diag = TRUE)]^2)

  out <- list(
    V = V,
    S = S,
    I = I,
    N = matrix(rep(100000, ncol(V)), nrow = 1L),
    m = NA
  )
  attr(out, "reported_min_eigen") <- min_eigen_reported
  attr(out, "smooth_max_abs_delta") <- max(abs(S - S_reported), na.rm = TRUE)
  attr(out, "traits") <- traits
  out
}

five_factor_model <- function(include_snp = FALSE) {
  base <- paste(
    "F1 =~ AN + OCD + TS + ANX",
    "F2 =~ SCZ + BIP",
    "F3 =~ ASD + ADHD + TS",
    "F4 =~ PTSD + MD + ANX",
    "F5 =~ OUD + CUD + AUD + NIC + ADHD",
    sep = "\n"
  )
  if (!include_snp) {
    return(base)
  }
  paste(
    base,
    "F1 ~ SNP",
    "F2 ~ SNP",
    "F3 ~ SNP",
    "F4 ~ SNP",
    "F5 ~ SNP",
    sep = "\n"
  )
}

public_model_check <- function(covstruc) {
  expected <- c(CFI = 0.971, SRMR = 0.063)
  timed <- time_expr({
    capture.output({
      result <- suppressWarnings(usermodel(
        covstruc = covstruc,
        model = five_factor_model(FALSE),
        std.lv = TRUE
      ))
    })
    result
  })

  fit <- unlist(timed$value$modelfit[1, names(expected)])
  data.frame(
    stage = "public_table_five_factor_model",
    backend = "GenomicSEM_usermodel_public_rounded",
    cores = 1L,
    elapsed_sec = timed$elapsed,
    n_snp = NA_integer_,
    rows = nrow(timed$value$results),
    cols = ncol(timed$value$results),
    checksum = numeric_checksum(timed$value),
    max_abs_diff_vs_old = NA_real_,
    equivalent_to_old = NA,
    max_abs_diff_vs_published = max(abs(as.numeric(fit) - expected), na.rm = TRUE),
    equivalent_to_published = max(abs(as.numeric(fit) - expected), na.rm = TRUE) < 5e-3,
    fast_path = NA_character_,
    fast_threads = NA_integer_,
    compared_numeric_cols = paste(names(expected), collapse = ","),
    note = sprintf(
      "public rounded S plus diagonal V; reported min eigen %.6g; nearPD max delta %.6g",
      attr(covstruc, "reported_min_eigen"),
      attr(covstruc, "smooth_max_abs_delta")
    ),
    stringsAsFactors = FALSE
  )
}

make_snp_fixture <- function(n_snp, covstruc, seed = 2025L) {
  set.seed(seed)
  traits <- attr(covstruc, "traits")
  maf <- runif(n_snp, 0.05, 0.49)

  snps <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    CHR = rep(11L, n_snp),
    BP = 112700000L + seq_len(n_snp),
    MAF = maf,
    A1 = rep("A", n_snp),
    A2 = rep("G", n_snp),
    check.names = FALSE
  )
  if (n_snp >= 5L) {
    snps$SNP[1:5] <- c("rs12107418", "rs9375188", "rs2970610", "rs542815", "rs7113596")
    snps$CHR[1:5] <- c(3L, 6L, 1L, 1L, 11L)
    snps$BP[1:5] <- c(48689787L, 98555272L, 44097530L, 61077199L, 112883761L)
    snps$MAF[1:5] <- c(0.210, 0.180, 0.350, 0.240, 0.120)
    snps$A1[1:5] <- c("A", "C", "G", "T", "A")
    snps$A2[1:5] <- c("G", "T", "A", "C", "C")
  }

  loading <- matrix(0, nrow = length(traits), ncol = 5L, dimnames = list(traits, paste0("F", 1:5)))
  loading["AN", "F1"] <- 0.655
  loading["OCD", "F1"] <- 0.958
  loading["TS", c("F1", "F3")] <- c(0.317, 0.343)
  loading["ANX", c("F1", "F4")] <- c(0.305, 0.766)
  loading["SCZ", "F2"] <- 0.837
  loading["BIP", "F2"] <- 0.776
  loading["ASD", "F3"] <- 0.379
  loading["ADHD", c("F3", "F5")] <- c(0.906, 0.338)
  loading["PTSD", "F4"] <- 0.824
  loading["MD", "F4"] <- 1.065
  loading["OUD", "F5"] <- 1.007
  loading["NIC", "F5"] <- 0.507
  loading["CUD", "F5"] <- 0.728
  loading["AUD", "F5"] <- 0.782

  latent_effects <- matrix(rnorm(n_snp * 5L, sd = 0.0015), n_snp, 5L)
  beta <- latent_effects %*% t(loading) + matrix(rnorm(n_snp * length(traits), sd = 0.0007), n_snp)
  se <- matrix(runif(n_snp * length(traits), min = 0.018, max = 0.04), n_snp)

  beta_df <- as.data.frame(beta, check.names = FALSE)
  names(beta_df) <- paste0("beta.", traits)
  se_df <- as.data.frame(se, check.names = FALSE)
  names(se_df) <- paste0("se.", traits)

  cbind(snps, beta_df, se_df)
}

run_usergwas_five_factor <- function(covstruc, SNPs, cores, use_new) {
  opts <- if (use_new) new_options() else old_options()
  with_options(opts, {
    timed <- time_expr({
      capture.output({
        result <- suppressWarnings(userGWAS(
          covstruc = covstruc,
          SNPs = SNPs,
          model = five_factor_model(TRUE),
          sub = paste0("F", 1:5, "~SNP"),
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

benchmark_model_workflow <- function(stage, run_fun, covstruc, SNPs, cores) {
  old <- run_fun(covstruc, SNPs, cores, FALSE)
  new <- run_fun(covstruc, SNPs, cores, TRUE)
  comparison <- compare_outputs(old$result, new$result, tolerance = 1e-4)

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
      note = NA_character_,
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
      note = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}

make_prep_fixture <- function(n_snp, traits, dir, seed = 2026L) {
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
  write.table(ref, file.path(dir, "reference.1000G.cdg2025_shape.txt"), row.names = FALSE, quote = FALSE)
  write.table(ref[c("SNP", "A1", "A2")], file.path(dir, "w_hm3.cdg2025_shape.snplist"), row.names = FALSE, quote = FALSE)

  files <- file.path(dir, paste0(traits, "_cdg2025_shape.txt"))
  for (i in seq_along(files)) {
    flip <- (seq_len(n_snp) + i) %% 7L == 0L
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
    trait_names = traits,
    ref = file.path(dir, "reference.1000G.cdg2025_shape.txt"),
    hm3 = file.path(dir, "w_hm3.cdg2025_shape.snplist")
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
      note = NA_character_,
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
      note = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

  covstruc <- make_grotzinger_2025_ldsc()
  SNPs <- make_snp_fixture(args$model_snps, covstruc)
  traits <- attr(covstruc, "traits")

  rows <- list(public_model_check(covstruc))
  for (core in args$cores) {
    rows[[length(rows) + 1L]] <- benchmark_model_workflow(
      "paper_shaped_2025_userGWAS_5factor_Q_SNP",
      run_usergwas_five_factor,
      covstruc,
      SNPs,
      core
    )
  }

  if (!args$skip_prep) {
    prep_dir <- tempfile("grotzinger-2025-prep-")
    fixture <- make_prep_fixture(args$prep_snps, traits, prep_dir)
    rows[[length(rows) + 1L]] <- benchmark_prep_workflow(
      "paper_shaped_2025_sumstats_14trait",
      benchmark_sumstats(fixture, FALSE),
      benchmark_sumstats(fixture, TRUE),
      args$prep_snps
    )
    rows[[length(rows) + 1L]] <- benchmark_prep_workflow(
      "paper_shaped_2025_munge_14trait",
      benchmark_munge(fixture, FALSE),
      benchmark_munge(fixture, TRUE),
      args$prep_snps
    )
  }

  results <- do.call(rbind, rows)
  csv_path <- file.path(args$out_dir, paste0("grotzinger_2025_nature_", timestamp, ".csv"))
  write.csv(results, csv_path, row.names = FALSE)

  cat("Wrote ", csv_path, "\n", sep = "")
  print(results, row.names = FALSE)
}

if (sys.nframe() == 0L) {
  main()
}
