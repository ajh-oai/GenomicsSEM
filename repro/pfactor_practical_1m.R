#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicSEM))
suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

source_manifest <- function() {
  data.frame(
    artifact = c("SCZ", "BIP", "MDD", "GenomicSEMPractical.RData"),
    url = c(
      "https://ndownloader.figshare.com/files/28198983",
      "https://ndownloader.figshare.com/files/28169301",
      "https://ndownloader.figshare.com/files/28169508",
      "https://utexas.app.box.com/index.php?rm=box_download_shared_file&shared_name=sounavy84gwygj0j2askcyaoo2ostbgu&file_id=f_819725878597"
    ),
    filename = c(
      "CLOZUK_PGC2noclo.METAL.assoc.dosage.fix.gz",
      "pgc.bip.full.2012-04.txt.gz",
      "MDD2018_ex23andMe.gz",
      "GenomicSEMPractical.RData"
    ),
    stringsAsFactors = FALSE
  )
}

parse_args <- function(args) {
  detected_cores <- parallel::detectCores(logical = TRUE)
  if (!is.finite(detected_cores)) {
    detected_cores <- 1L
  }

  out <- list(
    target_snps = 1000000L,
    candidate_multiplier = 1.35,
    data_dir = file.path("repro", "data", "pfactor_practical"),
    out_dir = file.path("repro", "results"),
    cores = c(1L, 4L, 16L),
    old_gwas_snps = 1000000L,
    q_snp = FALSE,
    skip_download = FALSE,
    reuse_subset = FALSE,
    reuse_sumstats = FALSE,
    threads = max(1L, detected_cores),
    save_outputs = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else NA_character_

    if (identical(key, "--target-snps")) {
      out$target_snps <- as.integer(value)
      i <- i + 2L
    } else if (identical(key, "--candidate-multiplier")) {
      out$candidate_multiplier <- as.numeric(value)
      i <- i + 2L
    } else if (identical(key, "--data-dir")) {
      out$data_dir <- value
      i <- i + 2L
    } else if (identical(key, "--out-dir")) {
      out$out_dir <- value
      i <- i + 2L
    } else if (identical(key, "--cores")) {
      out$cores <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
      i <- i + 2L
    } else if (identical(key, "--old-gwas-snps")) {
      out$old_gwas_snps <- as.integer(value)
      i <- i + 2L
    } else if (identical(key, "--threads")) {
      out$threads <- as.integer(value)
      i <- i + 2L
    } else if (identical(key, "--q-snp")) {
      out$q_snp <- TRUE
      i <- i + 1L
    } else if (identical(key, "--skip-download")) {
      out$skip_download <- TRUE
      i <- i + 1L
    } else if (identical(key, "--reuse-subset")) {
      out$reuse_subset <- TRUE
      i <- i + 1L
    } else if (identical(key, "--reuse-sumstats")) {
      out$reuse_sumstats <- TRUE
      i <- i + 1L
    } else if (identical(key, "--save-outputs")) {
      out$save_outputs <- TRUE
      i <- i + 1L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }

  out$target_snps <- as.integer(out$target_snps)
  out$old_gwas_snps <- as.integer(out$old_gwas_snps)
  out$threads <- max(1L, as.integer(out$threads))
  out$cores <- sort(unique(out$cores[is.finite(out$cores) & out$cores > 0L]))
  if (length(out$cores) == 0L) {
    out$cores <- 1L
  }
  if (!is.finite(out$candidate_multiplier) || out$candidate_multiplier < 1) {
    stop("--candidate-multiplier must be at least 1", call. = FALSE)
  }

  out
}

time_expr <- function(expr) {
  gc()
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = proc.time()[["elapsed"]] - start)
}

with_options <- function(opts, code) {
  old <- options(opts)
  on.exit(options(old), add = TRUE)
  force(code)
}

with_cwd <- function(path, code) {
  old <- getwd()
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  on.exit(setwd(old), add = TRUE)
  setwd(path)
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

numeric_checksum <- function(x) {
  values <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  sum(values[is.finite(values)], na.rm = TRUE)
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
      value <- suppressWarnings(max(abs(a - b), na.rm = TRUE))
      if (is.finite(value)) value else 0
    }, old_df[numeric_cols], new_df[numeric_cols])
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

download_sources <- function(data_dir, skip_download, require_raw = TRUE) {
  manifest <- source_manifest()
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(manifest))) {
    dest <- file.path(data_dir, manifest$filename[[i]])
    if (file.exists(dest) && file.info(dest)$size > 0) {
      next
    }
    if (!require_raw && manifest$artifact[[i]] != "GenomicSEMPractical.RData") {
      next
    }
    if (skip_download) {
      stop("Missing ", dest, " and --skip-download was set", call. = FALSE)
    }
    message("Downloading ", manifest$artifact[[i]], " -> ", dest)
    utils::download.file(manifest$url[[i]], dest, mode = "wb", quiet = FALSE)
  }

  paths <- setNames(file.path(data_dir, manifest$filename), manifest$artifact)
  attr(paths, "manifest") <- manifest
  paths
}

valid_allele <- function(x) {
  nchar(x) == 1L & toupper(x) %chin% c("A", "C", "G", "T")
}

dedupe_first <- function(dt, snp_col = "SNP") {
  dt[!duplicated(dt[[snp_col]])]
}

read_bip_by_position <- function(path, select, names, threads) {
  out <- fread(
    path,
    skip = 1L,
    header = FALSE,
    select = select,
    nThread = threads,
    fill = TRUE,
    showProgress = TRUE
  )
  setnames(out, names)
  out
}

prepare_public_subset <- function(paths, args) {
  subset_dir <- file.path(args$data_dir, sprintf("subset_%s", args$target_snps))
  dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)

  out <- list(
    ref = file.path(subset_dir, "reference.1000G.public_pfactor.txt"),
    files = c(
      SCZ = file.path(subset_dir, "SCZ_public_pfactor.txt.gz"),
      BIP = file.path(subset_dir, "BIP_public_pfactor.txt.gz"),
      MDD = file.path(subset_dir, "MDD_public_pfactor.txt.gz")
    ),
    metadata = file.path(subset_dir, "subset_metadata.rds")
  )

  if (args$reuse_subset && file.exists(out$ref) && all(file.exists(out$files)) && file.exists(out$metadata)) {
    return(out)
  }

  setDTthreads(args$threads)
  message("Reading public SNP IDs and BIP allele-frequency reference columns")
  bip_ref <- read_bip_by_position(paths[["BIP"]], c(1L, 4L, 5L, 11L), c("SNP", "A1", "A2", "MAF"), args$threads)
  bip_ref[, `:=`(A1 = toupper(A1), A2 = toupper(A2), MAF = suppressWarnings(as.numeric(MAF)))]
  bip_ref <- bip_ref[
    valid_allele(A1) & valid_allele(A2) & A1 != A2 & is.finite(MAF) & MAF >= 0.01 & MAF <= 0.99
  ]
  bip_ref <- dedupe_first(bip_ref)

  scz_ids <- fread(paths[["SCZ"]], select = "SNP", nThread = args$threads, fill = TRUE, showProgress = TRUE)
  mdd_ids <- fread(paths[["MDD"]], select = "SNP", nThread = args$threads, fill = TRUE, showProgress = TRUE)

  common <- intersect(intersect(bip_ref$SNP, scz_ids$SNP), mdd_ids$SNP)
  if (length(common) < args$target_snps) {
    stop("Only ", length(common), " SNPs overlap across public inputs; target was ", args$target_snps, call. = FALSE)
  }

  candidate_n <- min(length(common), ceiling(args$target_snps * args$candidate_multiplier))
  ref <- bip_ref[SNP %chin% common][seq_len(candidate_n)]
  ref <- ref[, .(SNP, A1, A2, MAF)]
  fwrite(ref, out$ref, sep = "\t")

  keep <- ref[, .(SNP)]
  write_compact <- function(dt, dest) {
    dt <- dedupe_first(dt)
    dt <- dt[match(keep$SNP, dt$SNP)]
    dt <- dt[!is.na(SNP)]
    fwrite(dt, dest, sep = "\t", compress = "gzip")
    invisible(nrow(dt))
  }

  message("Writing compact SCZ/BIP/MDD public-data subsets")
  scz <- fread(paths[["SCZ"]], select = c("SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P"), nThread = args$threads, fill = TRUE)
  scz[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  scz_rows <- write_compact(scz, out$files[["SCZ"]])

  bip <- read_bip_by_position(
    paths[["BIP"]],
    c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 11L),
    c("SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "INFO", "MAF"),
    args$threads
  )
  bip[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  bip_rows <- write_compact(bip, out$files[["BIP"]])

  mdd <- fread(paths[["MDD"]], select = c("SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "INFO"), nThread = args$threads, fill = TRUE)
  mdd[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  mdd_rows <- write_compact(mdd, out$files[["MDD"]])

  saveRDS(list(
    created_at = as.character(Sys.time()),
    target_snps = args$target_snps,
    candidate_snps = nrow(ref),
    overlap_snps = length(common),
    compact_rows = c(SCZ = scz_rows, BIP = bip_rows, MDD = mdd_rows),
    sources = attr(paths, "manifest")
  ), out$metadata)

  out
}

load_practical_covstruc <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env$PSYCH_COV
}

run_sumstats <- function(fixture, data_dir, use_new) {
  opts <- if (use_new) new_options() else old_options()
  work_dir <- file.path(data_dir, if (use_new) "sumstats_new" else "sumstats_old")
  with_cwd(work_dir, {
    with_options(opts, {
      time_expr({
        capture.output({
          result <- suppressWarnings(sumstats(
            files = normalizePath(unname(fixture$files)),
            ref = normalizePath(fixture$ref),
            trait.names = c("SCZ", "BIP", "MDD"),
            se.logit = c(TRUE, TRUE, TRUE),
            parallel = FALSE
          ))
        })
        result
      })
    })
  })
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

run_usergwas <- function(covstruc, snps, cores, use_new, q_snp) {
  opts <- if (use_new) new_options() else old_options()
  with_options(opts, {
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
    fast_path = attr(timed$value, "GenomicSEM.fast_path") %||% NA_character_,
    fast_threads = attr(timed$value, "GenomicSEM.fast_threads") %||% NA_integer_
  )
}

summary_row <- function(stage, backend, cores, elapsed, result, n_snp, comparison = NULL, fast_path = NA_character_, fast_threads = NA_integer_) {
  data.frame(
    stage = stage,
    backend = backend,
    cores = cores,
    elapsed_sec = elapsed,
    n_snp = n_snp,
    rows = count_rows(result),
    cols = count_cols(result),
    checksum = numeric_checksum(result),
    max_abs_diff_vs_old = if (is.null(comparison)) NA_real_ else comparison$max_abs_diff,
    equivalent_to_old = if (is.null(comparison)) NA else comparison$equivalent,
    fast_path = fast_path,
    fast_threads = fast_threads,
    compared_numeric_cols = if (is.null(comparison)) NA_character_ else comparison$numeric_cols,
    stringsAsFactors = FALSE
  )
}

benchmark_sumstats_pair <- function(fixture, args) {
  old_cache <- file.path(args$data_dir, sprintf("sumstats_old_%s.rds", args$target_snps))
  new_cache <- file.path(args$data_dir, sprintf("sumstats_new_%s.rds", args$target_snps))

  old <- if (args$reuse_sumstats && file.exists(old_cache)) {
    readRDS(old_cache)
  } else {
    result <- run_sumstats(fixture, args$data_dir, FALSE)
    saveRDS(result, old_cache)
    result
  }

  new <- if (args$reuse_sumstats && file.exists(new_cache)) {
    readRDS(new_cache)
  } else {
    result <- run_sumstats(fixture, args$data_dir, TRUE)
    saveRDS(result, new_cache)
    result
  }

  comparison <- compare_outputs(old$value, new$value, tolerance = 1e-10)
  rows <- rbind(
    summary_row("public_pfactor_sumstats", "old_r", 1L, old$elapsed, old$value, nrow(old$value), list(max_abs_diff = 0, equivalent = TRUE, numeric_cols = comparison$numeric_cols)),
    summary_row("public_pfactor_sumstats", "new_rust", 1L, new$elapsed, new$value, nrow(new$value), comparison)
  )

  list(old = old, new = new, rows = rows)
}

benchmark_usergwas_pair <- function(covstruc, snps, args, core) {
  n_old <- min(args$old_gwas_snps, nrow(snps))
  old_snps <- snps[seq_len(n_old), , drop = FALSE]

  old <- run_usergwas(covstruc, old_snps, core, FALSE, args$q_snp)
  new_compare <- run_usergwas(covstruc, old_snps, core, TRUE, args$q_snp)
  comparison <- compare_outputs(old$result, new_compare$result, tolerance = 1e-5)

  new_full <- NULL
  full_comparison <- NULL
  if (n_old == nrow(snps)) {
    new_full <- new_compare
    full_comparison <- comparison
  } else {
    new_full <- run_usergwas(covstruc, snps, core, TRUE, args$q_snp)
  }

  rows <- rbind(
    summary_row("public_pfactor_userGWAS_compare", "old_r_lavaan", core, old$elapsed, old$result, n_old, list(max_abs_diff = 0, equivalent = TRUE, numeric_cols = comparison$numeric_cols), old$fast_path, old$fast_threads),
    summary_row("public_pfactor_userGWAS_compare", "new_rust_binding", core, new_compare$elapsed, new_compare$result, n_old, comparison, new_compare$fast_path, new_compare$fast_threads),
    summary_row("public_pfactor_userGWAS_full", "new_rust_binding", core, new_full$elapsed, new_full$result, nrow(snps), full_comparison, new_full$fast_path, new_full$fast_threads)
  )

  list(old = old, new_compare = new_compare, new_full = new_full, rows = rows)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  args$data_dir <- normalizePath(args$data_dir, mustWork = FALSE)
  args$out_dir <- normalizePath(args$out_dir, mustWork = FALSE)
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

  paths <- download_sources(args$data_dir, args$skip_download, require_raw = !args$reuse_subset)
  fixture <- prepare_public_subset(paths, args)
  covstruc <- load_practical_covstruc(paths[["GenomicSEMPractical.RData"]])

  sumstats_result <- benchmark_sumstats_pair(fixture, args)
  snps <- sumstats_result$new$value
  if (nrow(snps) < args$target_snps) {
    stop("sumstats produced only ", nrow(snps), " SNPs after QC; target was ", args$target_snps, call. = FALSE)
  }
  snps <- snps[seq_len(args$target_snps), , drop = FALSE]

  rows <- list(sumstats_result$rows)
  usergwas_results <- list()
  for (core in args$cores) {
    result <- benchmark_usergwas_pair(covstruc, snps, args, core)
    rows[[length(rows) + 1L]] <- result$rows
    usergwas_results[[as.character(core)]] <- result
  }

  results <- do.call(rbind, rows)
  results$target_snps <- args$target_snps
  results$old_gwas_snps <- args$old_gwas_snps
  results$q_snp <- args$q_snp

  csv_path <- file.path(args$out_dir, paste0("pfactor_practical_1m_", timestamp, ".csv"))
  write.csv(results, csv_path, row.names = FALSE)

  source_path <- file.path(args$out_dir, paste0("pfactor_practical_1m_sources_", timestamp, ".csv"))
  write.csv(attr(paths, "manifest"), source_path, row.names = FALSE)

  if (args$save_outputs) {
    saveRDS(list(
      sumstats = sumstats_result,
      usergwas = usergwas_results,
      snps = snps,
      covstruc = covstruc
    ), file.path(args$out_dir, paste0("pfactor_practical_1m_outputs_", timestamp, ".rds")))
  }

  cat("Wrote ", csv_path, "\n", sep = "")
  cat("Wrote ", source_path, "\n", sep = "")
  print(results, row.names = FALSE)
}

main()
