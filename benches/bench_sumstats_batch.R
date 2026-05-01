args <- commandArgs(trailingOnly = TRUE)
n_snp <- if (length(args) >= 1) as.integer(args[[1]]) else 100000L
n_traits <- if (length(args) >= 2) as.integer(args[[2]]) else 4L
n_threads <- if (length(args) >= 3) as.integer(args[[3]]) else min(n_traits, 4L)
seed <- if (length(args) >= 4) as.integer(args[[4]]) else 1L

if (requireNamespace("GenomicSEM", quietly = TRUE)) {
  library(GenomicSEM)
} else {
  stop("GenomicSEM must be installed before running this benchmark", call. = FALSE)
}

make_ref <- function(n) {
  alleles <- c("A", "C", "G", "T")
  a1 <- sample(alleles, n, replace = TRUE)
  a2 <- sample(alleles, n, replace = TRUE)
  same <- a1 == a2
  a2[same] <- alleles[(match(a1[same], alleles) %% length(alleles)) + 1L]

  data.frame(
    SNP = paste0("rs", seq_len(n)),
    A1 = a1,
    A2 = a2,
    MAF = runif(n, 0.02, 0.49),
    stringsAsFactors = FALSE
  )
}

make_trait <- function(ref, trait_index) {
  n <- nrow(ref)
  flip <- (seq_len(n) + trait_index) %% 4L == 0L
  a1 <- ref$A1
  a2 <- ref$A2
  a1[flip] <- ref$A2[flip]
  a2[flip] <- ref$A1[flip]

  out <- data.frame(
    SNP = ref$SNP,
    A1 = a1,
    A2 = a2,
    BETA = rnorm(n, sd = 0.05),
    SE = runif(n, 0.02, 0.08),
    P = format.pval(runif(n, 1e-8, 0.999), digits = 8),
    N = sample(8000:12000, n, replace = TRUE),
    INFO = runif(n, 0.91, 1.0),
    MAF = ref$MAF,
    stringsAsFactors = FALSE
  )
  out$P[seq.int(trait_index, n, by = max(1L, n_traits * 2000L))] <- NA_character_
  out
}

with_temp_cwd <- function(code) {
  old <- getwd()
  tmp <- tempfile("genomicssem-sumstats-batch-")
  dir.create(tmp)
  on.exit({
    setwd(old)
    unlink(tmp, recursive = TRUE)
  })
  setwd(tmp)
  force(code)
}

time_expr <- function(expr) {
  gc()
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = proc.time()[["elapsed"]] - start)
}

checksum_output <- function(out) {
  numeric_cols <- vapply(out, is.numeric, logical(1))
  sum(as.matrix(out[, numeric_cols, drop = FALSE]), na.rm = TRUE) + nrow(out)
}

with_temp_cwd({
  set.seed(seed)
  ref <- make_ref(n_snp)
  write.table(ref, "ref.txt", row.names = FALSE, quote = FALSE)

  files <- paste0("trait", seq_len(n_traits), ".txt")
  for (i in seq_len(n_traits)) {
    write.table(make_trait(ref, i), files[i], row.names = FALSE, quote = FALSE)
  }

  old <- options(
    GenomicSEM.fast_sumstats_engine = getOption("GenomicSEM.fast_sumstats_engine"),
    GenomicSEM.fast_sumstats_threads = getOption("GenomicSEM.fast_sumstats_threads"),
    GenomicSEM.fast_table_read = getOption("GenomicSEM.fast_table_read"),
    GenomicSEM.fast_sumstats_qc = getOption("GenomicSEM.fast_sumstats_qc")
  )
  on.exit(options(old), add = TRUE)

  modes <- data.frame(
    backend = c("legacy_serial", "legacy_parallel", "native_batch_1t", "native_batch_threads"),
    fast_engine = c(FALSE, FALSE, TRUE, TRUE),
    parallel = c(FALSE, TRUE, FALSE, TRUE),
    threads = c(1L, n_threads, 1L, n_threads),
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(modes)), function(i) {
    mode <- modes[i, ]
    options(
      GenomicSEM.fast_sumstats_engine = mode$fast_engine,
      GenomicSEM.fast_sumstats_threads = mode$threads,
      GenomicSEM.fast_table_read = FALSE,
      GenomicSEM.fast_sumstats_qc = FALSE
    )
    timed <- time_expr({
      capture.output(out <- suppressWarnings(sumstats(
        files = files,
        ref = "ref.txt",
        trait.names = paste0("trait", seq_len(n_traits)),
        se.logit = rep(TRUE, n_traits),
        parallel = mode$parallel,
        cores = mode$threads
      )))
      out
    })
    data.frame(
      backend = mode$backend,
      n_snp = n_snp,
      n_traits = n_traits,
      threads = mode$threads,
      elapsed_sec = timed$elapsed,
      rows = nrow(timed$value),
      checksum = checksum_output(timed$value),
      stringsAsFactors = FALSE
    )
  })

  print(do.call(rbind, rows), row.names = FALSE)
})
