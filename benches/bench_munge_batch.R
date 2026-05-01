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

  data.frame(
    SNP = ref$SNP,
    A1 = a1,
    A2 = a2,
    BETA = rnorm(n, sd = 0.05),
    P = format.pval(runif(n, 1e-8, 0.999), digits = 8),
    N = sample(8000:12000, n, replace = TRUE),
    INFO = runif(n, 0.91, 1.0),
    MAF = runif(n, 0.02, 0.49),
    stringsAsFactors = FALSE
  )
}

with_temp_cwd <- function(code) {
  old <- getwd()
  tmp <- tempfile("genomicssem-munge-batch-")
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

checksum_outputs <- function(trait_names) {
  sum(vapply(trait_names, function(trait_name) {
    out <- read.table(gzfile(paste0(trait_name, ".sumstats.gz")), header = TRUE)
    sum(out$N, out$Z, na.rm = TRUE) + nrow(out)
  }, numeric(1)))
}

with_temp_cwd({
  set.seed(seed)
  ref <- make_ref(n_snp)
  write.table(ref, "hm3.txt", row.names = FALSE, quote = FALSE)

  files <- paste0("trait", seq_len(n_traits), ".txt")
  for (i in seq_len(n_traits)) {
    write.table(make_trait(ref, i), files[i], row.names = FALSE, quote = FALSE)
  }

  old <- options(
    GenomicSEM.fast_munge_engine = getOption("GenomicSEM.fast_munge_engine"),
    GenomicSEM.fast_munge_threads = getOption("GenomicSEM.fast_munge_threads"),
    GenomicSEM.fast_table_read = getOption("GenomicSEM.fast_table_read"),
    GenomicSEM.fast_munge_qc = getOption("GenomicSEM.fast_munge_qc")
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
      GenomicSEM.fast_munge_engine = mode$fast_engine,
      GenomicSEM.fast_munge_threads = mode$threads,
      GenomicSEM.fast_table_read = FALSE,
      GenomicSEM.fast_munge_qc = FALSE
    )
    trait_names <- paste0(mode$backend, "_trait", seq_len(n_traits))
    timed <- time_expr({
      capture.output(suppressWarnings(munge(
        files = files,
        hm3 = "hm3.txt",
        trait.names = trait_names,
        parallel = mode$parallel,
        cores = mode$threads,
        overwrite = TRUE,
        column.names = list(effect = "BETA")
      )))
    })
    checksum <- checksum_outputs(trait_names)
    data.frame(
      backend = mode$backend,
      n_snp = n_snp,
      n_traits = n_traits,
      threads = mode$threads,
      elapsed_sec = timed$elapsed,
      checksum = checksum,
      stringsAsFactors = FALSE
    )
  })

  print(do.call(rbind, rows), row.names = FALSE)
})
