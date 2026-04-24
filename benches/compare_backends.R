if (requireNamespace("GenomicSEM", quietly = TRUE)) {
  library(GenomicSEM)
  genomicssem_ns <- asNamespace("GenomicSEM")
} else {
  dyn.load("GenomicSEM.so")
  source("R/utils.R")
  genomicssem_ns <- globalenv()
}

.get_V_SNP <- get(".get_V_SNP", envir = genomicssem_ns)
.get_V_SNP_r <- get(".get_V_SNP_r", envir = genomicssem_ns)
.get_V_SNP_batch <- get(".get_V_SNP_batch", envir = genomicssem_ns)

args <- commandArgs(trailingOnly = TRUE)
n_snp <- if (length(args) >= 1) as.integer(args[[1]]) else 200000L
k <- if (length(args) >= 2) as.integer(args[[2]]) else 8L
gc_mode <- if (length(args) >= 3) args[[3]] else "standard"
thread_counts <- if (length(args) >= 4) {
  as.integer(strsplit(args[[4]], ",", fixed = TRUE)[[1]])
} else {
  c(1L, 2L, 4L, 8L, 16L)
}
thread_counts <- sort(unique(thread_counts[thread_counts > 0]))

set.seed(1)
SE_SNP <- matrix(runif(n_snp * k, 0.01, 0.5), n_snp, k)
I_LD <- matrix(runif(k * k, 0.05, 0.4), k, k)
I_LD <- 0.5 * (I_LD + t(I_LD))
diag(I_LD) <- runif(k, 1.0, 1.8)
coords <- which(!is.na(I_LD), arr.ind = TRUE)
varSNP <- runif(n_snp, 0.05, 0.95)

time_it <- function(expr) {
  gc()
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(elapsed = proc.time()[["elapsed"]] - start, value = value)
}

bench_r_loop <- function(use_rust) {
  options(GenomicSEM.use_rust = use_rust)
  result <- time_it({
    checksum <- 0
    for (i in seq_len(n_snp)) {
      v <- .get_V_SNP(SE_SNP, I_LD, varSNP, gc_mode, coords, k, i)
      checksum <- checksum + v[1, 1] + v[k, k]
    }
    checksum
  })
  data.frame(
    backend = if (use_rust) "r_binding_loop" else "old_r_loop",
    threads = 1L,
    n_snp = n_snp,
    k = k,
    gc = gc_mode,
    elapsed_sec = result$elapsed,
    checksum = result$value,
    data_source = "R"
  )
}

bench_r_batch <- function(n_threads) {
  options(GenomicSEM.use_rust = TRUE)
  result <- time_it({
    out <- .get_V_SNP_batch(SE_SNP, I_LD, varSNP, gc_mode, coords, k, n_threads = n_threads)
    sum(out[1, 1, ] + out[k, k, ])
  })
  data.frame(
    backend = "r_binding_batch",
    threads = n_threads,
    n_snp = n_snp,
    k = k,
    gc = gc_mode,
    elapsed_sec = result$elapsed,
    checksum = result$value,
    data_source = "R"
  )
}

bench_rust_cli <- function() {
  cmd <- c("run", "--release", "-q", "-p", "genomicssem-core", "--bin", "kernel_bench", "--", n_snp, k)
  output <- system2("cargo", cmd, stdout = TRUE, stderr = TRUE)
  csv <- output[grepl("^(backend,|rust_)", output)]
  if (length(csv) < 2) {
    warning("Could not parse Rust CLI benchmark output:\n", paste(output, collapse = "\n"))
    return(data.frame())
  }

  parsed <- read.csv(text = paste(csv, collapse = "\n"))
  parsed$gc <- gc_mode
  parsed$data_source <- "rust_cli"
  parsed <- parsed[, c("backend", "threads", "n_snp", "k", "gc", "elapsed_sec", "checksum", "data_source")]
  parsed
}

results <- do.call(
  rbind,
  c(
    list(bench_r_loop(FALSE), bench_r_loop(TRUE)),
    lapply(thread_counts, bench_r_batch),
    list(bench_rust_cli())
  )
)

print(results, row.names = FALSE)
