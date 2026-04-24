library(GenomicSEM)

genomicssem_ns <- asNamespace("GenomicSEM")
.get_V_SNP <- get(".get_V_SNP", envir = genomicssem_ns)
.get_V_SNP_r <- get(".get_V_SNP_r", envir = genomicssem_ns)
.get_V_SNP_batch <- get(".get_V_SNP_batch", envir = genomicssem_ns)

args <- commandArgs(trailingOnly = TRUE)
n_snp <- if (length(args) >= 1) as.integer(args[[1]]) else 200000L
k <- if (length(args) >= 2) as.integer(args[[2]]) else 8L
gc_mode <- if (length(args) >= 3) args[[3]] else "standard"

set.seed(1)
SE_SNP <- matrix(runif(n_snp * k, 0.01, 0.5), n_snp, k)
I_LD <- matrix(runif(k * k, 0.05, 0.4), k, k)
I_LD <- 0.5 * (I_LD + t(I_LD))
diag(I_LD) <- runif(k, 1.0, 1.8)
coords <- which(!is.na(I_LD), arr.ind = TRUE)
varSNP <- runif(n_snp, 0.05, 0.95)

bench <- function(use_rust) {
  options(GenomicSEM.use_rust = use_rust)
  gc()
  start <- proc.time()[["elapsed"]]
  checksum <- 0
  for (i in seq_len(n_snp)) {
    v <- .get_V_SNP(SE_SNP, I_LD, varSNP, gc_mode, coords, k, i)
    checksum <- checksum + v[1, 1] + v[k, k]
  }
  elapsed <- proc.time()[["elapsed"]] - start
  data.frame(
    backend = if (use_rust) "rust" else "r",
    n_snp = n_snp,
    k = k,
    gc = gc_mode,
    elapsed_sec = elapsed,
    checksum = checksum
  )
}

print(rbind(bench(FALSE), bench(TRUE)), row.names = FALSE)

if (is.loaded("genomicssem_get_v_snp_batch_call")) {
  thread_counts <- sort(unique(pmin(c(1L, 2L, 4L, 8L, 16L), parallel::detectCores())))
  batch_results <- lapply(thread_counts, function(n_threads) {
    gc()
    start <- proc.time()[["elapsed"]]
    out <- .get_V_SNP_batch(SE_SNP, I_LD, varSNP, gc_mode, coords, k, n_threads = n_threads)
    elapsed <- proc.time()[["elapsed"]] - start
    data.frame(
      backend = "rust_batch",
      n_threads = n_threads,
      n_snp = n_snp,
      k = k,
      gc = gc_mode,
      elapsed_sec = elapsed,
      checksum = sum(out[1, 1, ] + out[k, k, ])
    )
  })
  print(do.call(rbind, batch_results), row.names = FALSE)
}
