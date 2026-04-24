args <- commandArgs(trailingOnly = TRUE)
n_snp <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
k <- if (length(args) >= 2) as.integer(args[[2]]) else 6L
cores <- if (length(args) >= 3) {
  as.integer(strsplit(args[[3]], ",", fixed = TRUE)[[1]])
} else {
  c(1L, 2L, 4L, 8L, 16L)
}
workflow <- if (length(args) >= 4) args[[4]] else "userGWAS"
seed <- if (length(args) >= 5) as.integer(args[[5]]) else 1L
fast_diag_modes <- if (length(args) >= 6) {
  strsplit(args[[6]], ",", fixed = TRUE)[[1]]
} else {
  "TRUE"
}
fast_wls_modes <- if (length(args) >= 7) {
  strsplit(args[[7]], ",", fixed = TRUE)[[1]]
} else {
  "FALSE"
}
printwarn_modes <- if (length(args) >= 8) {
  strsplit(args[[8]], ",", fixed = TRUE)[[1]]
} else {
  "TRUE"
}
q_snp_modes <- if (length(args) >= 9) {
  strsplit(args[[9]], ",", fixed = TRUE)[[1]]
} else {
  "TRUE"
}
fast_commonfactor_modes <- if (length(args) >= 10) {
  strsplit(args[[10]], ",", fixed = TRUE)[[1]]
} else {
  "FALSE"
}
fast_usergwas_modes <- if (length(args) >= 11) {
  strsplit(args[[11]], ",", fixed = TRUE)[[1]]
} else {
  "FALSE"
}

cores <- sort(unique(cores[cores > 0]))
workflow <- match.arg(workflow, c("userGWAS", "commonfactorGWAS"))
fast_diag_modes <- fast_diag_modes %in% c("TRUE", "true", "1", "yes", "fast")
fast_wls_modes <- fast_wls_modes %in% c("TRUE", "true", "1", "yes", "fast")
printwarn_modes <- printwarn_modes %in% c("TRUE", "true", "1", "yes", "warn")
q_snp_modes <- q_snp_modes %in% c("TRUE", "true", "1", "yes", "q")
fast_commonfactor_modes <- fast_commonfactor_modes %in% c("TRUE", "true", "1", "yes", "fast")
fast_usergwas_modes <- fast_usergwas_modes %in% c("TRUE", "true", "1", "yes", "fast")
if (workflow == "commonfactorGWAS") {
  printwarn_modes <- NA
  fast_usergwas_modes <- NA
} else {
  fast_commonfactor_modes <- NA
}

if (requireNamespace("GenomicSEM", quietly = TRUE)) {
  library(GenomicSEM)
} else {
  dyn.load("GenomicSEM.so")
  for (file in list.files("R", full.names = TRUE, pattern = "\\.[Rr]$")) {
    source(file)
  }
}

source("benches/synthetic_inputs.R")
inputs <- make_synthetic_genomicsem_inputs(n_snp = n_snp, k = k, seed = seed)

run_workflow <- function(use_rust, n_core, fast_diag, fast_wls, printwarn, q_snp, fast_commonfactor, fast_usergwas) {
  options(GenomicSEM.use_rust = use_rust)
  options(GenomicSEM.fast_diag_inverse = fast_diag)
  options(GenomicSEM.fast_diagonal_wls = fast_wls)
  options(GenomicSEM.fast_commonfactor_fit = fast_commonfactor)
  options(GenomicSEM.fast_usergwas_fit = fast_usergwas)

  parallel <- n_core > 1L
  result <- NULL
  gc()
  start <- proc.time()[["elapsed"]]
  output <- capture.output({
    result <- suppressWarnings(
      if (workflow == "userGWAS") {
        userGWAS(
          covstruc = inputs$covstruc,
          SNPs = inputs$SNPs,
          model = inputs$model,
          estimation = "DWLS",
          parallel = parallel,
          cores = if (parallel) n_core else NULL,
          GC = "standard",
          fix_measurement = TRUE,
          Q_SNP = q_snp,
          printwarn = printwarn
        )
      } else {
        commonfactorGWAS(
          covstruc = inputs$covstruc,
          SNPs = inputs$SNPs,
          estimation = "DWLS",
          parallel = parallel,
          cores = if (parallel) n_core else NULL,
          GC = "standard",
          Q_SNP = q_snp
        )
      }
    )
  })
  elapsed <- proc.time()[["elapsed"]] - start

  data.frame(
    workflow = workflow,
    backend = if (use_rust) "rust_binding_workflow" else "old_r_workflow",
    cores = n_core,
    parallel = parallel,
    fast_diag_inverse = fast_diag,
    fast_diagonal_wls = fast_wls,
    printwarn = if (workflow == "userGWAS") printwarn else NA,
    Q_SNP = q_snp,
    fast_commonfactor_fit = if (workflow == "commonfactorGWAS") fast_commonfactor else NA,
    fast_usergwas_fit = if (workflow == "userGWAS") fast_usergwas else NA,
    n_snp = n_snp,
    k = k,
    elapsed_sec = elapsed,
    checksum = checksum_result(result),
    result_rows = count_result_rows(result),
    result_cols = count_result_cols(result),
    result_size_bytes = as.numeric(object.size(result)),
    stdout_lines = length(output),
    stringsAsFactors = FALSE
  )
}

results <- do.call(
  rbind,
  lapply(fast_diag_modes, function(fast_diag) {
    do.call(rbind, lapply(fast_wls_modes, function(fast_wls) {
      do.call(rbind, lapply(printwarn_modes, function(printwarn) {
        do.call(rbind, lapply(q_snp_modes, function(q_snp) {
          do.call(rbind, lapply(fast_commonfactor_modes, function(fast_commonfactor) {
            do.call(rbind, lapply(fast_usergwas_modes, function(fast_usergwas) {
              do.call(
                rbind,
                c(
                  lapply(cores, function(n_core) run_workflow(FALSE, n_core, fast_diag, fast_wls, printwarn, q_snp, fast_commonfactor, fast_usergwas)),
                  lapply(cores, function(n_core) run_workflow(TRUE, n_core, fast_diag, fast_wls, printwarn, q_snp, fast_commonfactor, fast_usergwas))
                )
              )
            }))
          }))
        }))
      }))
    }))
  })
)

print(results, row.names = FALSE)
