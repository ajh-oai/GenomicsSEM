args <- commandArgs(trailingOnly = TRUE)
n_snp <- if (length(args) >= 1) as.integer(args[[1]]) else 100L
k <- if (length(args) >= 2) as.integer(args[[2]]) else 12L
workflow <- if (length(args) >= 3) args[[3]] else "commonfactorGWAS"
seed <- if (length(args) >= 4) as.integer(args[[4]]) else 1L
out_dir <- if (length(args) >= 5) args[[5]] else "benches/profiles"
memory_profile <- if (length(args) >= 6) {
  args[[6]] %in% c("TRUE", "true", "1", "yes", "memory")
} else {
  FALSE
}

workflow <- match.arg(workflow, c("userGWAS", "commonfactorGWAS"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(GenomicSEM)
source("benches/synthetic_inputs.R")

profile_path <- file.path(out_dir, paste0(workflow, "_n", n_snp, "_k", k, ".out"))
summary_path <- file.path(out_dir, paste0(workflow, "_n", n_snp, "_k", k, "_summary.txt"))
self_csv_path <- file.path(out_dir, paste0(workflow, "_n", n_snp, "_k", k, "_self.csv"))
total_csv_path <- file.path(out_dir, paste0(workflow, "_n", n_snp, "_k", k, "_total.csv"))

inputs <- make_synthetic_genomicsem_inputs(n_snp = n_snp, k = k, seed = seed)
options(GenomicSEM.use_rust = TRUE)
options(GenomicSEM.fast_diag_inverse = TRUE)
options(GenomicSEM.fast_diagonal_wls = FALSE)

result <- NULL
invisible(gc())
start <- proc.time()[["elapsed"]]
Rprof(profile_path, interval = 0.001, memory.profiling = memory_profile)
output <- capture.output({
  result <- suppressWarnings(
    if (workflow == "userGWAS") {
      userGWAS(
        covstruc = inputs$covstruc,
        SNPs = inputs$SNPs,
        model = inputs$model,
        estimation = "DWLS",
        parallel = FALSE,
        GC = "standard",
        fix_measurement = TRUE,
        Q_SNP = FALSE,
        printwarn = FALSE
      )
    } else {
      commonfactorGWAS(
        covstruc = inputs$covstruc,
        SNPs = inputs$SNPs,
        estimation = "DWLS",
        parallel = FALSE,
        GC = "standard"
      )
    }
  )
})
Rprof(NULL)
elapsed <- proc.time()[["elapsed"]] - start

summary <- if (memory_profile) {
  summaryRprof(profile_path, memory = "both")
} else {
  summaryRprof(profile_path)
}
write.csv(summary$by.self, self_csv_path)
write.csv(summary$by.total, total_csv_path)

summary_lines <- c(
  paste("workflow:", workflow),
  paste("n_snp:", n_snp),
  paste("k:", k),
  paste("elapsed_sec:", elapsed),
  paste("checksum:", checksum_result(result)),
  paste("result_rows:", count_result_rows(result)),
  paste("result_cols:", count_result_cols(result)),
  paste("result_size_bytes:", as.numeric(object.size(result))),
  paste("stdout_lines:", length(output)),
  paste("memory_profile:", memory_profile),
  "",
  "Top by total:",
  capture.output(print(utils::head(summary$by.total, 30))),
  "",
  "Top by self:",
  capture.output(print(utils::head(summary$by.self, 30)))
)

writeLines(summary_lines, summary_path)
writeLines(summary_lines)
