\name{GenomicSEM-fast-options}
\alias{GenomicSEM-fast-options}
\title{Rust-backed fast path options}
\description{
Internal options controlling Rust-backed acceleration and high-throughput prep paths.
}
\details{
The Rust backend is used automatically for low-level matrix kernels when available and
\code{options(GenomicSEM.use_rust = TRUE)}. Additional model-fitting acceleration is opt-in:

\itemize{
  \item \code{GenomicSEM.fast_commonfactor_fit}: use the batched Rust DWLS solver for supported \code{commonfactorGWAS()} runs.
  \item \code{GenomicSEM.fast_usergwas_fit}: use the batched Rust RAM-matrix DWLS solver for supported \code{userGWAS()} runs.
  \item \code{GenomicSEM.fast_usergwas_native_setup}: compile and initialize supported \code{userGWAS()} fast-path models without the lavaan setup pass. Enabled by default when the native fast path is requested.
  \item \code{GenomicSEM.fast_table_read}: use the \code{data.table::fread()} prep-path parser for \code{munge()} and \code{sumstats()}, with automatic fallback to the legacy \code{read.table()} parser. Enabled by default.
  \item \code{GenomicSEM.fast_munge_engine}: use the native Rust streaming prep engine for supported \code{munge()} inputs, avoiding full GWAS table materialization in R. Enabled by default when the native library is available.
  \item \code{GenomicSEM.fast_munge_threads}: number of Rust worker threads for batched native \code{munge()} calls over multiple files. Defaults to one thread unless \code{parallel=TRUE}, in which case the existing \code{cores} argument is reused when supplied.
  \item \code{GenomicSEM.fast_munge_qc}: use the native Rust row-QC and Z-score transform kernel inside supported \code{munge()} runs. Enabled by default when the native library is available.
  \item \code{GenomicSEM.fast_sumstats_engine}: use the native Rust streaming prep engine for supported \code{sumstats()} inputs, avoiding full GWAS table materialization in R while preserving the existing R reference filtering and listwise merge. Enabled by default when the native library is available.
  \item \code{GenomicSEM.fast_sumstats_qc}: use the native Rust row-QC, allele flip, Z-score, and beta/SE transform kernel inside supported \code{sumstats()} runs. Enabled by default when the native library is available.
  \item \code{GenomicSEM.fast_ldsc_read}: use \code{data.table::fread()} for \code{ldsc()} and \code{s_ldsc()} chromosome, weight, M, and trait summary-stat table ingestion. Enabled by default.
  \item \code{GenomicSEM.fast_ldsc_blocks}: use the native Rust block crossproduct kernel inside \code{ldsc()} and \code{s_ldsc()} jackknife setup. Enabled by default when the native library is available.
  \item \code{GenomicSEM.fast_ldsc_threads}: number of Rust worker threads for the LDSC block crossproduct kernel. Defaults to up to 4 local cores.
  \item \code{GenomicSEM.fast_snp_join}: use the experimental fast SNP join helper in prep paths. Disabled by default because local benchmarks did not show a stable gain over the existing joins.
  \item \code{GenomicSEM.fast_diagnostics}: print whether a fast path was used or why a requested fast path fell back.
  \item \code{GenomicSEM.fast_strict}: error instead of falling back when a requested batched fast path is unavailable.
}

Supported batched paths currently require DWLS estimation and do not support MPI.
The \code{userGWAS()} batched path requires \code{smooth_check = FALSE}; supported
\code{sub=} extraction remains on the native path. Unsupported inputs use the existing
R/lavaan path unless \code{GenomicSEM.fast_strict = TRUE}.

The streaming \code{munge()} and \code{sumstats()} engines currently target simple
whitespace-delimited plain-text or gzip files with uniquely identified SNP, allele, effect,
P-value, and required N/SE columns. Unsupported prep features fall back to the existing R
implementation.

Returned objects include attributes \code{GenomicSEM.fast_path} and
\code{GenomicSEM.fast_threads}; when fallback occurs they may also include
\code{GenomicSEM.fast_fallback_reason}.
}
\examples{
\dontrun{
options(
  GenomicSEM.fast_commonfactor_fit = TRUE,
  GenomicSEM.fast_usergwas_fit = TRUE,
  GenomicSEM.fast_diagnostics = TRUE
)
}
}
