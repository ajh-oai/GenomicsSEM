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
  \item \code{GenomicSEM.fast_table_read}: use the \code{data.table::fread()} prep-path parser for \code{munge()} and \code{sumstats()}, with automatic fallback to the legacy \code{read.table()} parser. Enabled by default.
  \item \code{GenomicSEM.fast_snp_join}: use the experimental fast SNP join helper in prep paths. Disabled by default because local benchmarks did not show a stable gain over the existing joins.
  \item \code{GenomicSEM.fast_diagnostics}: print whether a fast path was used or why a requested fast path fell back.
  \item \code{GenomicSEM.fast_strict}: error instead of falling back when a requested batched fast path is unavailable.
}

Supported batched paths currently require DWLS estimation and do not support MPI.
The \code{userGWAS()} batched path also requires \code{sub = FALSE} and
\code{smooth_check = FALSE}. Unsupported inputs use the existing R/lavaan path unless
\code{GenomicSEM.fast_strict = TRUE}.

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
