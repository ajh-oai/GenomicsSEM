\name{GenomicSEM-fast-options}
\alias{GenomicSEM-fast-options}
\title{Rust-backed fast path options}
\description{
Internal options controlling Rust-backed acceleration for GWAS workflows.
}
\details{
The Rust backend is used automatically for low-level matrix kernels when available and
\code{options(GenomicSEM.use_rust = TRUE)}. Additional model-fitting acceleration is opt-in:

\itemize{
  \item \code{GenomicSEM.fast_commonfactor_fit}: use the batched Rust DWLS solver for supported \code{commonfactorGWAS()} runs.
  \item \code{GenomicSEM.fast_usergwas_fit}: use the batched Rust RAM-matrix DWLS solver for supported \code{userGWAS()} runs.
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
