#' Minimal compatibility fit object for the supported lavaanrust surface.
#'
#' @exportClass lavaan_rust_fit
methods::setClass(
  "lavaan_rust_fit",
  slots = c(
    ParTable = "data.frame",
    observed = "matrix",
    implied = "list",
    delta = "matrix",
    WLS.V = "matrix",
    fit = "numeric",
    cor.lv = "matrix",
    se = "list",
    converged = "logical",
    observed_names = "character",
    latent_name = "character",
    model_kind = "character"
  )
)
