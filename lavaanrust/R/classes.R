methods::setClass(
  "lavaan_rust_fit",
  slots = c(
    ParTable = "data.frame",
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
