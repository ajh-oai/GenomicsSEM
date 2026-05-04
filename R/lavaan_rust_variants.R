.with_lavaan_rust_backend <- function(fun) {
  if (!requireNamespace("lavaanrust", quietly = TRUE)) {
    stop(
      "The experimental lavaanrust package must be installed before using *_rust() wrappers.",
      call. = FALSE
    )
  }

  rust_env <- new.env(parent = environment(fun))
  rust_env$sem <- lavaanrust::sem_rust
  rust_env$lavaan <- lavaanrust::lavaan_rust
  rust_env$lavInspect <- lavaanrust::lavInspect_rust
  rust_env$inspect <- lavaanrust::inspect_rust
  rust_env$parTable <- lavaanrust::parTable_rust
  rust_env$fitted <- lavaanrust::fitted_rust
  rust_env$resid <- lavaanrust::resid_rust
  rust_env$lav_model_get_parameters <- lavaanrust::lav_model_get_parameters_rust
  rust_env$lav_func_jacobian_complex <- lavaanrust::lav_func_jacobian_complex_rust

  environment(fun) <- rust_env
  fun
}

# This reuses the exact original function body while rebinding only the lavaan
# surface underneath it.
commonfactor_rust <- .with_lavaan_rust_backend(commonfactor)
