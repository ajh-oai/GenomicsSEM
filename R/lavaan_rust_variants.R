.with_lavaan_rust_backend <- function(fun, helper_names = character()) {
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
  rust_env$class <- function(x) {
    if (methods::is(x, "lavaan_rust_fit")) {
      return(c("lavaan", "lavaan_rust_fit"))
    }

    base::class(x)
  }

  for (helper_name in helper_names) {
    helper <- get(helper_name, envir = environment(fun), inherits = TRUE)
    environment(helper) <- rust_env
    assign(helper_name, helper, envir = rust_env)
  }

  environment(fun) <- rust_env
  fun
}

# This reuses the exact original function body while rebinding only the lavaan
# surface underneath it.
commonfactor_rust <- .with_lavaan_rust_backend(commonfactor)
commonfactorGWAS_rust <- function(...) {
  dots <- list(...)
  parallel <- if ("parallel" %in% names(dots)) dots$parallel else formals(commonfactorGWAS)$parallel

  if (isTRUE(parallel)) {
    stop(
      "commonfactorGWAS_rust() currently supports parallel = FALSE only; the rust wrapper does not fall back to the original parallel worker path.",
      call. = FALSE
    )
  }

  rust_fun <- .with_lavaan_rust_backend(
    commonfactorGWAS,
    helper_names = c(".commonfactorGWAS_main", ".rearrange")
  )

  rust_fun(...)
}
