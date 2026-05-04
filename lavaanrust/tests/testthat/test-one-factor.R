test_that("one-factor DWLS slice matches frozen lavaan fixtures", {
  fixture <- one_factor_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  par_table <- lavaanrust::parTable_rust(fit)

  expect_s4_class(fit, "lavaan_rust_fit")
  expect_equal(par_table$est[1:3], fixture$loadings, tolerance = 2e-7)
  expect_equal(par_table$est[5:7], fixture$residuals, tolerance = 2e-7)
  expect_equal(lavaanrust::fitted_rust(fit)$cov, fixture$implied, tolerance = 2e-7)
  expect_equal(lavaanrust::lavInspect_rust(fit, "delta"), fixture$delta, tolerance = 2e-7)
  expect_true(lavaanrust::lavInspect_rust(fit, "converged"))
})

test_that("unsupported syntax falls back to lavaan", {
  fixture <- one_factor_fixture()
  model <- "F1 =~ A + B + C"
  fit <- suppressWarnings(lavaanrust::sem_rust(
    model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    se = "standard",
    sample.nobs = 2
  ))

  expect_s4_class(fit, "lavaan")
})
