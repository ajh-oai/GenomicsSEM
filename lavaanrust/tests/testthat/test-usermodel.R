test_that("marker-scaled one-factor usermodel slice matches frozen lavaan estimates", {
  fixture <- usermodel_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    std.lv = FALSE
  )

  expect_s4_class(fit, "lavaan_rust_fit")
  expect_equal(lavaanrust::parTable_rust(fit)$est, fixture$marker_est, tolerance = 2e-6)
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(6L, 6L))
  expect_equal(unname(lavaanrust::lavInspect_rust(fit, "cor.lv")), matrix(fixture$marker_est[[7L]]), tolerance = 2e-6)
})

test_that("simple one-factor usermodel slice also supports std.lv scaling", {
  fixture <- usermodel_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    std.lv = TRUE
  )

  expect_equal(lavaanrust::parTable_rust(fit)$est, fixture$std_lv_est, tolerance = 2e-6)
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(6L, 6L))
})

test_that("standardizedSolution_rust exposes the one-factor standardized rows usermodel needs", {
  fixture <- usermodel_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    std.lv = FALSE
  )
  standardized <- lavaanrust::standardizedSolution_rust(fit)

  expect_equal(standardized$est.std, fixture$standardized_est, tolerance = 2e-6)
  expect_equal(is.na(standardized$pvalue), c(rep(FALSE, 6L), TRUE))
})
