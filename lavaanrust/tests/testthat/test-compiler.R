test_that("lavaan_fast compiler reproduces unrestricted userGWAS implied covariance", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    lavaanrust::parTable_rust(fit),
    colnames(fixture$sample_cov)
  )

  expect_s3_class(compiled, "lavaan_fast_compiled")
  expect_equal(compiled$variable_names, c("SNP", "A", "B", "C", "F1"))
  expect_equal(
    lavaanrust:::.lavaan_fast_implied_covariance(compiled),
    lavaanrust::fitted_rust(fit)$cov,
    tolerance = 1e-10
  )
  expect_equal(
    lavaanrust:::.lavaan_fast_implied_jacobian(compiled),
    lavaanrust::lavInspect_rust(fit, "delta"),
    tolerance = 1e-10
  )
  rust_surfaces <- lavaanrust:::.lavaan_fast_implied_surfaces_rust(compiled)
  expect_equal(rust_surfaces$implied, lavaanrust::fitted_rust(fit)$cov, tolerance = 1e-10)
  expect_equal(rust_surfaces$delta, lavaanrust::lavInspect_rust(fit, "delta"), tolerance = 1e-10)
})

test_that("lavaan_fast compiler preserves fixed-measurement userGWAS structure", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$fixed_model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    lavaanrust::parTable_rust(fit),
    colnames(fixture$sample_cov)
  )

  expect_equal(
    lavaanrust:::.lavaan_fast_implied_covariance(compiled),
    lavaanrust::fitted_rust(fit)$cov,
    tolerance = 1e-10
  )
  expect_equal(
    lavaanrust:::.lavaan_fast_implied_jacobian(compiled),
    lavaanrust::lavInspect_rust(fit, "delta"),
    tolerance = 1e-10
  )
  rust_surfaces <- lavaanrust:::.lavaan_fast_implied_surfaces_rust(compiled)
  expect_equal(rust_surfaces$implied, lavaanrust::fitted_rust(fit)$cov, tolerance = 1e-10)
  expect_equal(rust_surfaces$delta, lavaanrust::lavInspect_rust(fit, "delta"), tolerance = 1e-10)
  expect_equal(compiled$free_ids, seq_len(6L))
})

test_that("lavaan_fast compiler rejects syntax outside the first RAM subset", {
  fixture <- user_gwas_fixture()
  par_table <- fixture$fixed_model
  par_table <- rbind(
    par_table,
    transform(
      par_table[1L, , drop = FALSE],
      lhs = "ghost",
      op = ":=",
      rhs = "F1~SNP",
      free = 0L
    )
  )

  expect_error(
    lavaanrust:::.lavaan_fast_compile_par_table(par_table, colnames(fixture$sample_cov)),
    "does not yet support operators"
  )
})
