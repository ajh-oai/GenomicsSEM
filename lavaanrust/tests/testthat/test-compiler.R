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

test_that("lavaan_fast generic DWLS optimizer reproduces fixed-measurement userGWAS fit", {
  fixture <- user_gwas_fixture()
  specialized <- lavaanrust::sem_rust(
    fixture$fixed_model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  par_table <- fixture$fixed_model
  free_rows <- which(par_table$free > 0L)
  par_table$free[] <- 0L
  par_table$free[free_rows] <- seq_along(free_rows)
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    par_table,
    colnames(fixture$sample_cov)
  )
  generic <- lavaanrust:::.lavaan_fast_fit_dwls_rust(
    compiled,
    fixture$sample_cov,
    fixture$wls_v
  )

  expect_equal(generic$estimates, fixture$fixed_est[free_rows], tolerance = 2e-6)
  expect_equal(generic$implied, lavaanrust::fitted_rust(specialized)$cov, tolerance = 1e-8)
  expect_equal(generic$delta, lavaanrust::lavInspect_rust(specialized, "delta"), tolerance = 1e-8)
})

test_that("lavaan_fast generic parameter-table path supports direct SNP effect", {
  fixture <- user_gwas_fixture()
  par_table <- fixture$fixed_model
  free_rows <- which(par_table$free > 0L)
  par_table$free[] <- 0L
  par_table$free[free_rows] <- seq_along(free_rows)
  extra <- transform(
    par_table[1L, , drop = FALSE],
    id = max(par_table$id) + 1L,
    lhs = "A",
    op = "~",
    rhs = "SNP",
    user = 1L,
    free = max(par_table$free) + 1L,
    ustart = NA_real_,
    plabel = paste0(".p", max(par_table$id) + 1L, "."),
    start = 0.02,
    est = 0.02,
    se = 0
  )
  par_table <- rbind(par_table, extra)
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    par_table,
    colnames(fixture$sample_cov)
  )
  sample_cov <- lavaanrust:::.lavaan_fast_implied_covariance(compiled)
  fit <- lavaanrust::sem_rust(
    par_table,
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  refit <- lavaanrust::lavaan_rust(
    sample.cov = sample_cov,
    WLS.V = fixture$wls_v,
    slotOptions = fit@Options,
    slotParTable = fit@ParTable,
    slotData = fit@Data,
    slotModel = fit@Model
  )

  expect_equal(fit@Model$model_kind, "ram_dwls_generic")
  expect_equal(lavaanrust::fitted_rust(fit)$cov, sample_cov, tolerance = 1e-8)
  expect_equal(lavaanrust::parTable_rust(fit)$est, par_table$est, tolerance = 1e-8)
  expect_equal(lavaanrust::fitted_rust(refit)$cov, sample_cov, tolerance = 1e-8)
})

test_that("lavaan_fast native Jacobian sums shared directed-edge parameters", {
  fixture <- user_gwas_fixture()
  par_table <- fixture$fixed_model
  free_rows <- which(par_table$free > 0L)
  par_table$free[] <- 0L
  par_table$free[free_rows] <- seq_along(free_rows)
  shared_free <- max(par_table$free) + 1L
  direct_rows <- rbind(
    transform(
      par_table[1L, , drop = FALSE],
      id = max(par_table$id) + 1L,
      lhs = "A",
      op = "~",
      rhs = "SNP",
      user = 1L,
      free = shared_free,
      ustart = NA_real_,
      plabel = paste0(".p", max(par_table$id) + 1L, "."),
      start = 0.02,
      est = 0.02,
      se = 0
    ),
    transform(
      par_table[1L, , drop = FALSE],
      id = max(par_table$id) + 2L,
      lhs = "B",
      op = "~",
      rhs = "SNP",
      user = 1L,
      free = shared_free,
      ustart = NA_real_,
      plabel = paste0(".p", max(par_table$id) + 2L, "."),
      start = 0.02,
      est = 0.02,
      se = 0
    )
  )
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    rbind(par_table, direct_rows),
    colnames(fixture$sample_cov)
  )
  rust_surfaces <- lavaanrust:::.lavaan_fast_implied_surfaces_rust(compiled)

  expect_equal(
    rust_surfaces$delta,
    lavaanrust:::.lavaan_fast_implied_jacobian(compiled),
    tolerance = 1e-10
  )
})
