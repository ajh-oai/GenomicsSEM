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
  expect_equal(compiled$free_labels, lavaanrust:::.lavaan_fast_free_labels(compiled))
  expect_equal(compiled$stat_names, lavaanrust:::.stat_names(compiled$observed_names))
  free_row_counts <- lengths(lapply(seq_along(compiled$free_ids), function(free_position) {
    which(compiled$free_index == free_position)
  }))
  expect_equal(compiled$free_row_offsets, c(0L, cumsum(free_row_counts)))
  expect_equal(
    compiled$free_row_indices,
    as.integer(unlist(lapply(seq_along(compiled$free_ids), function(free_position) {
      which(compiled$free_index == free_position)
    }), use.names = FALSE))
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
  rust_surfaces_flat <- lavaanrust:::.lavaan_fast_implied_surfaces_rust_flat(compiled)
  expect_equal(rust_surfaces$implied, lavaanrust::fitted_rust(fit)$cov, tolerance = 1e-10)
  expect_equal(rust_surfaces$delta, lavaanrust::lavInspect_rust(fit, "delta"), tolerance = 1e-10)
  expect_equal(
    matrix(rust_surfaces_flat$implied, nrow = compiled$n_observed),
    unname(rust_surfaces$implied),
    tolerance = 1e-10
  )
  expect_equal(
    matrix(rust_surfaces_flat$delta, nrow = compiled$n_stats),
    unname(rust_surfaces$delta),
    tolerance = 1e-10
  )
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

test_that("lavaan_fast string parser handles modifiers and repeated terms", {
  fixture <- usermodel_fixture()
  model <- paste(
    "F1 =~ NA*A + start(.1)*A + start(1.1)*l2*B + l3*C",
    "F1 ~~ 1*F1",
    "A ~~ rvA*A",
    "B ~~ rvB*B",
    "C ~~ rvC*C",
    "rvA > .001",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, fixture$sample_cov)

  expect_equal(
    paste0(parsed$lhs, parsed$op, parsed$rhs),
    c("F1=~A", "F1=~B", "F1=~C", "F1~~F1", "A~~A", "B~~B", "C~~C")
  )
  expect_equal(parsed$free, c(1L, 2L, 3L, 0L, 4L, 5L, 6L))
  expect_equal(parsed$label, c("", "l2", "l3", "", "rvA", "rvB", "rvC"))
  expect_equal(parsed$ustart, c(0.1, 1.1, NA, 1, NA, NA, NA))
  expect_equal(parsed$lower, c(NA, NA, NA, NA, 0.001, NA, NA))
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(parsed, colnames(fixture$sample_cov))
  expect_equal(compiled$free_lower_bounds, c(-Inf, -Inf, -Inf, 0.001, 1e-10, 1e-10))
})

test_that("lavaan_fast generic string path fits labeled direct-effect RAM models", {
  fixture <- user_gwas_fixture()
  model <- paste(
    "F1 =~ 1*A + start(1.14285707631494)*l2*B + start(0.888888778086765)*l3*C",
    "F1 ~ start(0.094108558629525)*gamma*SNP",
    "A ~ start(0.02)*direct*SNP",
    "B ~ direct*SNP",
    "A ~~ start(0.60624995424363)*rvA*A",
    "B ~~ start(0.585714285882577)*rvB*B",
    "C ~~ start(0.638888930371178)*rvC*C",
    "F1 ~~ start(0.390030348982832)*psi*F1",
    "SNP ~~ 0.42*SNP",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, fixture$sample_cov)
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    parsed,
    colnames(fixture$sample_cov)
  )
  sample_cov <- lavaanrust:::.lavaan_fast_implied_covariance(compiled)
  fit <- lavaanrust::sem_rust(
    model,
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )

  expect_equal(fit@Model$model_kind, "ram_dwls_generic")
  expect_equal(lavaanrust::fitted_rust(fit)$cov, sample_cov, tolerance = 1e-8)
  expect_equal(
    lavaanrust::parTable_rust(fit)$free[match(c("A~SNP", "B~SNP"), paste0(
      lavaanrust::parTable_rust(fit)$lhs,
      lavaanrust::parTable_rust(fit)$op,
      lavaanrust::parTable_rust(fit)$rhs
    ))],
    rep(lavaanrust::parTable_rust(fit)$free[which(lavaanrust::parTable_rust(fit)$label == "direct")[[1L]]], 2L)
  )
})

test_that("lavaan_fast generic string path enforces simple lower bounds", {
  sample_cov <- matrix(1, nrow = 1L, dimnames = list("A", "A"))
  fit <- lavaanrust::sem_rust(
    paste("A ~~ rvA*A", "rvA > 1.5", sep = "\n"),
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = matrix(1, nrow = 1L)
  )

  expect_equal(lavaanrust::parTable_rust(fit)$est, 1.5, tolerance = 1e-10)
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
