args <- commandArgs(trailingOnly = TRUE)
n_snp <- if (length(args) >= 1) as.integer(args[[1]]) else 50L
k <- if (length(args) >= 2) as.integer(args[[2]]) else 12L
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 1L

library(GenomicSEM)
source("benches/synthetic_inputs.R")

ns <- asNamespace("GenomicSEM")
.tryCatch.W.E <- get(".tryCatch.W.E", envir = ns)
.get_V_SNP <- get(".get_V_SNP", envir = ns)
.get_V_full <- get(".get_V_full", envir = ns)
.diag_inverse_from_values <- get(".diag_inverse_from_values", envir = ns)
.rearrange <- get(".rearrange", envir = ns)
.lavaan_resid <- get("resid", envir = ns)

inputs <- make_synthetic_genomicsem_inputs(n_snp = n_snp, k = k, seed = seed)
options(GenomicSEM.use_rust = TRUE)
options(GenomicSEM.fast_diag_inverse = TRUE)
options(GenomicSEM.fast_diagonal_wls = FALSE)

time_phase <- function(name, code) {
  start <- proc.time()[["elapsed"]]
  value <- force(code)
  elapsed <- proc.time()[["elapsed"]] - start
  phase_times[[name]] <<- phase_times[[name]] + elapsed
  value
}

SNPs <- data.frame(inputs$SNPs)
SNPs$A1 <- as.character(SNPs$A1)
SNPs$A2 <- as.character(SNPs$A2)
SNPs$SNP <- as.character(SNPs$SNP)
varSNP <- 2 * SNPs$MAF * (1 - SNPs$MAF)
varSNPSE2 <- (.0005)^2

V_LD <- as.matrix(inputs$covstruc[[1]])
S_LD <- as.matrix(inputs$covstruc[[2]])
I_LD <- as.matrix(inputs$covstruc[[3]])
diag(I_LD) <- ifelse(diag(I_LD) <= 1, 1, diag(I_LD))

beta_SNP <- SNPs[, grep("beta.", fixed = TRUE, colnames(SNPs))]
SE_SNP <- SNPs[, grep("se.", fixed = TRUE, colnames(SNPs))]
traits <- colnames(S_LD)
coords <- which(I_LD != "NA", arr.ind = TRUE)

model <- paste(
  paste0("F1 =~ ", paste(traits, collapse = " + ")),
  "F1 ~ SNP",
  paste0(traits, " ~ 0*SNP", collapse = "\n"),
  sep = "\n"
)

V_SNP <- .get_V_SNP(SE_SNP, I_LD, varSNP, "conserv", coords, k, 1)
V_Full <- .get_V_full(k, V_LD, varSNPSE2, V_SNP)
if (eigen(V_Full)$values[nrow(V_Full)] <= 0) {
  V_Full <- as.matrix(Matrix::nearPD(V_Full, corr = FALSE)$mat)
}
W <- solve(V_Full, tol = .Machine$double.eps)

S_SNP <- numeric(k + 1)
S_SNP[1] <- varSNP[1]
S_SNP[-1] <- varSNP[1] * as.numeric(beta_SNP[1, ])
S_Fullrun <- diag(k + 1)
S_Fullrun[2:(k + 1), 2:(k + 1)] <- S_LD
S_Fullrun[1:(k + 1), 1] <- S_SNP
S_Fullrun[1, 1:(k + 1)] <- S_SNP
dimnames(S_Fullrun) <- list(c("SNP", traits), c("SNP", traits))

reorder_model <- lavaan::sem(
  model,
  sample.cov = S_Fullrun,
  estimator = "DWLS",
  se = "standard",
  WLS.V = W,
  sample.nobs = 2,
  optim.dx.tol = .01,
  optim.force.converged = TRUE,
  control = list(iter.max = 1)
)
order <- .rearrange(k = k + 1, fit = reorder_model, names = rownames(S_Fullrun))

LavModel1 <- get(".commonfactorGWAS_main", envir = ns)(
  1,
  cores = 1,
  n = 1,
  S_LD,
  V_LD,
  I_LD,
  beta_SNP,
  SE_SNP,
  varSNP,
  varSNPSE2,
  "standard",
  coords,
  k,
  smooth_check = FALSE,
  model,
  toler = .Machine$double.eps,
  estimation = "DWLS",
  order,
  returnlavmodel = TRUE
)

phase_times <- c(
  build_v = 0,
  build_s = 0,
  main_lavaan = 0,
  main_sandwich_and_q_setup = 0,
  q_lavaan = 0,
  q_sandwich = 0
)

for (i in seq_len(n_snp)) {
  built_v <- time_phase("build_v", {
    V_SNP <- .get_V_SNP(SE_SNP, I_LD, varSNP, "standard", coords, k, i)
    V_Full <- .get_V_full(k, V_LD, varSNPSE2, V_SNP)
    if (eigen(V_Full)$values[nrow(V_Full)] <= 0) {
      V_Full <- as.matrix(Matrix::nearPD(V_Full, corr = FALSE)$mat)
    }
    V_Full_Reorder <- V_Full[order, order]
    W <- .diag_inverse_from_values(diag(V_Full_Reorder), toler = .Machine$double.eps)
    list(V_SNP = V_SNP, V_Full = V_Full, V_Full_Reorder = V_Full_Reorder, W = W)
  })

  S_Fullrun <- time_phase("build_s", {
    S_SNP <- numeric(k + 1)
    S_SNP[1] <- varSNP[i]
    S_SNP[-1] <- varSNP[i] * as.numeric(beta_SNP[i, ])
    S_Fullrun <- diag(k + 1)
    S_Fullrun[2:(k + 1), 2:(k + 1)] <- S_LD
    S_Fullrun[1:(k + 1), 1] <- S_SNP
    S_Fullrun[1, 1:(k + 1)] <- S_SNP
    dimnames(S_Fullrun) <- list(c("SNP", traits), c("SNP", traits))
    if (eigen(S_Fullrun)$values[nrow(S_Fullrun)] <= 0) {
      S_Fullrun <- as.matrix(Matrix::nearPD(S_Fullrun, corr = FALSE)$mat)
    }
    S_Fullrun
  })

  test <- time_phase("main_lavaan", {
    .tryCatch.W.E(lavaan::lavaan(
      sample.cov = S_Fullrun,
      WLS.V = built_v$W,
      ordered = NULL,
      sampling.weights = NULL,
      se = "standard",
      sample.mean = NULL,
      sample.th = NULL,
      sample.nobs = 2,
      group = NULL,
      cluster = NULL,
      constraints = "",
      NACOV = NULL,
      slotOptions = LavModel1@Options,
      slotParTable = LavModel1@ParTable,
      slotSampleStats = NULL,
      slotData = LavModel1@Data,
      slotModel = LavModel1@Model,
      slotCache = NULL,
      sloth1 = NULL
    ))
  })

  if (class(test$value)[1] != "lavaan" || grepl("solution has NOT", as.character(test$warning))) {
    next
  }
  Model1_Results <- test$value

  q_inputs <- time_phase("main_sandwich_and_q_setup", {
    S2.delt <- lavaan::lavInspect(Model1_Results, "delta")
    S2.W <- lavaan::lavInspect(Model1_Results, "WLS.V")
    lettuce <- S2.W %*% S2.delt
    bread <- solve(t(S2.delt) %*% lettuce, tol = .Machine$double.eps)
    Ohtt <- bread %*% t(lettuce) %*% built_v$V_Full_Reorder %*% lettuce %*% bread

    ModelQ <- lavaan::parTable(Model1_Results)
    ModelQ <- ModelQ[1:((k * 3) + 3), ]
    ModelQ$free <- c(rep(0, k + 1), 1:(k * 2), 0, 0)
    ModelQ$ustart <- ModelQ$est
    SNPresid <- .lavaan_resid(Model1_Results)$cov[k + 1, 1:k]
    for (t in seq_len(nrow(ModelQ))) {
      if (ModelQ$free[t] > 0 && ModelQ$free[t] <= k) {
        ModelQ$ustart[t] <- SNPresid[ModelQ$free[t]]
      }
    }
    ModelQ
  })

  testQ <- time_phase("q_lavaan", {
    .tryCatch.W.E(lavaan::sem(
      model = q_inputs,
      sample.cov = S_Fullrun,
      estimator = "DWLS",
      se = "standard",
      WLS.V = built_v$W,
      sample.nobs = 2,
      optim.dx.tol = .01
    ))
  })

  if (class(testQ$value)[1] != "lavaan") {
    next
  }
  ModelQ_Results <- testQ$value

  time_phase("q_sandwich", {
    S2.delt_Q <- lavaan::lavInspect(ModelQ_Results, "delta")
    S2.W_Q <- lavaan::lavInspect(ModelQ_Results, "WLS.V")
    lettuce_Q <- S2.W_Q %*% S2.delt_Q
    bread_Q <- solve(t(S2.delt_Q) %*% lettuce_Q, tol = .Machine$double.eps)
    Ohtt_Q <- bread_Q %*% t(lettuce_Q) %*% built_v$V_Full_Reorder %*% lettuce_Q %*% bread_Q
    V_eta <- Ohtt_Q[1:k, 1:k]
    eig <- eigen(V_eta)
    eta <- cbind(lavaan::inspect(ModelQ_Results, "list")[(k + 2):(2 * k + 1), 14])
    t(eta) %*% eig$vectors %*% .diag_inverse_from_values(eig$values) %*% t(eig$vectors) %*% eta
  })
}

phase_df <- data.frame(
  phase = names(phase_times),
  elapsed_sec = as.numeric(phase_times),
  pct = as.numeric(phase_times) / sum(phase_times),
  stringsAsFactors = FALSE
)
phase_df <- phase_df[order(-phase_df$elapsed_sec), ]
print(phase_df, row.names = FALSE)
