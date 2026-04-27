.LOG <- function(..., file, print = TRUE) {
  msg <- paste0(..., "\n")
  if (print) cat(msg)
  cat(msg, file = file, append = TRUE)
}

.get_renamed_colnames <- function(hold_names, userprovided, checkforsingle=c(), filename, N_provided, log.file,
                                  warnz=FALSE, warn_for_missing=c(), stop_on_missing=c(), utilfuncs=NULL) {
  interpreted_names <- list(
    SNP=c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID","PREDICTOR","SNP_ID", "VARIANTID", "VARIANT_ID", "RSIDS", "RS_ID"),
    A1=c("A1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF"),
    A2=c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT", "A0"),
    effect=c("OR","B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT","EST", "BETA1", "LOGOR"),
    INFO=c("INFO", "IMPINFO"),
    P=c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P"),
    N=c("N","WEIGHT","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES", "SAMPLESIZE", "NEFF", "N_EFF", "N_EFFECTIVE", "SUMNEFF"),
    MAF=c("MAF", "CEUAF", "FREQ1", "EAF", "FREQ1.HAPMAP", "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ", "FREQ.A1", "A1FREQ", "ALLELEFREQ","EFFECT_ALLELE_FREQUENCY"),
    Z=c("Z", "ZSCORE", "Z-SCORE", "ZSTATISTIC", "ZSTAT", "Z-STATISTIC"),
    SE=c("STDERR", "SE", "STDERRLOGOR", "SEBETA", "STANDARDERROR", "STANDARD_ERROR"),
    DIRECTION=c("DIRECTION", "DIREC", "DIRE", "SIGN")
  )
  full_names <- list(
    P="P-value",
    A1="effect allele",
    A2="other allele",
    effect="beta or effect",
    SNP="rs-id",
    SE="standard error",
    DIRECTION="direction"
  )
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
        assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  if (all(c("ALT", "REF") %in% hold_names)) {
    .LOG(paste0("Found REF and ALT columns in the summary statistic file ", filename, ". Please note that REF will be interpreted as A1 (effect allele) and ALT as A2 (other allele)"), print=TRUE, file=log.file)
  }
  if (N_provided) {
    interpreted_names[["N"]] <- NULL
  } else {
    if ("NEFF" %in% hold_names | "N_EFF" %in% hold_names | "N_EFFECTIVE" %in% hold_names | "SUMNEFF" %in% hold_names) {
      .LOG("Found an NEFF column for sample size. \n
Please note that this is likely effective sample size and should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.\n
Be aware that some NEFF columns reflect half of the effective sample size; the function will automatically double the column names if recognized [check above in .log file to determine if this is the case].
If the Neff value is halved in the summary stats, but not recognized by the munge function, this should be manually doubled prior to running munge.", file=log.file)
    }
  }
  for (col in names(interpreted_names)) {
    if (col %in% names(userprovided)) {
      .LOG("Interpreting the ",userprovided[[col]]," column as the ",col, " column, as requested",file=log.file)
      hold_names[ hold_names == toupper(userprovided[[col]]) ] <- col
    } else if (col %in% hold_names) {
      .LOG("Interpreting the ",col," column as the ",col, " column.",file=log.file)
    } else if (any(interpreted_names[[col]] %in% hold_names)) {
      .LOG("Interpreting the ", hold_names[ hold_names %in% interpreted_names[[col]] ], " column as the ",col," column.",file=log.file)
      hold_names[ hold_names %in% interpreted_names[[col]] ] <- col
    } else if ((col == "effect")){
      if (any(interpreted_names[["Z"]] %in% hold_names)) {
        if (!warnz) {
          .LOG("Interpreting the ", hold_names[hold_names %in% interpreted_names[["Z"]] ] , " column as the ",col," column.",file=log.file)
          hold_names[hold_names %in% interpreted_names[["Z"]] ] <- col
        } else {
          .LOG("There appears to be a Z-statistic column in the summary statistic file ", filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column.", print=FALSE, file=log.file)
          warning(paste0("There appears to be a Z-statistic column in the summary statistic file ", filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column."))
        }
      }
    } else {
      if (col %in% warn_for_missing) {
        .LOG('Cannot find ', col, ' column, try renaming it to ', col, ' in the summary statistics file for:',filename,file=log.file)
      } else if (col %in% stop_on_missing) {
        stop(paste0('Cannot find ', col, ' column, try renaming it to ', col, ' in the summary statistics file for:',filename))
      }
    }
  }
  # Print log and throw warning messages if multiple or no columns were found for those specified in checkforsingle
  if (length(checkforsingle) > 0) {
    for (col in checkforsingle) {
      if(sum(hold_names == col) == 0) {
        .LOG('Cannot find ',full_names[[col]],' column, try renaming it ', col, ' in the summary statistics file for:',filename,file=log.file)
        warning(paste0('Cannot find ',full_names[[col]],' column, try renaming it ', col, ' in the summary statistics file for:', filename))
      }
      if(sum(hold_names == col) > 1) {
        .LOG('Multiple columns are being interpreted as the ',full_names[[col]],' column, try renaming the column you dont want interpreted to ', col, '2 in the summary statistics file for:',filename,file=log.file)
        warning(paste0('Multiple columns are being interpreted as the ',full_names[[col]],' column, try renaming the column you dont want interpreted to ', col, '2 in the summary statistics file for:', filename))
      }
    }
  }
  return(hold_names)
}

#function to rearrange the sampling covariance matrix from original order to lavaan's order:
#'k' is the number of variables in the model
#'fit' is the fit function of the regression model
#'names' is a vector of variable names in the order you used
.rearrange <- function (k, fit, names) {
    order1 <- names
    order2 <- rownames(inspect(fit)[[1]]) #order of variables
    kst <- k*(k+1)/2
    covA <- matrix(NA, k, k)
    covA[lower.tri(covA, diag = TRUE)] <- 1:kst
    covA <- t(covA)
    covA[lower.tri(covA, diag = TRUE)] <- 1:kst
    colnames(covA) <- rownames(covA) <- order1 #give A actual variable order from lavaan output
    #reorder A by order2
    covA <- covA[order2, order2] #rearrange rows/columns
    vec2 <- lav_matrix_vech(covA) #grab new vectorized order
    return(vec2)
}

##modification of trycatch that allows the results of a failed run to still be saved
.tryCatch.W.E <- function(expr) {
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler), warning = W)
}

.get_V_full <- function(k, V_LD, varSNPSE2, V_SNP) {
    ##create shell of full sampling covariance matrix
    V_Full<-diag(((k+1)*(k+2))/2)

    ##input the ld-score regression region of sampling covariance from ld-score regression SEs
    V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD

    ##add in SE of SNP variance as first observation in sampling covariance matrix
    V_Full[1,1]<-varSNPSE2

    ##add in SNP region of sampling covariance matrix
    V_Full[2:(k+1),2:(k+1)]<-V_SNP
    return(V_Full)
}

.get_V_SNP <- function(SE_SNP, I_LD, varSNP, GC, coords, k, i) {
     V_SNP<-diag(k)
    #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
    if(GC == "conserv"){
      for (p in 1:nrow(coords)) {
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) {
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i])^2
        }
      }
    }

    if(GC == "standard"){
      for (p in 1:nrow(coords)) {
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) {
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*sqrt(I_LD[x,x])*varSNP[i])^2
        }
      }
    }

    if(GC == "none"){
      for (p in 1:nrow(coords)) {
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) {
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*varSNP[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*varSNP[i])^2
        }
      }
    }
    return(V_SNP)
}

         
.get_S_Full<-function(n_phenotypes,S_LD,varSNP,beta_SNP,TWAS,i){
  #create empty vector for S_SNP
  S_SNP <- vector(mode="numeric",length=n_phenotypes+1)

#enter SNP variance from reference panel as first observation
S_SNP[1] <- varSNP[i]

#enter SNP covariances (standardized beta * SNP variance from refference panel)
for (p in 1:n_phenotypes) {
  S_SNP[p+1] <- varSNP[i]*beta_SNP[i,p]
}

#create shell of the full S (observed covariance) matrix
S_Full <- diag(n_phenotypes+1)

##add the LD portion of the S matrix
S_Full[(2:(n_phenotypes+1)),(2:(n_phenotypes+1))] <- S_LD

##add in observed SNP variances as first row/column
S_Full[1:(n_phenotypes+1),1] <- S_SNP
S_Full[1,1:(n_phenotypes+1)] <- t(S_SNP)

##pull in variables names specified in LDSC function and name first column as SNP
if(TWAS){
  colnames(S_Full) <- c("Gene", colnames(S_LD))
} else {
  colnames(S_Full) <- c("SNP", colnames(S_LD))
}

##name rows like columns
rownames(S_Full) <- colnames(S_Full)

return(S_Full)
}


.get_Z_pre <- function(i, beta_SNP, SE_SNP, I_LD, GC) {
    if(GC == "conserv"){
        Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*diag(I_LD))
    }
    if(GC=="standard"){
        Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD)))
    }
    if(GC=="none"){
        Z_pre<-beta_SNP[i,]/SE_SNP[i,]
    }
    return(Z_pre)
}

.diag_inverse_from_values <- function(values, toler = .Machine$double.eps) {
  values <- as.numeric(values)
  if (!isTRUE(getOption("GenomicSEM.fast_diag_inverse", TRUE))) {
    fallback <- diag(length(values))
    diag(fallback) <- values
    return(solve(fallback, tol = toler))
  }

  scale <- max(abs(values))

  if (!is.finite(scale) || scale == 0 || any(!is.finite(values)) || min(abs(values)) / scale < toler) {
    fallback <- diag(length(values))
    diag(fallback) <- values
    return(solve(fallback, tol = toler))
  }

  out <- diag(length(values))
  diag(out) <- 1 / values
  out
}

.read_sumstats_table_fallback <- function(filename, p_as_character = FALSE) {
  if (p_as_character) {
    return(read.table(filename, header = TRUE, quote = "\"", fill = TRUE,
                      colClasses = c(P = "character"),
                      na.strings = c(".", "NA", "")))
  }

  read.table(filename, header = TRUE, quote = "\"", fill = TRUE,
             na.strings = c(".", "NA", ""))
}

.read_sumstats_table <- function(filename, p_as_character = FALSE) {
  if (!isTRUE(getOption("GenomicSEM.fast_table_read", TRUE))) {
    return(.read_sumstats_table_fallback(filename, p_as_character))
  }

  col_classes <- NULL
  if (p_as_character) {
    col_classes <- c(P = "character")
  }

  out <- tryCatch(
    suppressWarnings(fread(
      filename,
      header = TRUE,
      fill = TRUE,
      data.table = FALSE,
      check.names = TRUE,
      colClasses = col_classes,
      na.strings = c(".", "NA", ""),
      showProgress = FALSE
    )),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(.read_sumstats_table_fallback(filename, p_as_character))
  }

  if (p_as_character && "P" %in% names(out)) {
    out$P <- as.character(out$P)
  }

  out
}

.snp_inner_join_fallback <- function(x, y, by, mode) {
  if (identical(mode, "merge")) {
    return(merge(x, y, by = by, all.x = FALSE, all.y = FALSE))
  }

  suppressWarnings(inner_join(x, y, by = by))
}

.snp_inner_join <- function(x, y, by = "SNP", mode = "inner_join", sort_by_snp = FALSE) {
  if (!isTRUE(getOption("GenomicSEM.fast_snp_join", FALSE)) ||
      anyDuplicated(x[[by]]) || anyDuplicated(y[[by]])) {
    return(.snp_inner_join_fallback(x, y, by, mode))
  }

  y_match <- match(x[[by]], y[[by]])
  keep <- !is.na(y_match)
  x_out <- x[keep, , drop = FALSE]
  y_out <- y[y_match[keep], setdiff(names(y), by), drop = FALSE]

  common <- intersect(names(x_out), names(y_out))
  if (length(common) > 0L) {
    names(x_out)[match(common, names(x_out))] <- paste0(common, ".x")
    names(y_out)[match(common, names(y_out))] <- paste0(common, ".y")
  }

  out <- cbind.data.frame(x_out, y_out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  if (sort_by_snp && nrow(out) > 1L) {
    out <- out[order(out[[by]]), , drop = FALSE]
    rownames(out) <- NULL
  }

  out
}

.ldsc_selected_chromosomes <- function(chr, select) {
  if (isFALSE(select)) {
    return(seq_len(chr))
  }
  if (identical(select, "ODD")) {
    return(seq(1, chr, 2))
  }
  if (identical(select, "EVEN")) {
    return(seq(2, chr, 2))
  }
  if (is.numeric(select)) {
    return(select)
  }
  stop("select must be one of the following values: FALSE, 'ODD', 'EVEN', or a numeric vector of chromosome numbers.")
}

.ldsc_read_table <- function(path) {
  if (!isTRUE(getOption("GenomicSEM.fast_ldsc_read", TRUE))) {
    return(suppressMessages(read_delim(
      path,
      delim = "\t",
      escape_double = FALSE,
      trim_ws = TRUE,
      progress = FALSE
    )))
  }

  fread(
    path,
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )
}

.ldsc_read_m <- function(path) {
  if (!isTRUE(getOption("GenomicSEM.fast_ldsc_read", TRUE))) {
    return(suppressMessages(read_csv(path, col_names = FALSE)))
  }

  fread(
    path,
    header = FALSE,
    data.table = FALSE,
    showProgress = FALSE
  )
}

.ldsc_read_chromosome_tables <- function(root, suffix, chromosomes) {
  tables <- lapply(chromosomes, function(i) {
    .ldsc_read_table(file.path(root, paste0(i, suffix)))
  })
  as.data.frame(rbindlist(tables, use.names = TRUE, fill = TRUE))
}

.ldsc_read_m_files <- function(root, chromosomes) {
  tables <- lapply(chromosomes, function(i) {
    .ldsc_read_m(file.path(root, paste0(i, ".l2.M_5_50")))
  })
  as.data.frame(rbindlist(tables, use.names = TRUE, fill = TRUE))
}

.ldsc_read_file_list <- function(files) {
  tables <- lapply(files, .ldsc_read_table)
  as.data.frame(rbindlist(tables, use.names = TRUE, fill = TRUE))
}

.ldsc_read_m_file_list <- function(files) {
  tables <- lapply(files, .ldsc_read_m)
  as.data.frame(rbindlist(tables, use.names = TRUE, fill = TRUE))
}

.ldsc_block_products_r <- function(weighted.LD, weighted.chi, n.blocks) {
  weighted.LD <- as.matrix(weighted.LD)
  weighted.chi <- as.numeric(weighted.chi)

  n.snps <- nrow(weighted.LD)
  n.annot <- ncol(weighted.LD)
  select.from <- floor(seq(from = 1, to = n.snps, length.out = n.blocks + 1))
  select.from <- select.from[seq_len(n.blocks)]
  select.to <- if (n.blocks == 1L) n.snps else c(select.from[2:n.blocks] - 1, n.snps)

  xty.block.values <- matrix(NA_real_, nrow = n.blocks, ncol = n.annot)
  xtx.block.values <- matrix(NA_real_, nrow = n.annot * n.blocks, ncol = n.annot)
  for (i in seq_len(n.blocks)) {
    rows <- select.from[i]:select.to[i]
    out.rows <- ((i - 1L) * n.annot + 1L):(i * n.annot)
    weighted.LD.block <- weighted.LD[rows, , drop = FALSE]
    xty.block.values[i, ] <- crossprod(weighted.LD.block, weighted.chi[rows])
    xtx.block.values[out.rows, ] <- crossprod(weighted.LD.block)
  }
  colnames(xty.block.values) <- colnames(weighted.LD)
  colnames(xtx.block.values) <- colnames(weighted.LD)

  xty <- as.matrix(colSums(xty.block.values))
  xtx <- matrix(NA_real_, nrow = n.annot, ncol = n.annot)
  colnames(xtx) <- colnames(weighted.LD)
  for (i in seq_len(n.annot)) {
    xtx[i, ] <- colSums(xtx.block.values[seq(from = i, to = nrow(xtx.block.values), by = n.annot), , drop = FALSE])
  }

  list(
    xty.block.values = xty.block.values,
    xtx.block.values = xtx.block.values,
    xty = xty,
    xtx = xtx,
    delete.from = seq(from = 1, to = nrow(xtx.block.values), by = n.annot),
    delete.to = seq(from = n.annot, to = nrow(xtx.block.values), by = n.annot)
  )
}

.ldsc_block_products <- function(weighted.LD, weighted.chi, n.blocks) {
  if (!.genomicssem_use_rust() || !isTRUE(getOption("GenomicSEM.fast_ldsc_blocks", TRUE))) {
    return(.ldsc_block_products_r(weighted.LD, weighted.chi, n.blocks))
  }

  weighted.LD <- as.matrix(weighted.LD)
  weighted.chi <- as.numeric(weighted.chi)
  n.blocks <- as.integer(n.blocks)
  n.threads <- getOption("GenomicSEM.fast_ldsc_threads", NA_integer_)[1]
  if (is.null(n.threads) || is.na(n.threads)) {
    cores <- detectCores()
    if (is.na(cores)) {
      cores <- 1L
    }
    n.threads <- min(4L, max(1L, cores - 1L))
  }
  n.threads <- as.integer(n.threads)

  if (is.na(n.blocks) || n.blocks <= 0L || is.na(n.threads) || n.threads <= 0L ||
      nrow(weighted.LD) == 0L || ncol(weighted.LD) == 0L || n.blocks > nrow(weighted.LD)) {
    return(.ldsc_block_products_r(weighted.LD, weighted.chi, n.blocks))
  }

  out <- .Call(
    "genomicssem_ldsc_block_products_call",
    weighted.LD,
    weighted.chi,
    n.blocks,
    n.threads,
    PACKAGE = "GenomicSEM"
  )

  annot.names <- colnames(weighted.LD)
  if (!is.null(annot.names)) {
    colnames(out$xty.block.values) <- annot.names
    colnames(out$xtx.block.values) <- annot.names
    rownames(out$xty) <- annot.names
    colnames(out$xtx) <- annot.names
  }

  out
}

.get_V_full_r <- .get_V_full
.get_V_SNP_r <- .get_V_SNP
.get_S_Full_r <- .get_S_Full
.get_Z_pre_r <- .get_Z_pre

.genomicssem_native_state <- new.env(parent = emptyenv())
.genomicssem_native_state$rust_loaded <- NA

.genomicssem_use_rust <- function() {
  if (!isTRUE(getOption("GenomicSEM.use_rust", TRUE))) {
    return(FALSE)
  }

  if (is.na(.genomicssem_native_state$rust_loaded)) {
    .genomicssem_native_state$rust_loaded <- is.loaded("genomicssem_get_v_snp_call")
  }

  .genomicssem_native_state$rust_loaded
}

.fast_diagnostics_enabled <- function() {
  isTRUE(getOption("GenomicSEM.fast_diagnostics", FALSE))
}

.fast_strict_enabled <- function() {
  isTRUE(getOption("GenomicSEM.fast_strict", FALSE))
}

.fast_note <- function(workflow, message) {
  if (.fast_diagnostics_enabled()) {
    message(sprintf("GenomicSEM fast path [%s]: %s", workflow, message))
  }
}

.fast_fallback <- function(workflow, reason) {
  fallback_message <- sprintf("GenomicSEM fast path [%s] fallback: %s", workflow, reason)
  if (.fast_strict_enabled()) {
    stop(fallback_message, call. = FALSE)
  }
  if (.fast_diagnostics_enabled()) {
    message(fallback_message)
  }
  invisible(NULL)
}

.set_fast_path_attr <- function(x, path, threads = NA_integer_, reason = NULL) {
  attr(x, "GenomicSEM.fast_path") <- path
  attr(x, "GenomicSEM.fast_threads") <- threads
  if (!is.null(reason)) {
    attr(x, "GenomicSEM.fast_fallback_reason") <- reason
  }
  x
}

.allele_code <- function(x) {
  match(as.character(x), c("A", "C", "G", "T"))
}

.munge_qc_fast <- function(file, info.filter, maf.filter) {
  if (!.genomicssem_use_rust() || !isTRUE(getOption("GenomicSEM.fast_munge_qc", TRUE))) {
    return(NULL)
  }
  required <- c("A1.x", "A2.x", "A1.y", "A2.y", "effect", "P")
  if (!all(required %in% colnames(file))) {
    return(NULL)
  }

  p <- suppressWarnings(as.numeric(file$P))
  effect <- suppressWarnings(as.numeric(file$effect))
  if (any(is.finite(p) & (p < 0 | p > 1))) {
    return(NULL)
  }

  info <- if ("INFO" %in% colnames(file)) suppressWarnings(as.numeric(file$INFO)) else numeric(0)
  maf <- if ("MAF" %in% colnames(file)) suppressWarnings(as.numeric(as.character(file$MAF))) else numeric(0)

  out <- .Call(
    "genomicssem_munge_qc_call",
    as.integer(.allele_code(file$A1.x)),
    as.integer(.allele_code(file$A2.x)),
    as.integer(.allele_code(file$A1.y)),
    as.integer(.allele_code(file$A2.y)),
    as.numeric(effect),
    as.numeric(p),
    as.numeric(info),
    as.numeric(maf),
    as.numeric(info.filter),
    as.numeric(maf.filter),
    PACKAGE = "GenomicSEM"
  )

  n <- as.integer(out$n)
  if (n == 0L) {
    out$keep <- integer(0)
    out$z <- numeric(0)
  } else {
    out$keep <- out$keep[seq_len(n)]
    out$z <- out$z[seq_len(n)]
  }
  out
}

.sumstats_qc_fast <- function(file, info.filter, OLS, beta, linprob, se.logit) {
  if (!.genomicssem_use_rust() || !isTRUE(getOption("GenomicSEM.fast_sumstats_qc", TRUE))) {
    return(NULL)
  }
  required <- c("A1.x", "A2.x", "A1.y", "A2.y", "effect", "P")
  if (!all(required %in% colnames(file))) {
    return(NULL)
  }
  if (!("MAF.x" %in% colnames(file) || "MAF" %in% colnames(file))) {
    return(NULL)
  }
  if (!("N" %in% colnames(file)) && (OLS || linprob)) {
    return(NULL)
  }
  if (!("SE" %in% colnames(file)) && !(OLS && is.character(beta))) {
    return(NULL)
  }

  p <- suppressWarnings(as.numeric(file$P))
  if (any(is.finite(p) & (p < 0 | p > 1))) {
    return(NULL)
  }

  n <- if ("N" %in% colnames(file)) suppressWarnings(as.numeric(file$N)) else rep(NA_real_, nrow(file))
  se <- if ("SE" %in% colnames(file)) suppressWarnings(as.numeric(file$SE)) else rep(NA_real_, nrow(file))
  maf_ref <- if ("MAF.x" %in% colnames(file)) {
    suppressWarnings(as.numeric(file$MAF.x))
  } else {
    suppressWarnings(as.numeric(file$MAF))
  }
  maf_file <- if ("MAF.y" %in% colnames(file)) suppressWarnings(as.numeric(file$MAF.y)) else numeric(0)
  info <- if ("INFO" %in% colnames(file)) suppressWarnings(as.numeric(file$INFO)) else numeric(0)

  out <- .Call(
    "genomicssem_sumstats_qc_call",
    as.integer(.allele_code(file$A1.x)),
    as.integer(.allele_code(file$A2.x)),
    as.integer(.allele_code(file$A1.y)),
    as.integer(.allele_code(file$A2.y)),
    as.numeric(file$effect),
    as.numeric(se),
    as.numeric(p),
    as.numeric(n),
    as.numeric(maf_ref),
    as.numeric(maf_file),
    as.numeric(info),
    as.numeric(info.filter),
    as.logical(OLS),
    as.logical(is.character(beta)),
    as.logical(linprob),
    as.logical(se.logit),
    PACKAGE = "GenomicSEM"
  )

  out_n <- as.integer(out$n)
  if (out_n == 0L) {
    out$keep <- integer(0)
    out$beta <- numeric(0)
    out$se <- numeric(0)
  } else {
    rows <- seq_len(out_n)
    out$keep <- out$keep[rows]
    out$beta <- out$beta[rows]
    out$se <- out$se[rows]
  }
  out
}

.get_V_full <- function(k, V_LD, varSNPSE2, V_SNP) {
  if (.genomicssem_use_rust()) {
    return(.Call(
      "genomicssem_get_v_full_call",
      as.integer(k),
      as.matrix(V_LD),
      as.numeric(varSNPSE2),
      as.matrix(V_SNP),
      PACKAGE = "GenomicSEM"
    ))
  }

  .get_V_full_r(k, V_LD, varSNPSE2, V_SNP)
}

.get_V_SNP <- function(SE_SNP, I_LD, varSNP, GC, coords, k, i) {
  if (.genomicssem_use_rust()) {
    return(.Call(
      "genomicssem_get_v_snp_call",
      as.matrix(SE_SNP),
      as.matrix(I_LD),
      as.numeric(varSNP),
      as.character(GC),
      coords,
      as.integer(k),
      as.integer(i),
      PACKAGE = "GenomicSEM"
    ))
  }

  .get_V_SNP_r(SE_SNP, I_LD, varSNP, GC, coords, k, i)
}

.get_V_SNP_batch <- function(SE_SNP, I_LD, varSNP, GC, coords, k, n_threads = 1L) {
  if (!.genomicssem_use_rust()) {
    out <- array(NA_real_, dim = c(k, k, nrow(SE_SNP)))
    for (i in seq_len(nrow(SE_SNP))) {
      out[, , i] <- .get_V_SNP_r(SE_SNP, I_LD, varSNP, GC, coords, k, i)
    }
    return(out)
  }

  .Call(
    "genomicssem_get_v_snp_batch_call",
    as.matrix(SE_SNP),
    as.matrix(I_LD),
    as.numeric(varSNP),
    as.character(GC),
    coords,
    as.integer(k),
    as.integer(n_threads),
    PACKAGE = "GenomicSEM"
  )
}

.commonfactor_fast_start <- function(fit, k) {
  co <- coef(fit)
  start <- numeric(2 * k + 2)
  loading_idx <- grep("=~", names(co), fixed = TRUE)
  if (length(loading_idx) != k - 1L) {
    return(NULL)
  }
  start[seq_len(k - 1L)] <- co[loading_idx]

  factor_regression <- co[["F1~SNP"]]
  if (is.null(factor_regression) || is.na(factor_regression)) {
    return(NULL)
  }
  start[k] <- factor_regression

  for (j in seq_len(k)) {
    trait_name <- rownames(inspect(fit)[[1]])[j + 1L]
    residual_name <- paste0(trait_name, "~~", trait_name)
    residual <- co[[residual_name]]
    if (is.null(residual) || is.na(residual)) {
      return(NULL)
    }
    start[k + j] <- residual
  }

  factor_var <- co[["F1~~F1"]]
  snp_var <- co[["SNP~~SNP"]]
  if (is.null(factor_var) || is.null(snp_var) || is.na(factor_var) || is.na(snp_var)) {
    return(NULL)
  }
  start[2 * k + 1L] <- factor_var
  start[2 * k + 2L] <- snp_var
  start
}

.commonfactor_fast_start_from_cov <- function(S_LD) {
  S_LD <- as.matrix(S_LD)
  k <- ncol(S_LD)
  if (k < 1L || nrow(S_LD) != k) {
    return(NULL)
  }

  eig <- tryCatch(eigen(S_LD, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig) || !all(is.finite(eig$values)) || !all(is.finite(eig$vectors))) {
    return(NULL)
  }

  v <- eig$vectors[, 1L]
  if (abs(v[1L]) < sqrt(.Machine$double.eps)) {
    v[1L] <- ifelse(v[1L] < 0, -sqrt(.Machine$double.eps), sqrt(.Machine$double.eps))
  }

  lambdas <- v / v[1L]
  psi <- max(eig$values[1L] * v[1L]^2, sqrt(.Machine$double.eps))
  theta <- diag(S_LD) - lambdas^2 * psi
  theta[!is.finite(theta) | theta <= 0] <- sqrt(.Machine$double.eps)

  start <- numeric(2 * k + 2L)
  if (k > 1L) {
    start[seq_len(k - 1L)] <- lambdas[2:k]
  }
  start[k] <- 0
  start[k + seq_len(k)] <- theta
  start[2L * k + 1L] <- psi
  start[2L * k + 2L] <- 1
  start
}

.commonfactor_fit_fast <- function(S_Fullrun, V_Full_Reorder, W_diag, start, k, max_iter = 100L, tol = 1e-10) {
  if (!.genomicssem_use_rust() || is.null(start)) {
    return(NULL)
  }

  out <- tryCatch(
    .Call(
      "genomicssem_fit_commonfactor_main_call",
      as.integer(k),
      as.matrix(S_Fullrun),
      as.matrix(V_Full_Reorder),
      as.numeric(W_diag),
      as.numeric(start),
      as.integer(max_iter),
      as.numeric(tol),
      PACKAGE = "GenomicSEM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(NULL)
  }

  q <- 2 * k + 2L
  list(
    par = out[seq_len(q)],
    se = out[q + seq_len(q)],
    objective = out[2L * q + 1L],
    converged = isTRUE(out[2L * q + 2L] == 1),
    iterations = as.integer(out[2L * q + 3L])
  )
}

.commonfactor_q_fit_fast <- function(S_Fullrun, V_Full_Reorder, W_diag, fixed, start, k, max_iter = 100L, tol = 1e-10) {
  if (!.genomicssem_use_rust() || is.null(fixed) || is.null(start)) {
    return(NULL)
  }

  out <- tryCatch(
    .Call(
      "genomicssem_fit_commonfactor_q_call",
      as.integer(k),
      as.matrix(S_Fullrun),
      as.matrix(V_Full_Reorder),
      as.numeric(W_diag),
      as.numeric(fixed),
      as.numeric(start),
      as.integer(max_iter),
      as.numeric(tol),
      PACKAGE = "GenomicSEM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(NULL)
  }

  q <- 2 * k
  list(
    par = out[seq_len(q)],
    gamma_cov = matrix(out[q + seq_len(k * k)], nrow = k, ncol = k),
    objective = out[q + k * k + 1L],
    converged = isTRUE(out[q + k * k + 2L] == 1),
    iterations = as.integer(out[q + k * k + 3L])
  )
}

.commonfactor_batch_fit_fast <- function(S_LD, V_LD, I_LD, beta_SNP, SE_SNP, varSNP, GC, coords, varSNPSE2,
                                         start, k, n_threads = 1L, max_iter = 100L, max_iter_q = 500L,
                                         tol = 1e-10) {
  if (!.genomicssem_use_rust() || is.null(start)) {
    return(NULL)
  }

  out <- tryCatch(
    .Call(
      "genomicssem_fit_commonfactor_batch_call",
      as.integer(k),
      as.matrix(S_LD),
      as.matrix(V_LD),
      as.matrix(I_LD),
      as.matrix(beta_SNP),
      as.matrix(SE_SNP),
      as.numeric(varSNP),
      as.character(GC),
      coords,
      as.numeric(varSNPSE2),
      as.numeric(start),
      as.integer(max_iter),
      as.integer(max_iter_q),
      as.numeric(tol),
      as.integer(n_threads),
      PACKAGE = "GenomicSEM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(NULL)
  }

  out <- matrix(out, ncol = 7L, byrow = TRUE)
  colnames(out) <- c("est", "se_c", "Q", "main_converged", "q_converged", "main_iterations", "q_iterations")
  if (!all(out[, "main_converged"] == 1) || !all(out[, "q_converged"] == 1) ||
      !all(is.finite(out[, c("est", "se_c", "Q"), drop = FALSE]))) {
    return(NULL)
  }

  out
}

.sem_fast_numeric_value <- function(x, fallback = 0) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), x, fallback)
}

.sem_fast_compile <- function(ptable, observed_names) {
  unsupported_ops <- setdiff(unique(ptable$op), c("=~", "~", "~~"))
  if (length(unsupported_ops) > 0L) {
    return(list(supported = FALSE, reason = paste("unsupported ops:", paste(unsupported_ops, collapse = ","))))
  }

  observed_names <- as.character(observed_names)
  latent_names <- unique(ptable$lhs[ptable$op == "=~"])
  latent_names <- setdiff(latent_names, observed_names)
  total_names <- c(observed_names, latent_names)

  if (length(total_names) == 0L || anyNA(match(c(ptable$lhs, ptable$rhs), total_names))) {
    return(list(supported = FALSE, reason = "could not map model variables"))
  }

  free_ids <- sort(unique(ptable$free[ptable$free > 0]))
  if (length(free_ids) == 0L) {
    return(list(supported = FALSE, reason = "no free parameters"))
  }

  free_map <- seq_along(free_ids)
  names(free_map) <- as.character(free_ids)
  ptable$free_fast <- ifelse(ptable$free > 0, free_map[as.character(ptable$free)], 0L)

  total_n <- length(total_names)
  b_fixed <- matrix(0, total_n, total_n, dimnames = list(total_names, total_names))
  psi_fixed <- matrix(0, total_n, total_n, dimnames = list(total_names, total_names))
  b_free <- matrix(0L, total_n, total_n, dimnames = list(total_names, total_names))
  psi_free <- matrix(0L, total_n, total_n, dimnames = list(total_names, total_names))
  start <- rep(0, length(free_ids))

  value_for_row <- function(row) {
    est <- .sem_fast_numeric_value(ptable$est[row], NA_real_)
    if (is.finite(est)) {
      return(est)
    }
    ustart <- .sem_fast_numeric_value(ptable$ustart[row], NA_real_)
    if (is.finite(ustart)) {
      return(ustart)
    }
    if (ptable$op[row] == "~~" && ptable$lhs[row] == ptable$rhs[row]) {
      return(1)
    }
    0
  }

  for (row in seq_len(nrow(ptable))) {
    op <- ptable$op[row]
    lhs <- ptable$lhs[row]
    rhs <- ptable$rhs[row]
    value <- value_for_row(row)
    free <- as.integer(ptable$free_fast[row])

    if (free > 0L && start[free] == 0) {
      start[free] <- value
    }

    if (op == "=~") {
      target_row <- rhs
      target_col <- lhs
      if (free > 0L) {
        b_free[target_row, target_col] <- free
      } else {
        b_fixed[target_row, target_col] <- value
      }
    } else if (op == "~") {
      if (free > 0L) {
        b_free[lhs, rhs] <- free
      } else {
        b_fixed[lhs, rhs] <- value
      }
    } else if (op == "~~") {
      if (free > 0L) {
        psi_free[lhs, rhs] <- free
        psi_free[rhs, lhs] <- free
      } else {
        psi_fixed[lhs, rhs] <- value
        psi_fixed[rhs, lhs] <- value
      }
    }
  }

  start[!is.finite(start)] <- 0
  list(
    supported = TRUE,
    ptable = ptable,
    observed_names = observed_names,
    latent_names = latent_names,
    total_names = total_names,
    b_fixed = b_fixed,
    psi_fixed = psi_fixed,
    b_free = b_free,
    psi_free = psi_free,
    start = start
  )
}

.sem_fit_fast <- function(S_Fullrun, V_Full_Reorder, W_diag, spec, max_iter = 100L, tol = 1e-10) {
  if (!.genomicssem_use_rust() || is.null(spec) || !isTRUE(spec$supported)) {
    return(NULL)
  }

  S_Fullrun <- as.matrix(S_Fullrun)[spec$observed_names, spec$observed_names, drop = FALSE]

  out <- tryCatch(
    .Call(
      "genomicssem_fit_generic_sem_call",
      as.integer(length(spec$observed_names)),
      as.integer(length(spec$total_names)),
      as.matrix(S_Fullrun),
      as.matrix(V_Full_Reorder),
      as.numeric(W_diag),
      as.matrix(spec$b_fixed),
      as.matrix(spec$psi_fixed),
      matrix(as.integer(spec$b_free), nrow = nrow(spec$b_free), ncol = ncol(spec$b_free)),
      matrix(as.integer(spec$psi_free), nrow = nrow(spec$psi_free), ncol = ncol(spec$psi_free)),
      as.numeric(spec$start),
      as.integer(max_iter),
      as.numeric(tol),
      PACKAGE = "GenomicSEM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(NULL)
  }

  q <- length(spec$start)
  obs_n <- length(spec$observed_names)
  implied <- matrix(out[2L * q + seq_len(obs_n * obs_n)], nrow = obs_n, ncol = obs_n)
  dimnames(implied) <- list(spec$observed_names, spec$observed_names)

  list(
    par = out[seq_len(q)],
    se = out[q + seq_len(q)],
    implied = implied,
    objective = out[2L * q + obs_n * obs_n + 1L],
    converged = isTRUE(out[2L * q + obs_n * obs_n + 2L] == 1),
    iterations = as.integer(out[2L * q + obs_n * obs_n + 3L])
  )
}

.sem_fast_q_snp_info <- function(spec, model, S_LD, TWAS, Q_SNP) {
  if (!Q_SNP) {
    return(list(
      lv = character(),
      indices = matrix(integer(), nrow = 1L, ncol = 0L),
      lengths = integer(),
      df = integer()
    ))
  }

  lv <- spec$latent_names
  lines_SNP <- strsplit(model, "\n")[[1]]
  lines_SNP <- str_replace_all(lines_SNP, fixed(" "), "")
  if(TWAS){
    lines_SNP <- lines_SNP[grepl("Gene", lines_SNP)]
  }else{
    lines_SNP <- lines_SNP[grepl("SNP", lines_SNP)]
  }
  lv <- lv[lv %in% gsub(" ~.*|~.*", "", lines_SNP)]

  if (length(lv) == 0L) {
    return(list(
      lv = character(),
      indices = matrix(integer(), nrow = 1L, ncol = 0L),
      lengths = integer(),
      df = integer()
    ))
  }

  trait_names <- colnames(S_LD)
  indicators <- lapply(lv, function(factor_name) {
    subset(spec$ptable$rhs, spec$ptable$lhs == factor_name & spec$ptable$op == "=~")
  })
  df <- vapply(indicators, function(x) length(x) - 1L, integer(1))
  idx <- lapply(indicators, function(x) match(x, trait_names))
  lengths <- vapply(idx, function(x) if (length(x) > 0L && !is.na(x[[1L]])) length(x) else 0L, integer(1))
  max_len <- max(lengths, 1L)
  index_matrix <- matrix(0L, nrow = max_len, ncol = length(lv))
  for (j in seq_along(idx)) {
    if (lengths[[j]] > 0L && all(is.finite(idx[[j]]))) {
      index_matrix[seq_len(lengths[[j]]), j] <- as.integer(idx[[j]])
    }
  }

  list(lv = lv, indices = index_matrix, lengths = as.integer(lengths), df = as.integer(df))
}

.sem_fit_batch_fast <- function(S_LD, V_LD, I_LD, beta_SNP, SE_SNP, varSNP, GC, coords, varSNPSE2,
                                order, spec, observed_original_names, q_snp_info,
                                n_threads = 1L, max_iter = 100L, tol = 1e-10) {
  if (!.genomicssem_use_rust() || is.null(spec) || !isTRUE(spec$supported)) {
    return(NULL)
  }

  spec_to_original <- match(spec$observed_names, observed_original_names)
  if (anyNA(spec_to_original)) {
    return(NULL)
  }

  q <- length(spec$start)
  q_snp_n <- length(q_snp_info$lv)
  out <- tryCatch(
    .Call(
      "genomicssem_fit_generic_sem_batch_call",
      as.integer(length(spec$observed_names)),
      as.integer(length(spec$total_names)),
      as.matrix(S_LD),
      as.matrix(V_LD),
      as.matrix(I_LD),
      as.matrix(beta_SNP),
      as.matrix(SE_SNP),
      as.numeric(varSNP),
      as.character(GC),
      coords,
      as.numeric(varSNPSE2),
      matrix(as.integer(order), ncol = 1L),
      matrix(as.integer(spec_to_original), ncol = 1L),
      as.matrix(spec$b_fixed),
      as.matrix(spec$psi_fixed),
      matrix(as.integer(spec$b_free), nrow = nrow(spec$b_free), ncol = ncol(spec$b_free)),
      matrix(as.integer(spec$psi_free), nrow = nrow(spec$psi_free), ncol = ncol(spec$psi_free)),
      as.numeric(spec$start),
      matrix(as.integer(q_snp_info$indices), nrow = nrow(q_snp_info$indices), ncol = ncol(q_snp_info$indices)),
      matrix(as.integer(q_snp_info$lengths), ncol = 1L),
      as.integer(max_iter),
      as.numeric(tol),
      as.integer(n_threads),
      PACKAGE = "GenomicSEM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(NULL)
  }

  out_cols <- 2L * q + 1L + q_snp_n + 2L
  out <- matrix(out, ncol = out_cols, byrow = TRUE)
  converged_col <- 2L * q + 1L + q_snp_n + 1L
  if (!all(out[, converged_col] == 1) || !all(is.finite(out[, seq_len(2L * q + 1L), drop = FALSE]))) {
    return(NULL)
  }

  q_snp <- if (q_snp_n > 0L) {
    out[, (2L * q + 2L):(2L * q + 1L + q_snp_n), drop = FALSE]
  } else {
    matrix(numeric(), nrow = nrow(out), ncol = 0L)
  }
  colnames(q_snp) <- q_snp_info$lv

  list(
    par = out[, seq_len(q), drop = FALSE],
    se = out[, q + seq_len(q), drop = FALSE],
    chisq = out[, 2L * q + 1L],
    q_snp = q_snp,
    iterations = as.integer(out[, out_cols])
  )
}

.userGWAS_batch_results_fast <- function(batch_fit, spec, SNPs, TWAS, printwarn, Q_SNP, q_snp_info, df, npar, model) {
  f <- nrow(batch_fit$par)
  out <- vector(mode = "list", length = f)
  free_fast <- as.integer(spec$ptable$free_fast)
  free_rows <- free_fast > 0L
  warn_names <- if(printwarn) c("error","warning") else character()

  for (i in seq_len(f)) {
    Model_Output <- spec$ptable
    Model_Output$est[free_rows] <- batch_fit$par[i, free_fast[free_rows]]

    SE <- rep(NA_real_, nrow(Model_Output))
    SE[free_rows] <- batch_fit$se[i, free_fast[free_rows]]

    unstand <- subset(Model_Output, Model_Output$plabel != "" & Model_Output$free > 0)[,c("lhs","op","rhs","free","label","est")]
    unstand2 <- cbind(unstand, SE[Model_Output$plabel != "" & Model_Output$free > 0])
    colnames(unstand2)[7] <- "SE"

    other <- subset(Model_Output, (Model_Output$plabel == "" & Model_Output$op != ":=") | (Model_Output$free == 0 & Model_Output$plabel != ""))[,c("lhs","op","rhs","free","label","est")]
    other$SE <- rep(NA, nrow(other))

    if(nrow(other) > 0){
      final <- rbind(unstand2,other)
    }else{
      final <- unstand2
    }

    final$index <- as.numeric(row.names(final))
    final <- final[order(final$index), ]
    final$index <- NULL

    if(class(final$SE) != "factor"){
      final$Z_Estimate <- final$est/final$SE
      final$Pval_Estimate <- 2*pnorm(abs(final$Z_Estimate),lower.tail=FALSE)
    }else{
      final$SE <- as.character(final$SE)
      final$Z_Estimate <- NA
      final$Pval_Estimate <- NA
    }

    Q <- batch_fit$chisq[[i]]
    if(!(is.na(Q))){
      final$chisq <- rep(Q,nrow(final))
      final$chisq_df <- df
      final$chisq_pval <- pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
      final$AIC <- rep(Q + 2*npar,nrow(final))
    }else{
      final$chisq <- rep(NA, nrow(final))
      final$chisq_df <- rep(NA,nrow(final))
      final$chisq_pval <- rep(NA,nrow(final))
      final$AIC <- rep(NA, nrow(final))
    }

    if(Q_SNP){
      final$Q_SNP <- rep(NA,nrow(final))
      final$Q_SNP_df <- rep(NA,nrow(final))
      final$Q_SNP_pval <- rep(NA,nrow(final))
      if(length(q_snp_info$lv) > 0L){
        for(r in seq_len(nrow(final))){
          for(h in seq_along(q_snp_info$lv)){
            if(final$lhs[r] == q_snp_info$lv[h] & ((final$rhs[r] == "Gene" & TWAS) | (final$rhs[r] == "SNP" & !TWAS))) {
              final$Q_SNP[r] <- batch_fit$q_snp[i, h]
              final$Q_SNP_df[r] <- q_snp_info$df[h]
              final$Q_SNP_pval[r] <- pchisq(final$Q_SNP[r],final$Q_SNP_df[r],lower.tail=FALSE)
            }
          }
        }
      }
    }

    if(printwarn){
      final$error <- 0
      final$warning <- 0
    }

    final2 <- cbind(SNPs[i,],final,row.names=NULL)
    final2 <- subset(final2, final2$op != "da")
    final2$est <- ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)

    if(TWAS){
      if(Q_SNP){
        new_names <- c("Gene","Panel","HSQ", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","Q_SNP","Q_SNP_df","Q_SNP_pval", warn_names)
      }else{
        new_names <- c("Gene","Panel","HSQ", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC", warn_names)
      }
    }else{
      if(Q_SNP){
        new_names <- c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","Q_SNP","Q_SNP_df","Q_SNP_pval", warn_names)
      }else{
        new_names <- c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC", warn_names)
      }
    }
    colnames(final2) <- new_names
    out[[i]] <- final2
  }

  out
}

.get_S_Full <- function(n_phenotypes, S_LD, varSNP, beta_SNP, TWAS, i) {
  if (.genomicssem_use_rust()) {
    S_Full <- .Call(
      "genomicssem_get_s_full_call",
      as.integer(n_phenotypes),
      as.matrix(S_LD),
      as.numeric(varSNP),
      as.matrix(beta_SNP),
      as.integer(i),
      PACKAGE = "GenomicSEM"
    )

    first_name <- if (TWAS) "Gene" else "SNP"
    S_names <- c(first_name, colnames(S_LD))
    dimnames(S_Full) <- list(S_names, S_names)
    return(S_Full)
  }

  .get_S_Full_r(n_phenotypes, S_LD, varSNP, beta_SNP, TWAS, i)
}

.get_Z_pre <- function(i, beta_SNP, SE_SNP, I_LD, GC) {
  if (.genomicssem_use_rust()) {
    Z_pre <- .Call(
      "genomicssem_get_z_pre_call",
      as.integer(i),
      as.matrix(beta_SNP),
      as.matrix(SE_SNP),
      as.matrix(I_LD),
      as.character(GC),
      PACKAGE = "GenomicSEM"
    )

    names(Z_pre) <- colnames(as.matrix(beta_SNP))
    return(Z_pre)
  }

  .get_Z_pre_r(i, beta_SNP, SE_SNP, I_LD, GC)
}
