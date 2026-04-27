args <- commandArgs(trailingOnly = TRUE)
n_snp <- if (length(args) >= 1) as.integer(args[[1]]) else 100000L
n_traits <- if (length(args) >= 2) as.integer(args[[2]]) else 2L
n_annot <- if (length(args) >= 3) as.integer(args[[3]]) else 1L
n_blocks <- if (length(args) >= 4) as.integer(args[[4]]) else 200L
seed <- if (length(args) >= 5) as.integer(args[[5]]) else 1L

if (requireNamespace("GenomicSEM", quietly = TRUE)) {
  library(GenomicSEM)
} else {
  stop("GenomicSEM must be installed before running this benchmark", call. = FALSE)
}

genomicssem_ns <- asNamespace("GenomicSEM")
.ldsc_block_products <- get(".ldsc_block_products", envir = genomicssem_ns)

with_temp_cwd <- function(code) {
  old <- getwd()
  tmp <- tempfile("genomicssem-prep-bench-")
  dir.create(tmp)
  on.exit({
    setwd(old)
    unlink(tmp, recursive = TRUE)
  })
  setwd(tmp)
  force(code)
}

make_ref <- function(n) {
  alleles <- c("A", "C", "G", "T")
  a1 <- sample(alleles, n, replace = TRUE)
  a2 <- sample(alleles, n, replace = TRUE)
  same <- a1 == a2
  a2[same] <- alleles[(match(a1[same], alleles) %% length(alleles)) + 1L]

  data.frame(
    SNP = paste0("rs", seq_len(n)),
    A1 = a1,
    A2 = a2,
    MAF = runif(n, 0.02, 0.49),
    stringsAsFactors = FALSE
  )
}

make_trait <- function(ref, trait_index) {
  n <- nrow(ref)
  flip <- (seq_len(n) + trait_index) %% 4L == 0L
  a1 <- ref$A1
  a2 <- ref$A2
  a1[flip] <- ref$A2[flip]
  a2[flip] <- ref$A1[flip]

  data.frame(
    SNP = ref$SNP,
    A1 = a1,
    A2 = a2,
    BETA = rnorm(n, sd = 0.05),
    SE = runif(n, 0.02, 0.08),
    P = format.pval(runif(n, 1e-8, 0.999), digits = 8),
    N = sample(8000:12000, n, replace = TRUE),
    INFO = runif(n, 0.91, 1.0),
    MAF = ref$MAF,
    stringsAsFactors = FALSE
  )
}

time_expr <- function(expr) {
  gc()
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  list(value = value, elapsed = elapsed)
}

checksum_df <- function(x) {
  numeric_cols <- vapply(x, is.numeric, logical(1))
  sum(as.matrix(x[, numeric_cols, drop = FALSE]), na.rm = TRUE)
}

bench_prep <- function() {
  with_temp_cwd({
    set.seed(seed)
    ref <- make_ref(n_snp)
    write.table(ref, "ref.txt", row.names = FALSE, quote = FALSE)
    write.table(ref[, c("SNP", "A1", "A2")], "hm3.txt", row.names = FALSE, quote = FALSE)

    trait_files <- paste0("trait", seq_len(n_traits), ".txt")
    trait_names <- paste0("trait", seq_len(n_traits))
    for (i in seq_len(n_traits)) {
      write.table(make_trait(ref, i), trait_files[i], row.names = FALSE, quote = FALSE)
    }

    old_read <- getOption("GenomicSEM.fast_table_read")
    old_join <- getOption("GenomicSEM.fast_snp_join")
    on.exit(options(GenomicSEM.fast_table_read = old_read, GenomicSEM.fast_snp_join = old_join), add = TRUE)

    rows <- list()

    modes <- data.frame(
      fast_table_read = c(FALSE, TRUE, FALSE, TRUE),
      fast_snp_join = c(FALSE, FALSE, TRUE, TRUE)
    )

    for (mode_i in seq_len(nrow(modes))) {
      fast_read <- modes$fast_table_read[mode_i]
      fast_join <- modes$fast_snp_join[mode_i]
      options(GenomicSEM.fast_table_read = fast_read)
      options(GenomicSEM.fast_snp_join = fast_join)

      sumstats_result <- time_expr({
        capture.output({
          result <- suppressWarnings(sumstats(
            files = trait_files,
            ref = "ref.txt",
            trait.names = trait_names,
            se.logit = rep(TRUE, n_traits),
            parallel = FALSE
          ))
        })
        result
      })
      rows[[length(rows) + 1L]] <- data.frame(
        workflow = "sumstats",
        fast_table_read = fast_read,
        fast_snp_join = fast_join,
        n_snp = n_snp,
        n_traits = n_traits,
        n_annot = NA_integer_,
        n_blocks = NA_integer_,
        elapsed_sec = sumstats_result$elapsed,
        checksum = checksum_df(sumstats_result$value),
        stringsAsFactors = FALSE
      )

      munge_result <- time_expr({
        capture.output({
          suppressWarnings(munge(
            files = trait_files,
            hm3 = "hm3.txt",
          trait.names = paste0(trait_names, "_", mode_i),
            parallel = FALSE,
            overwrite = TRUE,
            column.names = list(effect = "BETA")
          ))
        })
      })
      munged_rows <- 0L
      checksum <- 0
      for (trait_name in paste0(trait_names, "_", mode_i)) {
        munged <- read.table(gzfile(paste0(trait_name, ".sumstats.gz")), header = TRUE)
        munged_rows <- munged_rows + nrow(munged)
        checksum <- checksum + checksum_df(munged)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        workflow = "munge",
        fast_table_read = fast_read,
        fast_snp_join = fast_join,
        n_snp = n_snp,
        n_traits = n_traits,
        n_annot = NA_integer_,
        n_blocks = NA_integer_,
        elapsed_sec = munge_result$elapsed,
        checksum = checksum + munged_rows,
        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, rows)
  })
}

bench_ldsc_blocks <- function() {
  set.seed(seed + 1000L)
  weighted.LD <- matrix(rnorm(n_snp * (n_annot + 1L)), nrow = n_snp)
  weighted.chi <- rnorm(n_snp)

  old_loop <- function() {
    select.from <- floor(seq(from = 1, to = nrow(weighted.LD), length.out = n_blocks + 1))
    select.to <- c(select.from[2:n_blocks] - 1, nrow(weighted.LD))
    xty.block.values <- matrix(data = NA_real_, nrow = n_blocks, ncol = ncol(weighted.LD))
    xtx.block.values <- matrix(data = NA_real_, nrow = ncol(weighted.LD) * n_blocks, ncol = ncol(weighted.LD))
    replace.from <- seq(from = 1, to = nrow(xtx.block.values), by = ncol(weighted.LD))
    replace.to <- seq(from = ncol(weighted.LD), to = nrow(xtx.block.values), by = ncol(weighted.LD))
    for (i in seq_len(n_blocks)) {
      rows <- select.from[i]:select.to[i]
      xty.block.values[i, ] <- t(t(weighted.LD[rows, ]) %*% weighted.chi[rows])
      xtx.block.values[replace.from[i]:replace.to[i], ] <- crossprod(weighted.LD[rows, ])
    }
    list(xty = xty.block.values, xtx = xtx.block.values)
  }

  old <- time_expr(old_loop())
  fast <- time_expr(.ldsc_block_products(weighted.LD, weighted.chi, n_blocks))

  data.frame(
    workflow = c("ldsc_block_old_loop", "ldsc_block_fast"),
    fast_table_read = NA,
    fast_snp_join = NA,
    n_snp = n_snp,
    n_traits = NA_integer_,
    n_annot = n_annot,
    n_blocks = n_blocks,
    elapsed_sec = c(old$elapsed, fast$elapsed),
    checksum = c(sum(old$value$xty, old$value$xtx), sum(fast$value$xty.block.values, fast$value$xtx.block.values)),
    stringsAsFactors = FALSE
  )
}

results <- rbind(bench_prep(), bench_ldsc_blocks())
print(results, row.names = FALSE)
