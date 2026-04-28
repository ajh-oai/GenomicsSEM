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
.ldsc_block_products_r <- get(".ldsc_block_products_r", envir = genomicssem_ns)
.ldsc_read_table <- get(".ldsc_read_table", envir = genomicssem_ns)
.ldsc_read_chromosome_tables <- get(".ldsc_read_chromosome_tables", envir = genomicssem_ns)
.ldsc_read_m_files <- get(".ldsc_read_m_files", envir = genomicssem_ns)
.ldsc_read_file_list <- get(".ldsc_read_file_list", envir = genomicssem_ns)
.ldsc_read_m_file_list <- get(".ldsc_read_m_file_list", envir = genomicssem_ns)

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
    old_munge_qc <- getOption("GenomicSEM.fast_munge_qc")
    old_sumstats_qc <- getOption("GenomicSEM.fast_sumstats_qc")
    old_munge_engine <- getOption("GenomicSEM.fast_munge_engine")
    old_sumstats_engine <- getOption("GenomicSEM.fast_sumstats_engine")
    on.exit(options(
      GenomicSEM.fast_table_read = old_read,
      GenomicSEM.fast_snp_join = old_join,
      GenomicSEM.fast_munge_qc = old_munge_qc,
      GenomicSEM.fast_sumstats_qc = old_sumstats_qc,
      GenomicSEM.fast_munge_engine = old_munge_engine,
      GenomicSEM.fast_sumstats_engine = old_sumstats_engine
    ), add = TRUE)

    rows <- list()

    modes <- data.frame(
      fast_table_read = c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE),
      fast_snp_join = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
      fast_prep_qc = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
      fast_prep_engine = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
    )

    for (mode_i in seq_len(nrow(modes))) {
      fast_read <- modes$fast_table_read[mode_i]
      fast_join <- modes$fast_snp_join[mode_i]
      fast_prep_qc <- modes$fast_prep_qc[mode_i]
      fast_prep_engine <- modes$fast_prep_engine[mode_i]
      options(GenomicSEM.fast_table_read = fast_read)
      options(GenomicSEM.fast_snp_join = fast_join)
      options(GenomicSEM.fast_munge_qc = fast_prep_qc)
      options(GenomicSEM.fast_sumstats_qc = fast_prep_qc)
      options(GenomicSEM.fast_munge_engine = fast_prep_engine)
      options(GenomicSEM.fast_sumstats_engine = fast_prep_engine)

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
        fast_prep_qc = fast_prep_qc,
        fast_prep_engine = fast_prep_engine,
        fast_ldsc_read = NA,
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
        fast_prep_qc = fast_prep_qc,
        fast_prep_engine = fast_prep_engine,
        fast_ldsc_read = NA,
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

  old_opt <- getOption("GenomicSEM.fast_ldsc_blocks")
  old_threads_opt <- getOption("GenomicSEM.fast_ldsc_threads")
  on.exit(options(
    GenomicSEM.fast_ldsc_blocks = old_opt,
    GenomicSEM.fast_ldsc_threads = old_threads_opt
  ), add = TRUE)

  old <- time_expr(old_loop())
  r_helper <- time_expr(.ldsc_block_products_r(weighted.LD, weighted.chi, n_blocks))
  cores <- parallel::detectCores()
  if (is.na(cores)) {
    cores <- 1L
  }
  thread_counts <- unique(c(1L, min(4L, max(1L, cores - 1L))))
  fast <- lapply(thread_counts, function(n_threads) {
    options(GenomicSEM.fast_ldsc_blocks = TRUE, GenomicSEM.fast_ldsc_threads = n_threads)
    time_expr(.ldsc_block_products(weighted.LD, weighted.chi, n_blocks))
  })

  data.frame(
    workflow = c(
      "ldsc_block_old_loop",
      "ldsc_block_r_helper",
      paste0("ldsc_block_rust_", thread_counts, "t")
    ),
    fast_table_read = NA,
    fast_snp_join = NA,
    fast_prep_qc = NA,
    fast_prep_engine = NA,
    fast_ldsc_read = NA,
    n_snp = n_snp,
    n_traits = NA_integer_,
    n_annot = n_annot,
    n_blocks = n_blocks,
    elapsed_sec = c(old$elapsed, r_helper$elapsed, vapply(fast, function(x) x$elapsed, numeric(1))),
    checksum = c(
      sum(old$value$xty, old$value$xtx),
      sum(r_helper$value$xty.block.values, r_helper$value$xtx.block.values),
      vapply(fast, function(x) {
        sum(x$value$xty.block.values, x$value$xtx.block.values)
      }, numeric(1))
    ),
    stringsAsFactors = FALSE
  )
}

bench_ldsc_read <- function() {
  with_temp_cwd({
    set.seed(seed + 2000L)
    n_chrom <- 4L
    n_per_chrom <- ceiling(n_snp / n_chrom)
    for (chr in seq_len(n_chrom)) {
      n <- if (chr < n_chrom) n_per_chrom else n_snp - n_per_chrom * (n_chrom - 1L)
      snp_offset <- (chr - 1L) * n_per_chrom
      ld <- data.frame(
        CHR = chr,
        SNP = paste0("rs", snp_offset + seq_len(n)),
        BP = snp_offset + seq_len(n),
        CM = 0,
        MAF = runif(n, 0.05, 0.5),
        L2 = runif(n, 1, 20),
        stringsAsFactors = FALSE
      )
      write.table(ld, gzfile(paste0(chr, ".l2.ldscore.gz")), row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(data.frame(V1 = n), paste0(chr, ".l2.M_5_50"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    trait <- data.frame(
      SNP = paste0("rs", seq_len(n_snp)),
      N = sample(8000:12000, n_snp, replace = TRUE),
      Z = rnorm(n_snp),
      A1 = sample(c("A", "C", "G", "T"), n_snp, replace = TRUE),
      stringsAsFactors = FALSE
    )
    write.table(trait, gzfile("trait.sumstats.gz"), row.names = FALSE, quote = FALSE, sep = "\t")

    old <- getOption("GenomicSEM.fast_ldsc_read")
    on.exit(options(GenomicSEM.fast_ldsc_read = old), add = TRUE)

    rows <- list()
    for (fast in c(FALSE, TRUE)) {
      options(GenomicSEM.fast_ldsc_read = fast)
      timed <- time_expr({
        ld <- .ldsc_read_chromosome_tables(".", ".l2.ldscore.gz", seq_len(n_chrom))
        m <- .ldsc_read_m_files(".", seq_len(n_chrom))
        trait <- .ldsc_read_table("trait.sumstats.gz")
        list(ld = ld, m = m, trait = trait)
      })
      rows[[length(rows) + 1L]] <- data.frame(
        workflow = "ldsc_read",
        fast_table_read = NA,
        fast_snp_join = NA,
        fast_prep_qc = NA,
        fast_prep_engine = NA,
        fast_ldsc_read = fast,
        n_snp = n_snp,
        n_traits = 1L,
        n_annot = 1L,
        n_blocks = NA_integer_,
        elapsed_sec = timed$elapsed,
        checksum = checksum_df(timed$value$ld) + checksum_df(timed$value$m) + checksum_df(timed$value$trait),
        stringsAsFactors = FALSE
      )

      timed_list <- time_expr({
        ld <- .ldsc_read_file_list(paste0(seq_len(n_chrom), ".l2.ldscore.gz"))
        m <- .ldsc_read_m_file_list(paste0(seq_len(n_chrom), ".l2.M_5_50"))
        list(ld = ld, m = m)
      })
      rows[[length(rows) + 1L]] <- data.frame(
        workflow = "s_ldsc_read",
        fast_table_read = NA,
        fast_snp_join = NA,
        fast_prep_qc = NA,
        fast_prep_engine = NA,
        fast_ldsc_read = fast,
        n_snp = n_snp,
        n_traits = NA_integer_,
        n_annot = 1L,
        n_blocks = NA_integer_,
        elapsed_sec = timed_list$elapsed,
        checksum = checksum_df(timed_list$value$ld) + checksum_df(timed_list$value$m),
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  })
}

results <- rbind(bench_prep(), bench_ldsc_blocks(), bench_ldsc_read())
print(results, row.names = FALSE)
