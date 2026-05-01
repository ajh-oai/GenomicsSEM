library(GenomicSEM)

genomicssem_ns <- asNamespace("GenomicSEM")
.ldsc_block_products <- get(".ldsc_block_products", envir = genomicssem_ns)
.ldsc_block_products_r <- get(".ldsc_block_products_r", envir = genomicssem_ns)
.read_sumstats_table <- get(".read_sumstats_table", envir = genomicssem_ns)
.ldsc_read_chromosome_tables <- get(".ldsc_read_chromosome_tables", envir = genomicssem_ns)
.ldsc_read_m_files <- get(".ldsc_read_m_files", envir = genomicssem_ns)
.ldsc_read_file_list <- get(".ldsc_read_file_list", envir = genomicssem_ns)
.ldsc_read_m_file_list <- get(".ldsc_read_m_file_list", envir = genomicssem_ns)
.munge_fused_fast <- get(".munge_fused_fast", envir = genomicssem_ns)
.sumstats_fused_fast <- get(".sumstats_fused_fast", envir = genomicssem_ns)

with_temp_cwd <- function(code) {
  old <- getwd()
  tmp <- tempfile("genomicssem-prep-fast-")
  dir.create(tmp)
  on.exit({
    setwd(old)
    unlink(tmp, recursive = TRUE)
  })
  setwd(tmp)
  force(code)
}

check_block_products <- function() {
  set.seed(11)
  weighted.LD <- matrix(rnorm(300 * 6), nrow = 300, ncol = 6)
  weighted.chi <- rnorm(300)
  colnames(weighted.LD) <- paste0("A", seq_len(ncol(weighted.LD)))
  n.blocks <- 23

  old_opt <- getOption("GenomicSEM.fast_ldsc_blocks")
  on.exit(options(GenomicSEM.fast_ldsc_blocks = old_opt), add = TRUE)

  options(GenomicSEM.fast_ldsc_blocks = FALSE)
  fallback <- .ldsc_block_products(weighted.LD, weighted.chi, n.blocks)

  options(GenomicSEM.fast_ldsc_blocks = TRUE)
  fast <- .ldsc_block_products(weighted.LD, weighted.chi, n.blocks)
  direct_r <- .ldsc_block_products_r(weighted.LD, weighted.chi, n.blocks)

  select.from <- floor(seq(from = 1, to = nrow(weighted.LD), length.out = n.blocks + 1))
  select.to <- c(select.from[2:n.blocks] - 1, nrow(weighted.LD))
  xty.block.values <- matrix(NA_real_, nrow = n.blocks, ncol = ncol(weighted.LD))
  xtx.block.values <- matrix(NA_real_, nrow = ncol(weighted.LD) * n.blocks, ncol = ncol(weighted.LD))
  replace.from <- seq(from = 1, to = nrow(xtx.block.values), by = ncol(weighted.LD))
  replace.to <- seq(from = ncol(weighted.LD), to = nrow(xtx.block.values), by = ncol(weighted.LD))

  for (i in seq_len(n.blocks)) {
    rows <- select.from[i]:select.to[i]
    xty.block.values[i, ] <- t(t(weighted.LD[rows, ]) %*% weighted.chi[rows])
    xtx.block.values[replace.from[i]:replace.to[i], ] <- t(weighted.LD[rows, ]) %*% weighted.LD[rows, ]
  }

  stopifnot(identical(fallback, direct_r))
  stopifnot(max(abs(xty.block.values - fast$xty.block.values)) < 1e-10)
  stopifnot(max(abs(xtx.block.values - fast$xtx.block.values)) < 1e-10)
  stopifnot(max(abs(fallback$xty.block.values - fast$xty.block.values)) < 1e-10)
  stopifnot(max(abs(fallback$xtx.block.values - fast$xtx.block.values)) < 1e-10)
  stopifnot(identical(fallback$delete.from, fast$delete.from))
  stopifnot(identical(fallback$delete.to, fast$delete.to))
}

check_reader_preserves_p_character <- function() {
  path <- tempfile(fileext = ".txt")
  writeLines(c(
    "SNP A1 A2 effect SE P N",
    "rs1 A G 0.12 0.04 1e-320 10000"
  ), path)

  old <- getOption("GenomicSEM.fast_table_read")
  on.exit(options(GenomicSEM.fast_table_read = old))
  options(GenomicSEM.fast_table_read = TRUE)

  parsed <- .read_sumstats_table(path, p_as_character = TRUE)
  stopifnot(identical(parsed$P, "1e-320"))
}

make_ref <- function() {
  data.frame(
    SNP = paste0("rs", 1:6),
    A1 = c("A", "C", "G", "T", "A", "C"),
    A2 = c("G", "T", "A", "C", "T", "G"),
    MAF = c(0.20, 0.35, 0.45, 0.12, 0.02, 0.30),
    stringsAsFactors = FALSE
  )
}

make_sumstats <- function() {
  data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs_missing"),
    A1 = c("A", "T", "G", "C", "A", "G", "A"),
    A2 = c("G", "C", "A", "T", "T", "C", "G"),
    BETA = c(0.10, -0.05, 0.20, -0.12, 0.08, 0.03, 0.50),
    SE = c(0.03, 0.04, 0.05, 0.06, 0.04, 0.02, 0.10),
    P = c("0.001", "0.02", "0.15", "0.30", "0.80", "0.04", "0.90"),
    N = c(10000, 9500, 9000, 11000, 8700, 9200, 8000),
    INFO = c(0.99, 0.98, 0.95, 0.92, 0.99, 0.91, 0.99),
    MAF = c(0.20, 0.35, 0.45, 0.12, 0.02, 0.30, 0.22),
    stringsAsFactors = FALSE
  )
}

check_sumstats_fast_reader <- function() {
  with_temp_cwd({
    ref <- make_ref()
    ss <- make_sumstats()
    write.table(ref, "ref.txt", row.names = FALSE, quote = FALSE)
    write.table(ss, "trait.txt", row.names = FALSE, quote = FALSE)

    old_read <- getOption("GenomicSEM.fast_table_read")
    old_join <- getOption("GenomicSEM.fast_snp_join")
    old_sumstats_qc <- getOption("GenomicSEM.fast_sumstats_qc")
    old_sumstats_engine <- getOption("GenomicSEM.fast_sumstats_engine")
    on.exit(options(
      GenomicSEM.fast_table_read = old_read,
      GenomicSEM.fast_snp_join = old_join,
      GenomicSEM.fast_sumstats_qc = old_sumstats_qc,
      GenomicSEM.fast_sumstats_engine = old_sumstats_engine
    ), add = TRUE)

    options(GenomicSEM.fast_table_read = FALSE, GenomicSEM.fast_snp_join = FALSE,
            GenomicSEM.fast_sumstats_qc = FALSE, GenomicSEM.fast_sumstats_engine = FALSE)
    fallback <- suppressWarnings(sumstats(
      files = "trait.txt",
      ref = "ref.txt",
      trait.names = "trait",
      se.logit = TRUE,
      parallel = FALSE
    ))

    options(GenomicSEM.fast_table_read = TRUE, GenomicSEM.fast_snp_join = TRUE,
            GenomicSEM.fast_sumstats_qc = TRUE, GenomicSEM.fast_sumstats_engine = TRUE)
    fast <- suppressWarnings(sumstats(
      files = "trait.txt",
      ref = "ref.txt",
      trait.names = "trait",
      se.logit = TRUE,
      parallel = FALSE
    ))

    stopifnot(identical(names(fallback), names(fast)))
    stopifnot(max(abs(fallback$beta.trait - fast$beta.trait), na.rm = TRUE) < 1e-12)
    stopifnot(max(abs(fallback$se.trait - fast$se.trait), na.rm = TRUE) < 1e-12)
  })
}

check_sumstats_fused_engine <- function() {
  with_temp_cwd({
    ref <- make_ref()
    ss <- make_sumstats()
    write.table(ref, "ref.txt", row.names = FALSE, quote = FALSE)
    write.table(ss, gzfile("trait.txt.gz"), row.names = FALSE, quote = FALSE)

    old <- options(
      GenomicSEM.fast_table_read = FALSE,
      GenomicSEM.fast_snp_join = FALSE,
      GenomicSEM.fast_sumstats_qc = FALSE,
      GenomicSEM.fast_sumstats_engine = FALSE
    )
    on.exit(options(old), add = TRUE)

    fallback <- suppressWarnings(sumstats(
      files = "trait.txt.gz",
      ref = "ref.txt",
      trait.names = "trait",
      se.logit = TRUE,
      parallel = FALSE
    ))

    options(GenomicSEM.fast_sumstats_engine = TRUE)
    log.file <- file("sumstats_fused.log", open = "wt")
    fused <- .sumstats_fused_fast(
      filename = "trait.txt.gz",
      trait.name = "trait",
      N = NA_real_,
      keep.indel = FALSE,
      OLS = FALSE,
      beta = FALSE,
      info.filter = 0.6,
      linprob = FALSE,
      se.logit = TRUE,
      name.beta = "beta.trait",
      name.se = "se.trait",
      ref = ref,
      ref2 = "ref.txt",
      log.file = log.file,
      direct.filter = FALSE
    )
    close(log.file)

    stopifnot(!is.null(fused))
    stopifnot(identical(fallback$SNP, fused$SNP))
    stopifnot(max(abs(fallback$beta.trait - fused$beta.trait), na.rm = TRUE) < 1e-12)
    stopifnot(max(abs(fallback$se.trait - fused$se.trait), na.rm = TRUE) < 1e-12)
  })
}

check_sumstats_fused_batch_engine <- function() {
  with_temp_cwd({
    ref <- make_ref()
    ss1 <- make_sumstats()
    ss2 <- make_sumstats()
    ss2$P[ss2$SNP == "rs6"] <- NA_character_
    write.table(ref, "ref.txt", row.names = FALSE, quote = FALSE)
    write.table(ss1, gzfile("trait1.txt.gz"), row.names = FALSE, quote = FALSE)
    write.table(ss2, gzfile("trait2.txt.gz"), row.names = FALSE, quote = FALSE)

    old <- options(
      GenomicSEM.fast_table_read = FALSE,
      GenomicSEM.fast_snp_join = FALSE,
      GenomicSEM.fast_sumstats_qc = FALSE,
      GenomicSEM.fast_sumstats_engine = FALSE,
      GenomicSEM.fast_sumstats_threads = NULL
    )
    on.exit(options(old), add = TRUE)

    fallback <- suppressWarnings(sumstats(
      files = c("trait1.txt.gz", "trait2.txt.gz"),
      ref = "ref.txt",
      trait.names = c("trait1", "trait2"),
      se.logit = c(TRUE, TRUE),
      parallel = FALSE
    ))

    options(
      GenomicSEM.fast_sumstats_engine = TRUE,
      GenomicSEM.fast_sumstats_threads = 2L
    )
    fast <- suppressWarnings(sumstats(
      files = c("trait1.txt.gz", "trait2.txt.gz"),
      ref = "ref.txt",
      trait.names = c("trait1", "trait2"),
      se.logit = c(TRUE, TRUE),
      parallel = TRUE,
      cores = 2L
    ))

    stopifnot(identical(names(fallback), names(fast)))
    stopifnot(identical(fallback$SNP, fast$SNP))
    stopifnot(!("rs6" %in% fast$SNP))
    numeric_cols <- setdiff(names(fallback), c("SNP", "A1", "A2"))
    stopifnot(max(abs(as.matrix(fallback[, numeric_cols]) - as.matrix(fast[, numeric_cols])), na.rm = TRUE) < 1e-12)
    stopifnot(file.exists("trait1_sumstats.log"))
    stopifnot(file.exists("trait2_sumstats.log"))
  })
}

check_munge_fast_reader <- function() {
  with_temp_cwd({
    ref <- make_ref()
    ss <- make_sumstats()
    write.table(ref[, c("SNP", "A1", "A2")], "hm3.txt", row.names = FALSE, quote = FALSE)
    write.table(ss, "trait.txt", row.names = FALSE, quote = FALSE)

    old_read <- getOption("GenomicSEM.fast_table_read")
    old_join <- getOption("GenomicSEM.fast_snp_join")
    old_munge_qc <- getOption("GenomicSEM.fast_munge_qc")
    old_munge_engine <- getOption("GenomicSEM.fast_munge_engine")
    on.exit(options(
      GenomicSEM.fast_table_read = old_read,
      GenomicSEM.fast_snp_join = old_join,
      GenomicSEM.fast_munge_qc = old_munge_qc,
      GenomicSEM.fast_munge_engine = old_munge_engine
    ), add = TRUE)

    options(GenomicSEM.fast_table_read = FALSE, GenomicSEM.fast_snp_join = FALSE,
            GenomicSEM.fast_munge_qc = FALSE, GenomicSEM.fast_munge_engine = FALSE)
    suppressWarnings(munge(
      files = "trait.txt",
      hm3 = "hm3.txt",
      trait.names = "trait_slow",
      parallel = FALSE,
      overwrite = TRUE,
      column.names = list(effect = "BETA")
    ))
    fallback <- read.table(gzfile("trait_slow.sumstats.gz"), header = TRUE)

    options(GenomicSEM.fast_table_read = TRUE, GenomicSEM.fast_snp_join = TRUE,
            GenomicSEM.fast_munge_qc = TRUE, GenomicSEM.fast_munge_engine = TRUE)
    suppressWarnings(munge(
      files = "trait.txt",
      hm3 = "hm3.txt",
      trait.names = "trait_fast",
      parallel = FALSE,
      overwrite = TRUE,
      column.names = list(effect = "BETA")
    ))
    fast <- read.table(gzfile("trait_fast.sumstats.gz"), header = TRUE)

    stopifnot(identical(fallback$SNP, fast$SNP))
    stopifnot(identical(fallback$N, fast$N))
    stopifnot(identical(fallback$A1, fast$A1))
    stopifnot(identical(fallback$A2, fast$A2))
    stopifnot(max(abs(fallback$Z - fast$Z), na.rm = TRUE) < 1e-5)
  })
}

check_munge_fused_engine <- function() {
  with_temp_cwd({
    ref <- make_ref()
    ss <- make_sumstats()
    hm3 <- ref[, c("SNP", "A1", "A2")]
    write.table(hm3, "hm3.txt", row.names = FALSE, quote = FALSE)
    write.table(ss, gzfile("trait.txt.gz"), row.names = FALSE, quote = FALSE)

    old <- options(
      GenomicSEM.fast_table_read = FALSE,
      GenomicSEM.fast_snp_join = FALSE,
      GenomicSEM.fast_munge_qc = FALSE,
      GenomicSEM.fast_munge_engine = FALSE
    )
    on.exit(options(old), add = TRUE)

    suppressWarnings(munge(
      files = "trait.txt.gz",
      hm3 = "hm3.txt",
      trait.names = "trait_slow",
      parallel = FALSE,
      overwrite = TRUE,
      column.names = list(effect = "BETA")
    ))
    fallback <- read.table(gzfile("trait_slow.sumstats.gz"), header = TRUE)

    options(GenomicSEM.fast_munge_engine = TRUE)
    log.file <- file("munge_fused.log", open = "wt")
    fused_result <- .munge_fused_fast(
      filename = "trait.txt.gz",
      trait.name = "trait_fused",
      N = NA_real_,
      ref = hm3,
      hm3 = "hm3.txt",
      info.filter = 0.9,
      maf.filter = 0.01,
      column.names = list(effect = "BETA"),
      overwrite = TRUE,
      log.file = log.file
    )
    close(log.file)

    stopifnot(!is.null(fused_result))
    fused <- read.table(gzfile("trait_fused.sumstats.gz"), header = TRUE)
    stopifnot(identical(fallback$SNP, fused$SNP))
    stopifnot(identical(fallback$N, fused$N))
    stopifnot(identical(fallback$A1, fused$A1))
    stopifnot(identical(fallback$A2, fused$A2))
    stopifnot(max(abs(fallback$Z - fused$Z), na.rm = TRUE) < 1e-5)
    stopifnot(!file.exists("trait_fused.sumstats"))
  })
}

check_munge_fused_batch_engine <- function() {
  with_temp_cwd({
    ref <- make_ref()
    ss1 <- make_sumstats()
    ss2 <- make_sumstats()
    hm3 <- ref[, c("SNP", "A1", "A2")]
    write.table(hm3, "hm3.txt", row.names = FALSE, quote = FALSE)
    write.table(ss1, gzfile("trait1.txt.gz"), row.names = FALSE, quote = FALSE)
    write.table(ss2, gzfile("trait2.txt.gz"), row.names = FALSE, quote = FALSE)

    old <- options(
      GenomicSEM.fast_table_read = FALSE,
      GenomicSEM.fast_snp_join = FALSE,
      GenomicSEM.fast_munge_qc = FALSE,
      GenomicSEM.fast_munge_engine = FALSE,
      GenomicSEM.fast_munge_threads = NULL
    )
    on.exit(options(old), add = TRUE)

    suppressWarnings(munge(
      files = c("trait1.txt.gz", "trait2.txt.gz"),
      hm3 = "hm3.txt",
      trait.names = c("trait1_slow", "trait2_slow"),
      parallel = FALSE,
      overwrite = TRUE,
      column.names = list(effect = "BETA")
    ))
    fallback1 <- read.table(gzfile("trait1_slow.sumstats.gz"), header = TRUE)
    fallback2 <- read.table(gzfile("trait2_slow.sumstats.gz"), header = TRUE)

    options(
      GenomicSEM.fast_munge_engine = TRUE,
      GenomicSEM.fast_munge_threads = 2L
    )
    suppressWarnings(munge(
      files = c("trait1.txt.gz", "trait2.txt.gz"),
      hm3 = "hm3.txt",
      trait.names = c("trait1_fast", "trait2_fast"),
      parallel = TRUE,
      cores = 2L,
      overwrite = TRUE,
      column.names = list(effect = "BETA")
    ))
    fast1 <- read.table(gzfile("trait1_fast.sumstats.gz"), header = TRUE)
    fast2 <- read.table(gzfile("trait2_fast.sumstats.gz"), header = TRUE)

    stopifnot(identical(fallback1$SNP, fast1$SNP))
    stopifnot(identical(fallback2$SNP, fast2$SNP))
    stopifnot(identical(fallback1$N, fast1$N))
    stopifnot(identical(fallback2$N, fast2$N))
    stopifnot(max(abs(fallback1$Z - fast1$Z), na.rm = TRUE) < 1e-5)
    stopifnot(max(abs(fallback2$Z - fast2$Z), na.rm = TRUE) < 1e-5)
    stopifnot(file.exists("trait1_fast_munge.log"))
    stopifnot(file.exists("trait2_fast_munge.log"))
    stopifnot(!file.exists("trait1_fast.sumstats"))
    stopifnot(!file.exists("trait2_fast.sumstats"))
  })
}

check_ldsc_fast_reader <- function() {
  with_temp_cwd({
    for (chr in 1:2) {
      ld <- data.frame(
        CHR = chr,
        SNP = paste0("rs", chr, "_", 1:3),
        BP = chr * 1000 + 1:3,
        CM = 0,
        MAF = c(0.1, 0.2, 0.3),
        L2 = c(1.2, 1.3, 1.4),
        stringsAsFactors = FALSE
      )
      write.table(ld, gzfile(paste0(chr, ".l2.ldscore.gz")), row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(data.frame(V1 = chr + 10), paste0(chr, ".l2.M_5_50"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    old <- getOption("GenomicSEM.fast_ldsc_read")
    on.exit(options(GenomicSEM.fast_ldsc_read = old), add = TRUE)

    options(GenomicSEM.fast_ldsc_read = FALSE)
    fallback_ld <- .ldsc_read_chromosome_tables(".", ".l2.ldscore.gz", 1:2)
    fallback_m <- .ldsc_read_m_files(".", 1:2)
    fallback_ld_list <- .ldsc_read_file_list(paste0(1:2, ".l2.ldscore.gz"))
    fallback_m_list <- .ldsc_read_m_file_list(paste0(1:2, ".l2.M_5_50"))

    options(GenomicSEM.fast_ldsc_read = TRUE)
    fast_ld <- .ldsc_read_chromosome_tables(".", ".l2.ldscore.gz", 1:2)
    fast_m <- .ldsc_read_m_files(".", 1:2)
    fast_ld_list <- .ldsc_read_file_list(paste0(1:2, ".l2.ldscore.gz"))
    fast_m_list <- .ldsc_read_m_file_list(paste0(1:2, ".l2.M_5_50"))

    stopifnot(isTRUE(all.equal(as.data.frame(fallback_ld), as.data.frame(fast_ld), check.attributes = FALSE)))
    stopifnot(identical(as.numeric(as.matrix(fallback_m)), as.numeric(as.matrix(fast_m))))
    stopifnot(isTRUE(all.equal(as.data.frame(fallback_ld_list), as.data.frame(fast_ld_list), check.attributes = FALSE)))
    stopifnot(identical(as.numeric(as.matrix(fallback_m_list)), as.numeric(as.matrix(fast_m_list))))
  })
}

check_ldsc_full_fast_reader <- function() {
  with_temp_cwd({
    set.seed(5)
    for (chr in 1:2) {
      n <- 300
      idx <- (chr - 1L) * n + seq_len(n)
      ld <- data.frame(
        CHR = chr,
        SNP = paste0("rs", idx),
        BP = idx,
        CM = 0,
        MAF = runif(n, 0.05, 0.49),
        L2 = runif(n, 1, 30),
        stringsAsFactors = FALSE
      )
      write.table(ld, gzfile(paste0(chr, ".l2.ldscore.gz")), row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(data.frame(V1 = n), paste0(chr, ".l2.M_5_50"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    traits <- paste0("trait", 1:2, ".sumstats.gz")
    for (trait_i in 1:2) {
      sumstats <- data.frame(
        SNP = paste0("rs", 1:600),
        N = sample(9000:12000, 600, TRUE),
        Z = rnorm(600),
        A1 = sample(c("A", "C", "G", "T"), 600, TRUE),
        stringsAsFactors = FALSE
      )
      write.table(sumstats, gzfile(traits[trait_i]), row.names = FALSE, quote = FALSE, sep = "\t")
    }

    old <- getOption("GenomicSEM.fast_ldsc_read")
    on.exit(options(GenomicSEM.fast_ldsc_read = old), add = TRUE)

    run_ldsc <- function(fast) {
      options(GenomicSEM.fast_ldsc_read = fast)
      capture.output({
        out <- suppressWarnings(ldsc(
          traits = traits,
          sample.prev = c(NA, NA),
          population.prev = c(NA, NA),
          ld = ".",
          wld = ".",
          trait.names = c("t1", "t2"),
          sep_weights = FALSE,
          chr = 2,
          n.blocks = 10,
          ldsc.log = paste0("ldsc_", fast)
        ))
      })
      out
    }

    fallback <- run_ldsc(FALSE)
    fast <- run_ldsc(TRUE)

    stopifnot(isTRUE(all.equal(fallback$S, fast$S, tolerance = 1e-10, check.attributes = FALSE)))
    stopifnot(isTRUE(all.equal(fallback$V, fast$V, tolerance = 1e-10, check.attributes = FALSE)))
    stopifnot(isTRUE(all.equal(fallback$I, fast$I, tolerance = 1e-10, check.attributes = FALSE)))
  })
}

check_block_products()
check_reader_preserves_p_character()
check_sumstats_fast_reader()
check_sumstats_fused_engine()
check_sumstats_fused_batch_engine()
check_munge_fast_reader()
check_munge_fused_engine()
check_munge_fused_batch_engine()
check_ldsc_fast_reader()
check_ldsc_full_fast_reader()

cat("prep fast path tests passed\n")
