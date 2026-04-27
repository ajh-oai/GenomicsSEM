library(GenomicSEM)

genomicssem_ns <- asNamespace("GenomicSEM")
.ldsc_block_products <- get(".ldsc_block_products", envir = genomicssem_ns)
.read_sumstats_table <- get(".read_sumstats_table", envir = genomicssem_ns)

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

  fast <- .ldsc_block_products(weighted.LD, weighted.chi, n.blocks)

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

  stopifnot(max(abs(xty.block.values - fast$xty.block.values)) < 1e-10)
  stopifnot(max(abs(xtx.block.values - fast$xtx.block.values)) < 1e-10)
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

    old <- getOption("GenomicSEM.fast_table_read")
    on.exit(options(GenomicSEM.fast_table_read = old), add = TRUE)

    options(GenomicSEM.fast_table_read = FALSE)
    fallback <- suppressWarnings(sumstats(
      files = "trait.txt",
      ref = "ref.txt",
      trait.names = "trait",
      se.logit = TRUE,
      parallel = FALSE
    ))

    options(GenomicSEM.fast_table_read = TRUE)
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

check_munge_fast_reader <- function() {
  with_temp_cwd({
    ref <- make_ref()
    ss <- make_sumstats()
    write.table(ref[, c("SNP", "A1", "A2")], "hm3.txt", row.names = FALSE, quote = FALSE)
    write.table(ss, "trait.txt", row.names = FALSE, quote = FALSE)

    old <- getOption("GenomicSEM.fast_table_read")
    on.exit(options(GenomicSEM.fast_table_read = old), add = TRUE)

    options(GenomicSEM.fast_table_read = FALSE)
    suppressWarnings(munge(
      files = "trait.txt",
      hm3 = "hm3.txt",
      trait.names = "trait_slow",
      parallel = FALSE,
      overwrite = TRUE,
      column.names = list(effect = "BETA")
    ))
    fallback <- read.table(gzfile("trait_slow.sumstats.gz"), header = TRUE)

    options(GenomicSEM.fast_table_read = TRUE)
    suppressWarnings(munge(
      files = "trait.txt",
      hm3 = "hm3.txt",
      trait.names = "trait_fast",
      parallel = FALSE,
      overwrite = TRUE,
      column.names = list(effect = "BETA")
    ))
    fast <- read.table(gzfile("trait_fast.sumstats.gz"), header = TRUE)

    stopifnot(identical(fallback, fast))
  })
}

check_block_products()
check_reader_preserves_p_character()
check_sumstats_fast_reader()
check_munge_fast_reader()

cat("prep fast path tests passed\n")
