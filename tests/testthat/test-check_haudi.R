test_that("HAUDI works", {
  y_file <- system.file("extdata", "y.txt", package = "HAUDI")
  y <- as.numeric(readLines(y_file))

  bk_file <- system.file("extdata", "fbm.bk", package = "HAUDI")
  bk_file <- sub(x = bk_file, pattern = ".bk", replacement = "")
  code_dosage <- rep(NA_real_, 256)
  code_dosage[1:201] <- seq(0, 2, length.out = 201)
  fbm <- bigstatsr::FBM.code256(
    nrow = 40, ncol = 1200,
    code = code_dosage,
    backingfile = bk_file,
    create_bk = FALSE,
    is_read_only = TRUE
  )

  info_file <- system.file("extdata", "fbm_info.tsv", package = "HAUDI")
  fbm_info <- read.table(info_file, header = TRUE, sep = "\t")
  fbm_info$chrom <- as.character(fbm_info$chrom)

  testthat::expect_no_error({
    haudi(
      fbm_obj = fbm, y = y, fbm_info = fbm_info,
      gamma = 1, family = "gaussian", K = 3
    )
  })
})
