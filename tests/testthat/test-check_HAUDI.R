test_that("HAUDI works", {
  y_file <- system.file("extdata", "y.txt", package = "HAUDI")
  y <- as.numeric(readLines(y_file))

  rds_file <- system.file("extdata", "fbm.rds", package = "HAUDI")
  fbm <- readRDS(rds_file)

  info_file <- system.file("extdata", "fbm_info.tsv", package = "HAUDI")
  fbm_info <- read.table(info_file, header = TRUE, sep = "\t")
  fbm_info$chrom <- as.character(fbm_info$chrom)

  testthat::expect_no_error({
    haudi(fbm_obj = fbm, y = y, fbm_info = fbm_info, gamma = 1, family = "gaussian", K = 3)
  })
})
