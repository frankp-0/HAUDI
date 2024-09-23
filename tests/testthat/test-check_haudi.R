test_that("HAUDI works", {
  vcf_file <- system.file("extdata", "2way.lanc.vcf.gz", package = "HAUDI")
  fbm <- make_fbm(
    vcf_file = vcf_file,
    fbm_pref = tempfile(),
    chunk_size = 100,
    rds = NULL,
    min_ac = 10,
    geno_format = "GT",
    anc_names = c("pop_1", "pop_2")
  )

  y_file <- system.file("extdata", "y.txt", package = "HAUDI")
  y <- as.numeric(readLines(y_file))
  testthat::expect_no_error({
    haudi(
      fbm_obj = fbm$FBM, y = y, fbm_info = fbm$info,
      gamma = 1, family = "gaussian", K = 3
    )
  })
})
