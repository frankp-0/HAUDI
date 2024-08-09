test_that("GAUDI works", {
  vcf_file <- system.file("extdata", "target.flare.vcf.gz", package = "HAUDI")
  fbm <- make_fbm(
    vcf_file = vcf_file,
    fbm_pref = tempfile(),
    chunk_size = 100,
    rds = NULL,
    min_ac = 1,
    geno_format = "GT",
    anc_names = c("pop_1", "pop_2")
  )

  y_file <- system.file("extdata", "y.txt", package = "HAUDI")
  y <- as.numeric(readLines(y_file))
  testthat::expect_no_error({
    gaudi(
      fbm_obj = fbm$FBM, y = y, fbm_info = fbm$info,
      gamma_vec = 1, k = 3
    )
  })
})


vcf_file <- system.file("extdata", "target.flare.vcf.gz", package = "HAUDI")
vcf_file <- vcf_file
fbm_pref <- tempfile()
chunk_size <- 100
rds <- NULL
min_ac <- 1
geno_format <- "GT"
anc_names <- c("pop_1", "pop_2")
