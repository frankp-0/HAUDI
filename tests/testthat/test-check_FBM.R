test_that("FBM genotype data is correct", {
  vcf_file <- system.file("extdata", "target.flare.vcf.gz", package = "HAUDI")
  fbm_list <- make_fbm(
    vcf_file = vcf_file,
    fbm_pref = tempfile(),
    chunk_size = 100,
    rds = NULL,
    min_ac = 1,
    geno_format = "GT",
    anc_names = c("pop_1", "pop_2", "pop_3")
  )

  fbm_rds <- system.file("extdata", "fbm.rds", package = "HAUDI")
  fbm_test <- readRDS(fbm_rds)
  expect_equal(fbm_list$FBM[, ], fbm_test[, ])
})

test_that("FBM info is correct", {
  vcf_file <- system.file("extdata", "target.flare.vcf.gz", package = "HAUDI")
  fbm_list <- make_fbm(
    vcf_file = vcf_file,
    fbm_pref = tempfile(),
    chunk_size = 100,
    rds = NULL,
    min_ac = 1,
    geno_format = "GT",
    anc_names = c("pop_1", "pop_2", "pop_3")
  )

  fbm_info_file <- system.file("extdata", "fbm_info.tsv", package = "HAUDI")
  fbm_info <- read.table(fbm_info_file, sep = "\t", header = TRUE)
  fbm_info$chrom <- as.character(fbm_info$chrom)

  expect_equal(fbm_list$info, fbm_info)
})
