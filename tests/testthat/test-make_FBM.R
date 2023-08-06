test_that("making FBM works", {
    anc_vcf_file <- system.file("extdata", "test.anc.vcf.gz", package="HAUDI")
    gt_anc <- matrix(
        scan(system.file("extdata", "test.geno.anc.txt", package="HAUDI")),
        byrow = TRUE,
        ncol = 12)
    FBM_pref <- tempfile()
    make_FBM(anc_vcf_file, FBM_pref, 10)
    a <- readRDS(paste0(FBM_pref, ".rds"))
    expect_equal(a$geno[], gt_anc)
})
