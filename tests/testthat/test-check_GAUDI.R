test_that("GAUDI works", {
    vcf_file <- system.file("extdata", "target.flare.vcf.gz", package="HAUDI")    
    FBM_list <- make_FBM(vcf_file=vcf_file,
                         FBM_pref=tempfile(),
                         chunk_size=100,
                         rds=NULL,
                         minAC=1,
                         geno_format="GT",
                         anc_names=c("pop_1", "pop_2"))

    y_file <- system.file("extdata", "y.txt", package="HAUDI")
    y <- as.numeric(readLines(y_file))
    testthat::expect_no_error({
        GAUDI(FBM_obj=FBM_list$FBM, y=y, FBM_info=FBM_list$info, gamma_vec=1, K=3)
    })
})
