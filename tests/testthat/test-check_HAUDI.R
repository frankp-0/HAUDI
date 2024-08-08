test_that("HAUDI works", {
    y_file <- system.file("extdata", "y.txt", package="HAUDI")    
    y <- as.numeric(readLines(y_file))

    rds_file <- system.file("extdata", "fbm.rds", package="HAUDI")
    FBM <- readRDS(rds_file)

    info_file <- system.file("extdata", "fbm_info.tsv", package="HAUDI")
    FBM_info <- read.table(info_file, header=T, sep='\t')
    FBM_info$chrom <- as.character(FBM_info$chrom)
    
    testthat::expect_no_error({
        HAUDI(FBM_obj=FBM, y=y, FBM_info=FBM_info, gamma=1, family="gaussian", K=3)
    })
})
