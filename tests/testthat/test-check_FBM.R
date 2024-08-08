test_that("FBM genotype data is correct", {
    vcf_file <- system.file("extdata", "target.flare.vcf.gz", package="HAUDI")
    FBM_list <- make_FBM(vcf_file=vcf_file,
                         FBM_pref=tempfile(),
                         chunk_size=100,
                         rds=NULL,
                         minAC=1,
                         geno_format="GT",
                         anc_names=c("pop_1", "pop_2", "pop_3"))

    FBM_rds <- system.file("extdata", "fbm.rds", package="HAUDI")
    FBM_test <- readRDS(FBM_rds)
    expect_equal(FBM_list$FBM[,], FBM_test[,])
})

test_that("FBM info is correct", {
    vcf_file <- system.file("extdata", "target.flare.vcf.gz", package="HAUDI")    
    FBM_list <- make_FBM(vcf_file=vcf_file,
                         FBM_pref=tempfile(),
                         chunk_size=100,
                         rds=NULL,
                         minAC=1,
                         geno_format="GT",
                         anc_names=c("pop_1", "pop_2", "pop_3"))

    FBM_info_file <- system.file("extdata", "fbm_info.tsv", package="HAUDI")    
    FBM_info <- read.table(FBM_info_file, sep='\t', header=T)
    FBM_info$chrom <- as.character(FBM_info$chrom)
    
    expect_equal(FBM_list$info, FBM_info)
})
