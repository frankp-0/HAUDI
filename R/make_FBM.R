mat_substr_to_int <- function(mat, start, end){
    mat <- apply(mat, 2, function(x) {
        x <- substr(x,start,end)
        mode(x) <- "integer"
        return(x)})
    return(mat)
}

interleave_matrices <- function(L){
    return(do.call(rbind, L)[order(sequence(sapply(L, nrow))),])
}


initialize_FBM <- function(FBM_pref, nrow){
    code_dosage <- rep(NA_real_, 256)
    code_dosage[1:201] <- seq(0, 2, length.out = 201)
    backing_file <- paste0(FBM_pref, ".bk")
    if(file.exists(backing_file)){file.remove(backing_file)}
    anc_FBM <- FBM.code256(nrow = nrow,
                           ncol = 0,
                           backingfile = FBM_pref,
                           code = code_dosage)
    return(anc_FBM)
}

get_anc_gt <- function(vcf, anc){
    gt1 <- mat_substr_to_int(assays(vcf)$GT, 1, 1)
    gt2 <- mat_substr_to_int(assays(vcf)$GT, 3, 3)
    anc_gt <- lapply(anc, function(a) {
        (gt1 * (assays(vcf)$AN1 == a)) + (gt2 * (assays(vcf)$AN2 == a))})
    anc_gt %<>% append(list(gt1 + gt2)) %>% interleave_matrices()
    anc_gt <- t(anc_gt)
    return(anc_gt)
}

##' Make FBM containing ancestry-genotypes and save in an RDS file with supplementary info
##'
##' @param anc_vcf_file a file path to a VCF file with phased genotypes and estimated local ancestry, i.e. the output of flare (https://doi.org/10.1101/2022.08.02.502540)
##' @param FBM_pref a file path to store FBM and RDS file to. Omit the file extension
##' @param chunk_size an integer indicating the max number of VCF records to read at a time
##' @param minAC an integer indicating the minimum allele count (per-ancestry) to retain
##' @author Frank Ockerman
##' @importFrom VariantAnnotation VcfFile scanVcfHeader samples readVcf ref alt
##' @importFrom SummarizedExperiment assays
##' @importFrom GenomeInfoDb seqnames
##' @importFrom IRanges ranges
##' @importFrom bigstatsr FBM.code256 sub_bk
##' @importFrom magrittr `%>%` `%<>%`
##' @export
make_FBM <- function(anc_vcf_file, FBM_pref, chunk_size, minAC){
    ## Reference and open VCF file
    tab <- VcfFile(anc_vcf_file, yieldSize = chunk_size)
    open(tab)    

    ## Extract VCF info
    anc_names <- names(scanVcfHeader(tab)@header@listData$ANCESTRY)
    anc <- as.integer(unlist(scanVcfHeader(tab)@header@listData$ANCESTRY))
    samples <- samples(scanVcfHeader(tab))

    ## Initialize FBM
    anc_FBM <- initialize_FBM(FBM_pref, length(samples))

    ## Iterate through VCF
    chrom <- pos <- rsid <- ref <- alt <- anc_snp <- c()
    chunk <- 1
    while(nrow(vcf <- readVcf(tab))){
        print(paste0("Processing chunk ", chunk))

        ## Add SNP info
        chrom <- c(chrom, as.vector(seqnames(vcf)))
        pos <- c(pos, as.vector(ranges(vcf)@start))
        ref <- c(ref, as.vector(ref(vcf)))
        alt <- c(alt, as.vector(unlist(alt(vcf))))
        rsid <- c(rsid, rownames(vcf))
        anc_snp <- c(anc_snp, rep(c(anc, "all"), length.out = (length(anc)+1)*nrow(vcf)))

        ## Get ancestry-genotype matrix
        anc_gt <- get_anc_gt(vcf, anc)

        ## filter minAC
        idx_remove <- colSums(anc_gt) < minAC
        anc_gt <- anc_gt[, !idx_remove]
        chrom <- chrom[!idx_remove]
        pos <- pos[!idx_remove]
        ref <- ref[!idx_remove]
        alt <- alt[!idx_remove]
        rsid <- rsid[!idx_remove]
        anc_snp <- anc_snp[!idx_remove]

        ## Make ancestry-genotye raw type
        nc <- ncol(anc_gt)
        anc_gt <- round(100*anc_gt) %>% as.raw() %>% matrix(., ncol = nc)

        ## Add chunk to FBM
        FBM_nc <- ncol(anc_FBM)
        anc_FBM$add_columns(nc)
        anc_FBM[,(FBM_nc+1):ncol(anc_FBM)] <- anc_gt
        chunk <- chunk+1
    }
    close(tab)
    
    ## Save FBM and SNP names to RDS file
    rds <- sub_bk(anc_FBM$backingfile, ".rds")
    rds_obj <- list(chrom = chrom,
                    pos = pos,
                    ref = ref,
                    alt = alt,
                    rsid = rsid,
                    anc = anc_snp,
                    samples = samples,
                    geno = anc_FBM)
    saveRDS(object = rds_obj, file = rds)
}
