mat_substr_to_int <- function(mat, start, end){
    mat <- apply(mat, 2, function(x) {
        x <- substr(x,start,end)
        mode(x) <- "integer"
        return(x)})
    return(mat)
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

get_anc_gt <- function(gt1, gt2, anc, rf_mat, tract_idx){
    anc1 <- rf_mat[tract_idx, seq(from=1, to=ncol(rf_mat), by=2)]
    anc2 <- rf_mat[tract_idx, seq(from=2, to=ncol(rf_mat), by=2)]    
    anc_gt <- lapply(anc, function(a) {
        (gt1 * (anc1 == a)) + (gt2 * (anc2 == a))})
    anc_gt %<>% append(list(gt1 + gt2)) %>% interleave_matrices()
    anc_gt <- t(anc_gt)
    return(anc_gt)
}

 
##' Make FBM containing ancestry-genotypes and return list with containing FBM and SNP info
##'
##' @param vcf_file a file path to a VCF file
##' @param rf_file a file path to an RFMix output file
##' @param FBM_pref a file path to store FBM. Omit the file extension
##' @param chunk_size an integer indicating the max number of VCF records to read at a time
##' @param minAC an integer indicating the minimum allele count (per-ancestry) to retain
##' @author Frank Ockerman
##' @importFrom VariantAnnotation VcfFile scanVcfHeader ScanVcfParam samples readVcf ref alt
##' @importFrom data.table fread
##' @importFrom SummarizedExperiment assays
##' @importFrom GenomeInfoDb seqnames
##' @importFrom IRanges ranges
##' @importFrom bigstatsr FBM.code256 sub_bk
##' @importFrom magrittr `%>%` `%<>%` 
##' @export
make_rf_FBM <- function(vcf_file, rf_file, FBM_pref, chunk_size, rds=NULL, minAC=1){
    ## Read ancestry tracts
    rf <- fread(file = rf_file, sep = '\t', header = T)
    rf_mat <- as.matrix(rf[,4:ncol(rf)])
    
    ## Reference and open VCF file
    tab <- VcfFile(vcf_file, yieldSize = chunk_size)
    open(tab)
    param <- ScanVcfParam(c("ALT"), geno="HDS")

    ## Extract VCF info
    rf_header <- readLines(gzfile(rf_file), n=1) %>% strsplit(., " +") %>% `[[`(1)
    rf_header <- rf_header[grepl("=", rf_header)]
    rf_codes <- lapply(rf_header, function(x) strsplit(x, "=")[[1]])
    anc_names <- sapply(rf_codes, function(x) x[1])
    anc <- as.integer(sapply(rf_codes, function(x) x[2]))
    samples <- samples(scanVcfHeader(tab))

    ## Initialize FBM
    if(is.null(rds)){
        anc_FBM <- initialize_FBM(FBM_pref, length(samples))
        h_obj <- list(chrom = c(),
                        pos = c(),
                        ref = c(),
                        alt = c(),
                        rsid = c(),
                        anc = c(),
                        anc_ref = c(),
                        samples = c(),
                        geno = anc_FBM)                
    } else {
        h_obj <- readRDS(rds)
        }


    ## Iterate through VCF
    chunk <- 1
    while(nrow(vcf <- readVcf(tab, param=param))){
        print(paste0("Processing chunk ", chunk))

        ## Remove variants not in tracts file
        tmp_pos <- as.vector(ranges(vcf)@start)
        tmp_chm <- as.vector(seqnames(vcf))
        keep <- sapply(1:length(tmp_pos), function(i) {
            any(tmp_pos[i] <= rf$epos & tmp_pos[i] >= rf$spos & tmp_chm[i] == rf$`#chm`)
        })
        if(sum(keep) < 2){chunk <- chunk+1; next}        
        vcf <- vcf[keep]

        ## Remove variants below minAC threshold
        gt1 <- assays(vcf)$HDS[,,1]
        gt2 <- assays(vcf)$HDS[,,2]        
        idx_remove <- rowSums(gt1+gt2) < minAC
        gt1 <- gt1[!idx_remove,]
        gt2 <- gt2[!idx_remove,]        
        vcf <- vcf[!idx_remove]
        if(nrow(vcf) < 2){chunk <- chunk+1; next}

        ## Add SNP info        
        ref <- rep(as.vector(ref(vcf)), each=length(anc)+1)
        alt <- rep(as.vector(unlist(alt(vcf))), each=length(anc)+1)
        rsid <- rep(rownames(vcf), each=length(anc)+1)
        anc_snp <- c(rep(c(anc, "all"), length.out = (length(anc)+1)*nrow(vcf)))

        ## Create ancestry matrix
        chrom <- rep(as.vector(seqnames(vcf)))
        pos <- rep(as.vector(ranges(vcf)@start))        
        tract_idx <- sapply(1:nrow(vcf), function(i) {
            which(rf$`#chm` == chrom[i] & rf$spos <= pos[i] & rf$epos >= pos[i])
        })
        chrom <- rep(as.vector(seqnames(vcf)), each=length(anc)+1)
        pos <- rep(as.vector(ranges(vcf)@start), each=length(anc)+1)        

        ## Get ancestry-genotype matrix
        anc_gt <- get_anc_gt(gt1, gt2, anc, rf_mat, tract_idx)

        ## Assign reference anc per-SNP
        gt_sum <- colSums(anc_gt)
        n1 <- length(anc)+1
        anc_ref <- rep(FALSE, length(gt_sum))
        anc_max <- sapply(1:length(rsid), function(i){
            idx.ref <- which.max(gt_sum[((i-1)*n1+1):(i*n1-1)]) + (i-1)*n1
            anc_ref[idx.ref] <<- TRUE
        })

        ## Make ancestry-genotye raw type
        nc <- ncol(anc_gt)
        anc_gt <- round(100*anc_gt) %>% as.raw() %>% matrix(., ncol = nc)

        ## Add chunk to FBM
        FBM_nc <- ncol(h_obj$geno)
        h_obj$geno$add_columns(nc)
        h_obj$geno[,(FBM_nc+1):ncol(h_obj$geno)] <- anc_gt
        h_obj$chrom <- c(h_obj$chrom, chrom)
        h_obj$pos <- c(h_obj$pos, pos)
        h_obj$ref <- c(h_obj$ref, ref)
        h_obj$alt <- c(h_obj$alt, alt)
        h_obj$rsid <- c(h_obj$rsid, rsid)
        h_obj$anc <- c(h_obj$anc, anc_snp)
        h_obj$anc_ref <- c(h_obj$anc_ref, anc_ref)
        chunk <- chunk+1
    }
    close(tab)
    return(h_obj)
}
