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


interpolate_ancestry <- function(x, x_prev = NULL, pos, pos_prev = NULL) {
    if (all(is.na(x)) & is.null(x_prev)) {
        stop("No local ancestry information in segment")
    } else if (all(!is.na(x))) {
        return(x)
    }

    if (all(!is.na(x))) {
        return(x)
    } else if (all(is.na(x))) {
        if (is.null(x_prev)) {
            stop("No local ancestry information in segment")
        } else {
            x[1:length(x)] <- x_prev
            return(x)
        }
    } else {
        x <- c(x_prev,x)
        pos <- c(pos_prev, pos)
        
        idx_NA <- which(is.na(x))
        idx_no_NA <- which(!is.na(x))

        lower <- findInterval(idx_NA, idx_no_NA, all.inside = T)
        upper <- pmin(lower + 1, length(idx_no_NA))
        dist_lower <- abs(pos[idx_NA] - pos[idx_no_NA[lower]])
        dist_upper <- abs(pos[idx_no_NA[upper]] - pos[idx_NA])

        x[idx_NA] <- ifelse(dist_lower < dist_upper,
                            x[idx_no_NA[lower]],
                            x[idx_no_NA[upper]])
        return(x[(1 + length(pos_prev)):length(x)])
    }
}

##' Make FBM containing ancestry-genotypes and return list containing FBM and SNP info
##'
##' @param vcf_file File path to a VCF file
##' @param FBM_pref File path to store FBM. Omit the file extension
##' @param chunk_size Integer indicating the max number of VCF records to read at a time
##' @param rds Optional file path for an existing RDS  file to append data to
##' @param minAC Integer indicating the minimum allele count (per-ancestry) to retain
##' @param geno_format either "HDS" or "GT", specifying the VCF format field to use
##' @author Frank Ockerman
##' @importFrom VariantAnnotation VcfFile scanVcfHeader ScanVcfParam samples readVcf ref alt
##' @importFrom data.table fread
##' @importFrom SummarizedExperiment assays
##' @importFrom GenomeInfoDb seqnames
##' @importFrom IRanges ranges
##' @importFrom bigstatsr FBM.code256 sub_bk
##' @export
make_FBM <- function(vcf_file, FBM_pref, chunk_size, rds=NULL,
                        minAC=1, geno_format="HDS", anc_names){
    ## Reference and open VCF file
    tab <- VcfFile(vcf_file, yieldSize = chunk_size)
    open(tab)
    param <- ScanVcfParam(c("ALT"), geno=c(geno_format, "AN1", "AN2"))

    ## Initialize FBM
    anc_FBM <- initialize_FBM(FBM_pref, length(samples(scanVcfHeader(tab))))
    FBM_info <- list(chrom = c(),
                     pos = c(),
                     ref = c(),
                     alt = c(),
                     rsid = c(),
                     anc = c(),
                     anc_ref = c())
    
    ## Iterate through VCF
    chunk <- 1
    AN1_prev <- AN2_prev <- pos_prev <- NULL
    
    while(nrow(vcf <- readVcf(tab, param=param))){
        print(paste0("Processing chunk ", chunk))

        ## Remove variants below minAC threshold
        if (geno_format == "HDS") {
            gt1 <- assays(vcf)$HDS[,,1]
            gt2 <- assays(vcf)$HDS[,,2]                    
        } else if (geno_format == "GT") {
            gt1 <- assays(vcf)$GT
            gt2 <- apply(substr(gt1, 3, 3), 2, FUN=as.numeric)
            gt1 <- apply(substr(gt1, 1, 1), 2, FUN=as.numeric)
        } else {
            stop("Please specify either HDS or GT for geno_format")
        }
        idx_remove <- rowSums(gt1+gt2) < minAC
        gt1 <- gt1[!idx_remove,]
        gt2 <- gt2[!idx_remove,]        
        vcf <- vcf[!idx_remove]
        
        if(nrow(vcf) < 2){chunk <- chunk+1; next}

        ## Add SNP info
        chrom <- rep(as.vector(seqnames(vcf)), each=length(anc_names)+1)
        ref <- rep(as.vector(ref(vcf)), each=length(anc_names)+1)
        alt <- rep(as.vector(unlist(alt(vcf))), each=length(anc_names)+1)
        rsid <- rep(rownames(vcf), each=length(anc_names)+1)
        anc_snp <- c(rep(c(anc_names, "all"), length.out = (length(anc_names)+1)*nrow(vcf)))

        ## Get ancestries
        AN1 <- assays(vcf)$AN1
        AN2 <- assays(vcf)$AN2

        ## Interpolate ancestries
        pos <- as.vector(ranges(vcf)@start)
        AN1 <- sapply(1:ncol(AN1), function(j) {
            interpolate_ancestry(AN1[,j], AN1_prev[,j], pos, pos_prev)            
        })
        AN2 <- sapply(1:ncol(AN2), function(j) {
            interpolate_ancestry(AN2[,j], AN2_prev[,j], pos, pos_prev)            
        })        

        AN1_prev <- AN1[nrow(AN1), , drop = F]
        AN2_prev <- AN2[nrow(AN2), , drop = F]
        pos_prev <- pos[length(pos)]

        ## Get ancestry-genotype matrix
        anc_gt <- lapply(0:(length(anc_names)-1), function(anc) {
            gt1 * (AN1 == anc) + gt2 * (AN2 == anc)
        })
        anc_gt <- append(x = anc_gt, values = list(gt1 + gt2)) |>
            interleave_matrices() |> t()

        ## Assign reference anc per-SNP
        gt_sum <- colSums(anc_gt)
        n1 <- length(anc_names)+1
        anc_ref <- rep(FALSE, length(gt_sum))
        anc_max <- sapply(1:nrow(vcf), function(i){
            idx.ref <- which.max(gt_sum[((i-1)*n1+1):(i*n1-1)]) + (i-1)*n1
            anc_ref[idx.ref] <<- TRUE
        })

        ## Make ancestry-genotye raw type
        nc <- ncol(anc_gt)
        anc_gt <- round(100*anc_gt) |> as.raw() |> matrix(ncol = nc)

        ## Add chunk to FBM
        FBM_nc <- ncol(anc_FBM)
        anc_FBM$add_columns(nc)
        anc_FBM[,(FBM_nc+1):ncol(anc_FBM)] <- anc_gt
        FBM_info$chrom <- c(FBM_info$chrom, chrom)
        FBM_info$pos <- c(FBM_info$pos, pos)
        FBM_info$ref <- c(FBM_info$ref, ref)
        FBM_info$alt <- c(FBM_info$alt, alt)
        FBM_info$rsid <- c(FBM_info$rsid, rsid)
        FBM_info$anc <- c(FBM_info$anc, anc_snp)
        FBM_info$anc_ref <- c(FBM_info$anc_ref, anc_ref)
        chunk <- chunk+1
    }
    close(tab)
    return(list(FBM=anc_FBM,
                info=data.frame(FBM_info)))
}
