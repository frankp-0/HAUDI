##' Convert FLARE output to RFMix output
##'
##' @title make_Tracts
##' @param anc_vcf_file Path to a VCF file with phased genotypes and estimated local ancestry, i.e. the output of flare (https://doi.org/10.1101/2022.08.02.502540)
##' @param chunk_size Integer indicating the max number of VCF records to read at a time
##' @param out_file File path for storing output
##' @author Frank Ockerman
##' @importFrom VariantAnnotation VcfFile scanVcfHeader samples readVcf ref alt
##' @importFrom SummarizedExperiment assays
##' @importFrom GenomeInfoDb seqnames
##' @importFrom IRanges ranges
##' @import data.table
##' @export
make_Tracts <- function(anc_vcf_file, chunk_size, out_file){
    ## open FLARE VCF
    tab <- VcfFile(anc_vcf_file, yieldSize = chunk_size)
    open(tab)    

    ## write RFMix header
    anc <- unlist(scanVcfHeader(tab)@header@listData$ANCESTRY)
    header <- paste0(
        "#Subpopulation order/codes: ",
        paste(c(paste(names(anc), anc, sep="=")), collapse="    ")
    )
    write(header, file=out_file, append=F)

    ## write column names
    samples <- samples(scanVcfHeader(tab))
    samp_cols <- paste(rep(samples, each = 2), 0:1, sep=".")
    info_cols <- c("#chm", "spos", "epos")
    cols <- c(info_cols, samp_cols)
    write(paste(cols, collapse='\t'), file=out_file, append=T, sep='\t')

    last <- NULL
    chunk <- 1
    ## iterate through FLARE VCF in chunks
    while(nrow(vcf_chunk <- readVcf(tab))){
        print(paste0("Processing chunk: ", chunk, " (", chunk_size*(chunk-1), " SNPs processed)"))

        ## iterate through chromosomes in current chunk
        chrom_all <- as.vector(seqnames(vcf_chunk))
        for (chrom in unique(chrom_all)){
            vcf <- vcf_chunk[chrom_all == chrom]

            ## get matrix of ancestries
            mat <- t(interleave_matrices(list(
                t(assays(vcf)$AN1),
                t(assays(vcf)$AN2)
            )))

            ## start from last variant at previous chunk
            pos <- c(as.vector(ranges(vcf)@start))            
            if(!is.null(last)){
                if(last$chrom == chrom){
                    pos <- c(last$pos, pos)
                    mat <- rbind(last$mat, mat)
                }
            }

            ## switch points (ancestry changes for > 0 individuals)
            switch <- which(c(FALSE, sapply(2:nrow(mat), function(i) {!all(mat[i,] == mat[i-1,])})))

            ## take midpoint between positions when ancestry changes
            spos <- rep(0, nrow(mat))            
            for (i in 1:(nrow(mat)-1)){
                if (i %in% switch){
                    spos[i] <- ceiling(mean(c(pos[i], pos[i-1])))
                } else{
                    spos[i] <- pos[i]
                }
            }

            ## only keep switchpoints
            mat <- mat[switch,]
            
            ## record start pos and end pos
            spos <- spos[switch]
            epos <- spos[2:length(spos)]-1
            spos <- spos[1:(length(spos)-1)]

            ## create RFMix table and write result
            dt <- data.table(
                `#chm` = chrom,
                spos = spos,
                epos = epos)
            dt_anc <- data.table(mat[1:(nrow(mat)-1),])
            colnames(dt_anc) <- samp_cols
            dt <- cbind(dt, dt_anc)
            fwrite(x=dt, file=out_file, append=T, sep='\t', quote=F, row.names=F, col.names=F)

            ## record last variant from chunk
            last <- list(chrom=chrom, pos=pos[length(pos)], mat=mat[nrow(mat),])        
        }
        chunk <- chunk+1
    }
}
