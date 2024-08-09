## make_fbm.R -- functions for generating HAUDI-style File-Backed Matrices

## Combines list of matrices into single matrix,
## by interleaving their rows
interleave_matrices <- function(list_of_matrices) {
  idx_interleaved <- order(sequence(sapply(list_of_matrices, nrow)))
  val <- do.call(rbind, list_of_matrices)[idx_interleaved, ]
  return(val)
}


## Initializes a File-Backed Matrix (code 256)
initialize_fbm <- function(fbm_pref, nrow) {
  code_dosage <- rep(NA_real_, 256)
  code_dosage[1:201] <- seq(0, 2, length.out = 201)
  backing_file <- paste0(fbm_pref, ".bk")
  if (file.exists(backing_file)) {
    file.remove(backing_file)
  }
  anc_fbm <- FBM.code256(
    nrow = nrow,
    ncol = 0,
    backingfile = fbm_pref,
    code = code_dosage
  )
  return(anc_fbm)
}


interpolate_ancestry <- function(x, x_prev = NULL, pos, pos_prev = NULL) {
  if (all(is.na(x)) && is.null(x_prev)) {
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
      x[seq_along(x)] <- x_prev
      return(x)
    }
  } else {
    x <- c(x_prev, x)
    pos <- c(pos_prev, pos)

    idx_na <- which(is.na(x))
    idx_no_na <- which(!is.na(x))

    lower <- findInterval(idx_na, idx_no_na, all.inside = TRUE)
    upper <- pmin(lower + 1, length(idx_no_na))
    dist_lower <- abs(pos[idx_na] - pos[idx_no_na[lower]])
    dist_upper <- abs(pos[idx_no_na[upper]] - pos[idx_na])

    x[idx_na] <- ifelse(dist_lower < dist_upper,
      x[idx_no_na[lower]],
      x[idx_no_na[upper]]
    )
    return(x[(1 + length(pos_prev)):length(x)])
  }
}

##' Make HAUDI-style FBM and return list with FBM and SNP info
##'
##' @param vcf_file File path to a VCF file
##' @param fbm_pref File path to store FBM. Omit the file extension
##' @param chunk_size The maximum number of VCF records to read at a time
##' @param rds Optional file path for an existing RDS  file to append data to
##' @param min_ac Minimum allele count (per-ancestry) to retain
##' @param geno_format The VCF format field to use ("HDS" or "GT")
##' @param anc_names vector of labels for each population
##' @importFrom data.table fread
##' @importFrom SummarizedExperiment assays
##' @importFrom GenomeInfoDb seqnames
##' @importFrom IRanges ranges
##' @importFrom bigstatsr FBM.code256 sub_bk
##' @export
make_fbm <- function(vcf_file, fbm_pref, chunk_size, rds = NULL,
                     min_ac = 1, geno_format = "HDS", anc_names) {
  tab <- VcfFile(vcf_file, yieldSize = chunk_size)
  open(tab)
  param <- ScanVcfParam(c("ALT"), geno = c(geno_format, "AN1", "AN2"))
  ## Initialize FBM
  anc_fbm <- initialize_fbm(fbm_pref, length(samples(scanVcfHeader(tab))))
  fbm_info <- list(
    chrom = c(),
    pos = c(),
    ref = c(),
    alt = c(),
    rsid = c(),
    anc = c(),
    anc_ref = c()
  )

  ## Iterate through VCF
  chunk <- 1
  an1_prev <- an2_prev <- pos_prev <- NULL

  while (nrow(vcf <- readVcf(tab, param = param))) {
    print(paste0("Processing chunk ", chunk))

    ## Remove variants below min_ac threshold
    if (geno_format == "HDS") {
      gt1 <- assays(vcf)$HDS[, , 1]
      gt2 <- assays(vcf)$HDS[, , 2]
    } else if (geno_format == "GT") {
      gt1 <- assays(vcf)$GT
      gt2 <- apply(substr(gt1, 3, 3), 2, FUN = as.numeric)
      gt1 <- apply(substr(gt1, 1, 1), 2, FUN = as.numeric)
    } else {
      stop("Please specify either HDS or GT for geno_format")
    }
    idx_remove <- rowSums(gt1 + gt2) < min_ac
    gt1 <- gt1[!idx_remove, ]
    gt2 <- gt2[!idx_remove, ]
    vcf <- vcf[!idx_remove]

    if (nrow(vcf) < 2) {
      chunk <- chunk + 1
      next
    }

    ## Add SNP info
    chrom <- rep(as.vector(seqnames(vcf)), each = length(anc_names) + 1)
    ref <- rep(as.vector(ref(vcf)), each = length(anc_names) + 1)
    alt <- rep(as.vector(unlist(alt(vcf))), each = length(anc_names) + 1)
    rsid <- rep(rownames(vcf), each = length(anc_names) + 1)
    anc_snp <- c(rep(c(anc_names, "all"),
      length.out = (length(anc_names) + 1) * nrow(vcf)
    ))

    ## Get ancestries
    an1 <- assays(vcf)$AN1
    an2 <- assays(vcf)$AN2

    ## Interpolate ancestries
    pos <- as.vector(ranges(vcf)@start)
    an1 <- sapply(seq_len(ncol(an1)), function(j) {
      interpolate_ancestry(an1[, j], an1_prev[, j], pos, pos_prev)
    })
    an2 <- sapply(seq_len(ncol(an2)), function(j) {
      interpolate_ancestry(an2[, j], an2_prev[, j], pos, pos_prev)
    })

    an1_prev <- an1[nrow(an1), , drop = FALSE]
    an2_prev <- an2[nrow(an2), , drop = FALSE]
    pos_prev <- pos[length(pos)]

    ## Get ancestry-genotype matrix
    anc_gt <- lapply(0:(length(anc_names) - 1), function(anc) {
      gt1 * (an1 == anc) + gt2 * (an2 == anc)
    })
    anc_gt <- append(x = anc_gt, values = list(gt1 + gt2)) |>
      interleave_matrices() |>
      t()

    ## Assign reference anc per-SNP
    gt_sum <- colSums(anc_gt)
    n1 <- length(anc_names) + 1
    anc_ref <- rep(FALSE, length(gt_sum))

    for (i in seq_len(nrow(vcf))) {
      idx_ref <- gt_sum[((i - 1) * n1 + 1):(i * n1 - 1)] |>
        which.max() + (i - 1) * n1
      anc_ref[idx_ref] <- TRUE
    }

    ## Make ancestry-genotye raw type
    nc <- ncol(anc_gt)
    anc_gt <- round(100 * anc_gt) |>
      as.raw() |>
      matrix(ncol = nc)

    ## Add chunk to FBM
    fbm_nc <- ncol(anc_fbm)
    anc_fbm$add_columns(nc)
    anc_fbm[, (fbm_nc + 1):ncol(anc_fbm)] <- anc_gt
    fbm_info$chrom <- c(fbm_info$chrom, chrom)
    fbm_info$pos <- c(fbm_info$pos, pos)
    fbm_info$ref <- c(fbm_info$ref, ref)
    fbm_info$alt <- c(fbm_info$alt, alt)
    fbm_info$rsid <- c(fbm_info$rsid, rsid)
    fbm_info$anc <- c(fbm_info$anc, anc_snp)
    fbm_info$anc_ref <- c(fbm_info$anc_ref, anc_ref)
    chunk <- chunk + 1
  }
  close(tab)
  return(list(
    FBM = anc_fbm,
    info = data.frame(fbm_info)
  ))
}
