#' Match indices from (optional) subset and reference
#'
#' This utility function generates a vector of indices
#' corresponding to `reference`. When `ids`
#' or `idx` are provided, these are used
#' to subset the result.
#'
#' `ids` and `idx` may not both be provided
#'
#' @param ids A character vector which contains
#' a subset of `reference`
#' @param idx An integer vector of indices
#' representing a subset of `reference`
#' @param reference A vector of values
#' which we want indices from
#' @param label A string with what the `reference`
#' represents (intended to be "samples" or "variants")
#' @return A vector of sorted indices representing a
#' (possible) subset of `reference`
resolve_indices <- function(ids = NULL, idx = NULL,
                            reference, label = "items") {
  if (!is.null(idx) && !is.null(idx)) {
    stop(paste0(
      "Only one of `", label, "` or `idx_",
      label, "` may be provided"
    ))
  }
  if (is.null(ids) && is.null(idx)) {
    return(seq_len(length(reference)))
  }
  if (!is.null(ids)) {
    idx <- match(ids, reference) |>
      na.omit() |>
      as.vector()
  }
  if (!all(idx %in% seq_len(length(reference)))) {
    stop("Some `", label, "` indices exist outside range")
  }
  if (length(idx) == 0) {
    stop(paste0("No matches found with `", label, "`"))
  }
  return(sort(idx))
}

#' Check that plink2 files exist
#'
#' @param plink_prefix A string with the prefix
#' for the plink2 file paths
#' @return A list of plink2 file paths
verify_plink <- function(plink_prefix) {
  files <- paste0(plink_prefix, c(".pgen", ".pvar", ".psam")) |> as.list()
  for (file in files) {
    if (!file.exists(file)) {
      stop(paste0(file, " does not exist"))
    }
  }
  names(files) <- c("pgen", "pvar", "psam")
  return(files)
}

#' Creates a new File-Backed matrix for HAUDI
#'
#' This utility function creates an FBM.256 object
#' with the correct coding to represent genotype
#' values from 0 to 2.
#'
#' @param fbm_prefix A string with the prefix for
#' the intended FBM file path
#' @return An FBM object with an attribute "samples"
#' containing an ordered vector of samples
initialize_fbm <- function(fbm_prefix, samples) {
  code_dosage <- rep(NA_real_, 256)
  code_dosage[1:201] <- seq(0, 2, length.out = 201)
  backing_file <- paste0(fbm_prefix, ".bk")
  if (file.exists(backing_file)) {
    file.remove(backing_file)
  }
  fbm <- bigstatsr::FBM.code256(
    nrow = length(samples),
    ncol = 0,
    backingfile = fbm_prefix,
    code = code_dosage
  )
  attr(fbm, "samples") <- samples
  return(fbm)
}

#' Produces HAUDI matrix and variant info for a single chunk
#'
#' @param chunk An integer vector of indices in the pgen file
#' @param gr_tracts A `GRanges` object corresponding to the ancestry file
#' @param pgen A `pgen` object produced by `pgenlibr::NewPgen()`
#' @param pvar A data frame corresponding to the pvar file
#' @param min_ac An integer with the minimum allele count to
#' retain columns in the HAUDI matrix
#' @param anc_names An ordered vector of ancestry names
#' @param idx_samples An ordered integer vector of samples to
#' retain from the pgen
#' @return A list with the HAUDI matrix and variant info for the chunk
make_haudi_chunk <- function(chunk, gr_tracts, pgen, pvar,
                             min_ac, anc_names, idx_samples) {
  n_samp <- length(idx_samples)

  gr_target <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(values = pvar$`#CHROM`[chunk], lengths = 1),
    ranges = IRanges::IRanges(pvar$POS[chunk], pvar$POS[chunk]),
  )

  ## Find matching tract per-sample/hap for each variant in gr_target
  overlaps <- GenomicRanges::findOverlaps(gr_target, gr_tracts)
  gr_chunk <- gr_tracts[S4Vectors::subjectHits(overlaps)]

  ## Per-haplotype ancestry
  anc0 <- matrix(gr_chunk[gr_chunk$hap == 0]$ancestry, nrow = n_samp)
  anc1 <- matrix(gr_chunk[gr_chunk$hap == 1]$ancestry, nrow = n_samp)

  ## Per-haplotype alleles
  gen0 <- gen1 <- matrix(0, n_samp, length(chunk))
  buf <- pgenlibr::AlleleCodeBuf(pgen)
  for (i in seq_along(chunk)) {
    pgenlibr::ReadAlleles(pgen, acbuf = buf, variant_num = chunk[i])
    gen0[, i] <- buf[1, idx_samples]
    gen1[, i] <- buf[2, idx_samples]
  }

  ## Per-ancestry genotypes and total genotypes
  mat_haudi <- lapply(0:(length(anc_names) - 1), function(anc) {
    gen0 * (anc0 == anc) + gen1 * (anc1 == anc)
  }) |>
    do.call(what = "cbind") |>
    cbind(gen0 + gen1)

  dt_info <- data.table::data.table(
    chrom = rep(pvar[chunk, ]$`#CHROM`, length(anc_names) + 1),
    pos = rep(pvar[chunk, ]$POS, length(anc_names) + 1),
    id = rep(pvar[chunk, ]$ID, length(anc_names) + 1),
    ref = rep(pvar[chunk, ]$REF, length(anc_names) + 1),
    alt = rep(pvar[chunk, ]$ALT, length(anc_names) + 1),
    anc = c(rep(anc_names, each = length(chunk)), rep("all", length(chunk)))
  )
  dt_info$chrom <- as.character(dt_info$chrom)
  dt_info$ac <- colSums(mat_haudi)

  ## Filter by ancestry-specific allele count
  mat_haudi <- mat_haudi[, dt_info$ac >= min_ac]
  dt_info <- dt_info[ac >= min_ac, ]

  return(list(mat = mat_haudi, info = dt_info))
}


#' Append a new chromosome to the FBM
#'
#' @inheritParams make_fbm
#' @param ancestry_file A string with the file path for a single local
#' ancestry input file
#' @param plink_prefix A string with the prefix for a single set of plink2
#' files
#' @inherit make_fbm return
#' @export
add_to_fbm <- function(ancestry_file, ancestry_fmt, plink_prefix,
                       fbm_prefix, variants = NULL, idx_variants = NULL,
                       min_ac = 0, samples = NULL, idx_samples = NULL,
                       anc_names = NULL, chunk_size = 400, fbm = NULL) {
  ## TODO: add helpful output on how many samples, variants match
  ## TODO: check if GRanges or just data.table is best for matching regions

  ## Verify and read plink2 input files
  plink_files <- verify_plink(plink_prefix)
  pgen <- pgenlibr::NewPgen(plink_files$pgen)
  pvar <- data.table::fread(plink_files$pvar, skip = "#CHROM")
  psam <- data.table::fread(plink_files$psam, skip = "#IID")

  ## Get variant indices in pvar
  idx_variants <- resolve_indices(
    ids = variants, idx = idx_variants,
    reference = pvar$ID, label = "variants"
  )

  ## Get sample indices in psam
  idx_samples <- resolve_indices(
    ids = samples, idx = idx_samples,
    reference = psam$`#IID`, label = "samples"
  )

  ## Read ancestry tracts
  print(sprintf("Reading ancestry tracts for %s", ancestry_file))
  dt_tracts <- read_ancestry_tracts(
    file = ancestry_file, file_fmt = ancestry_fmt,
    extend_tracts = TRUE, plink_prefix = plink_prefix
  )
  print(sprintf("Finished reading ancestry tracts"))

  ## Subset samples matching with dt_tracts
  dt_tracts <- dt_tracts[sample %in% psam$`#IID`[idx_samples], ]
  dt_tracts$sample <- factor(dt_tracts$sample,
    levels = psam$`#IID`[idx_samples]
  )
  idx_samples <- idx_samples[
    which(psam$`#IID`[idx_samples] %in% dt_tracts$sample)
  ]
  data.table::setorder(dt_tracts, by = sample, hap, chrom, spos)

  ## Get ordered ancestries
  if (is.null(anc_names)) {
    anc_names <- unique(dt_tracts$ancestry)
  }

  ## Make GRanges object for ancestry tracts
  gr_tracts <- GenomicRanges::makeGRangesFromDataFrame(
    df = dt_tracts,
    seqnames.field = "chrom",
    start.field = "spos",
    end.field = "epos",
    keep.extra.columns = TRUE
  )

  ## Initialize FBM
  if (is.null(fbm)) {
    fbm <- initialize_fbm(fbm_prefix, psam$`#IID`[idx_samples])
  }

  ## Initialize info
  dt_info <- data.table::data.table()

  ## Loop through chunks
  chunks <- split(
    idx_variants,
    ceiling(seq_len(length(idx_variants)) / chunk_size)
  )
  i <- 1
  n_chunks <- length(chunks)
  for (chunk in chunks) {
    ## Get new HAUDI-style matrix and associated info
    haudi_chunk <- make_haudi_chunk(
      chunk = chunk, gr_tracts = gr_tracts,
      pgen = pgen, pvar = pvar, min_ac = min_ac,
      anc_names = anc_names, idx_samples = idx_samples
    )

    ## Append matrix to FBM
    nc_chunk <- ncol(haudi_chunk$mat)
    nc_fbm <- ncol(fbm)
    mat_haudi_raw <- (100 * haudi_chunk$mat) |>
      round() |>
      as.raw() |>
      matrix(ncol = nc_chunk)
    fbm$add_columns(nc_chunk)
    fbm[, (nc_fbm + 1):(nc_fbm + nc_chunk)] <- mat_haudi_raw

    ## Append info
    dt_info <- rbind(dt_info, haudi_chunk$info)
    print(sprintf("finished chunk %s out of %s", i, n_chunks))
    i <- i + 1
  }
  return(list(fbm = fbm, info = dt_info))
}


#' Make File-Backed Matrix input for HAUDI
#'
#' @description
#' Processes local ancestry files and corresponding plink2 files,
#' producing a file-backed matrix object and associated
#' information used for HAUDI.
#'
#' @details
#' Only one (or neither) of `variants` and `idx_variants`
#' may be provided. Only one (or neither) of `samples` and
#' `idx_samples` may be provided. Exactly one of `fbm_prefix`
#' or `fbm` may be provided. It is assumed that
#' the local ancestry input file and plink2 input
#' are sorted by chromosome and position.
#'
#' @param ancestry_files A string vector with file paths for
#' the local ancestry input
#' @param ancestry_fmt A string with the local ancestry format
#' of `ancestry_file`. Either "FLARE", "RFMix", or "lanc"
#' @param plink_prefixes A string vector with the prefixes for plink2 file paths
#' @param fbm_prefix A string with the prefix for the
#' file path where a new backing file for the FBM
#' will be written. This prefix will be appended with ".bk"
#' @param variants A character vector of variant IDs to
#' include in FBM
#' @param variant_idx An integer vector of variant indices,
#' corresponding to plink2 pvar input, to include in FBM
#' @param min_ac An integer with a minimum cut off for
#' allele counts in each column of the FBM
#' @param samples A character vector of samples to
#' include in FBM
#' @param idx_samples An integer vector of sample indices,
#' corresponding to plink2 psam input, to include in FBM
#' @param anc_names A character vector containing (ordered)
#' names for each population in `ancestry_file`. If not
#' provided, integer values are used
#' @param chunk_size An integer with the max number of
#' variants to load at a given time while reading the
#' plink2 pgen file
#' @param fbm An existing `FBM.code256` object to which
#' results are appended.
#'
#' @return A list containing:
#'   `fbm`: an `FBM.code256` object from `bigstatsr`,
#'   which links to a file-backed matrix used for
#'   fittting HAUDI models. Each column contains
#'   either ancestry-specific genotypes (e.g. number of
#'   alternate alleles where ancestry is also "pop0") or
#'   the total genotype. This object also contains a
#'   `samples` attribute ordered according to the samples
#'   in the matrix rows.
#'   `info` a `data.table` containing information
#'   about columns in the FBM.
#' }
#'
#' @importFrom pgenlibr NewPgen ReadAlleles AlleleCodeBuf
#' @importFrom data.table data.table setorder fread
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom bigstatsr FBM.code256
#' @importFrom S4Vectors Rle subjectHits
#' @importFrom IRanges IRanges
#' @export
make_fbm <- function(ancestry_files, ancestry_fmt, plink_prefixes,
                     fbm_prefix, variants = NULL, idx_variants = NULL,
                     min_ac = 0, samples = NULL, idx_samples = NULL,
                     anc_names = NULL, chunk_size = 400) {
  result <- add_to_fbm(
    ancestry_files[1], ancestry_fmt, plink_prefixes[1],
    fbm_prefix, variants, idx_variants, min_ac, samples,
    idx_samples, anc_names, chunk_size
  )
  dt_info <- result$info
  if (length(ancestry_files) > 1) {
    for (i in 2:length(ancestry_files)) {
      result <- add_to_fbm(
        ancestry_file = ancestry_files[i], ancestry_fmt = ancestry_fmt,
        plink_prefix = plink_prefixes[i], fbm_prefix = NULL,
        variants = variants, idx_variants = idx_variants,
        min_ac = min_ac, samples = samples,
        idx_samples = idx_samples, anc_names = anc_names,
        chunk_size = chunk_size, fbm = result$fbm
      )
      dt_info <- rbind(dt_info, result$info)
    }
  }
  return(list(fbm = result$fbm, info = dt_info))
}
