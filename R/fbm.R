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
#' @noRd
resolve_indices <- function(ids = NULL, idx = NULL,
                            reference, label = "items") {
  if (!is.null(ids) && !is.null(idx)) {
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
#' @noRd
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
#' @param pgen A `pgen` object produced by `pgenlibr::NewPgen()`
#' @param pvar A data frame corresponding to the pvar file
#' @param tracts A list of data.frames, where each element corresponds to a
#' sample and has columns `index` (representing indices
#' in a plink2 pvar file) and "anc0" and "anc1", corresponding to ancestries for
#' each haplotype
#' @param min_ac An integer with the minimum allele count to
#' retain columns in the HAUDI matrix
#' @param anc_names An ordered vector of ancestry names
#' @param idx_samples An ordered integer vector of samples to
#' retain from the pgen
#' @return A list with the HAUDI matrix and variant info for the chunk
#' @noRd
make_haudi_chunk <- function(chunk, pgen, pvar, tracts,
                             min_ac, anc_names, idx_samples) {
  ## Get ancestry matrices
  time_anc <- system.time({
    anc_mats <- query_tracts(chunk, tracts)
  })
  message(sprintf("Anc query time: %s", time_anc[3]))

  ## Per-haplotype alleles
  n_samp <- length(idx_samples)
  time_gen <- system.time({
    gen0 <- gen1 <- matrix(0, n_samp, length(chunk))
    buf <- pgenlibr::AlleleCodeBuf(pgen)
    for (i in seq_along(chunk)) {
      pgenlibr::ReadAlleles(pgen, acbuf = buf, variant_num = chunk[i])
      gen0[, i] <- buf[1, idx_samples]
      gen1[, i] <- buf[2, idx_samples]
    }
  })
  message(sprintf("Gen query time: %s", time_gen[3]))

  ## Per-ancestry genotypes and total genotypes
  time_haudi <- system.time({
    mat_haudi <- lapply(0:(length(anc_names) - 1), function(anc) {
      gen0 * (anc_mats$hap0 == anc) + gen1 * (anc_mats$hap1 == anc)
    }) |>
      do.call(what = "cbind") |>
      cbind(gen0 + gen1)
  })
  message(sprintf("HAUDI time: %s", time_haudi[3]))

  time_process <- system.time({
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

    ## Set reference ancestry
    dt_info[, row_id := .I]
    dt_info[, idx := rep(chunk, length(anc_names) + 1)]
    max_rows <- dt_info[anc != "all", .SD[which.max(ac)], by = idx][, row_id]
    dt_info[, anc_ref := FALSE]
    dt_info[max_rows, anc_ref := TRUE]
    dt_info$row_id <- dt_info$idx <- NULL


    ## Filter by ancestry-specific allele count
    mat_haudi <- mat_haudi[, dt_info$ac >= min_ac]
    dt_info <- dt_info[ac >= min_ac, ]
  })
  message(sprintf("Process time: %s", time_process[3]))

  return(list(mat = mat_haudi, info = dt_info))
}


#' Append a new chromosome to the FBM
#'
#' @inheritParams make_fbm
#' @param lanc_file A string with the file path for a single local
#' ancestry input file
#' @param plink_prefix A string with the prefix for a single set of plink2
#' files
#' @param idx_variants An integer vector with indices of variants to include
#' in the FBM
#' @inherit make_fbm return
#' @export
add_to_fbm <- function(lanc_file, plink_prefix,
                       fbm_prefix, variants = NULL, idx_variants = NULL,
                       min_ac = 0, samples = NULL, idx_samples = NULL,
                       anc_names = NULL, chunk_size = 400, fbm = NULL) {
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
  message(sprintf(
    "[%s] Reading ancestry tracts: %s",
    format(Sys.time(), "%H:%M:%S"), lanc_file
  ))
  tracts <- read_lanc(lanc_file)
  message(sprintf(
    "[%s] Finished reading ancestry tracts",
    format(Sys.time(), "%H:%M:%S")
  ))

  ## Subset samples matching with tracts
  tracts <- tracts[idx_samples]

  ## Get ordered ancestries
  if (is.null(anc_names)) {
    anc_names <- sapply(tracts, function(x) c(x$anc0, x$anc1)) |>
      unlist() |>
      unique() |>
      sort() |>
      as.character()
  }

  ## Initialize FBM
  if (is.null(fbm)) {
    fbm <- initialize_fbm(fbm_prefix, psam$`#IID`[idx_samples])
  }

  ## Initialize info
  dt_info <- data.table::data.table()

  ## Prepare chunks
  chunks <- split(
    idx_variants,
    ceiling(seq_len(length(idx_variants)) / chunk_size)
  )
  n_chunks <- length(chunks)
  start_time <- Sys.time()

  ## Initialize progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] Chunk :current/:total — ETA: :eta — :message",
    total = n_chunks,
    clear = FALSE, width = 70
  )
  pb$tick(0)

  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    chunk_start_time <- Sys.time()
    pb$message(sprintf("Processing chunk %i: %d variants", i, length(chunk)))

    ## Get new HAUDI-style matrix and associated info
    haudi_chunk <- make_haudi_chunk(
      chunk = chunk, pgen = pgen, pvar = pvar, tracts = tracts, min_ac = min_ac,
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

    elapsed <- difftime(Sys.time(), chunk_start_time, units = "secs")
    pb$message(sprintf(
      "Finished (%d/%d columns kept, %.1f sec)",
      nc_chunk, length(chunk) * (length(anc_names) + 1), elapsed
    ))
    pb$tick()
  }

  total_time <- difftime(Sys.time(), start_time, units = "mins")
  message(sprintf(
    "\n[%s] All chunks complete in %.1f minutes.",
    format(Sys.time(), "%H:%M:%S"), total_time
  ))

  return(list(fbm = fbm, info = dt_info))
}


#' Make File-Backed Matrix input for HAUDI
#'
#' @description
#' Processes .lanc local ancestry files and corresponding plink2 files,
#' producing a file-backed matrix object and associated
#' information used for HAUDI.
#'
#' @details
#' Only one (or neither) of `variants` and `idx_variants_list` may be provided.
#' `variants` is a single character with all variant IDs (across input files),
#' while `idx_variants_list` is a list where each element is a vector of indices
#' corresponding to a single set of plink2 files (e.g. each vector is
#' per-chromosome). Only one (or neither) of `samples` and `idx_samples` may
#' be provided. Exactly one of `fbm_prefix` or `fbm` may be provided.
#' It is assumed that the local ancestry input files and plink2 inputs
#' are sorted by chromosome and position.
#'
#' @param lanc_files A string vector with file paths for
#' the local ancestry input
#' @param plink_prefixes A string vector with the prefixes for plink2 file paths
#' @param fbm_prefix A string with the prefix for the
#' file path where a new backing file for the FBM
#' will be written. This prefix will be appended with ".bk"
#' @param variants A character vector with variant IDs to include in the FBM
#' @param idx_variants_list A list of integer vectors, one per plink2 input
#' file, containing indices of variants to include in the FBM.
#' @param min_ac An integer with a minimum cut off for
#' allele counts in each column of the FBM
#' @param samples A character vector of samples to
#' include in the FBM. It is assumed all plink2 files
#' have the same set of samples.
#' @param idx_samples An integer vector of sample indices,
#' corresponding to plink2 psam input, to include in the FBM.
#' It is assumed all plink2 files have the same set of samples.
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
#'
#' @importFrom pgenlibr NewPgen ReadAlleles AlleleCodeBuf
#' @importFrom data.table data.table setorder fread
#' @importFrom bigstatsr FBM.code256
#' @importFrom progress progress_bar
#' @export
make_fbm <- function(lanc_files, plink_prefixes,
                     fbm_prefix, variants = NULL, idx_variants_list = NULL,
                     min_ac = 0, samples = NULL, idx_samples = NULL,
                     anc_names = NULL, chunk_size = 400) {
  result <- add_to_fbm(
    lanc_files[1], plink_prefixes[1],
    fbm_prefix, variants, idx_variants_list[[1]],
    min_ac, samples, idx_samples, anc_names, chunk_size
  )
  dt_info <- result$info
  if (length(lanc_files) > 1) {
    for (i in 2:length(lanc_files)) {
      result <- add_to_fbm(
        lanc_file = lanc_files[i],
        plink_prefix = plink_prefixes[i], fbm_prefix = NULL,
        variants = variants, idx_variants = idx_variants_list[[i]],
        min_ac = min_ac, samples = samples,
        idx_samples = idx_samples, anc_names = anc_names,
        chunk_size = chunk_size, fbm = result$fbm
      )
      dt_info <- rbind(dt_info, result$info)
    }
  }
  return(list(fbm = result$fbm, info = dt_info))
}
