#' Check that ancestry tracts are valid
#'
#' @param dt_tracts A data.table with ancestry tracts
#' @noRd
validate_tracts <- function(dt_tracts) {
  ## check contiguity
  is_discontig <- dt_tracts[,
    {
      data.table::setorder(.SD, spos)
      !all(.SD$spos[-1] == .SD$epos[-nrow(.SD)] + 1)
    },
    by = .(sample, chrom)
  ]$V1
  if (any(is_discontig)) {
    stop("Ancestry tracts are not contiguous")
  }

  ## check min/max positions match
  min_unique <- dt_tracts[, .SD$spos[1],
    by = .(sample, chrom)
  ][, length(unique(V1)), by = .(chrom)]$V1 == 1
  if (!all(min_unique)) {
    stop("Not all samples have same start positions")
  }

  max_unique <- dt_tracts[, .SD$epos[.N],
    by = .(sample, chrom)
  ][, length(unique(V1)), by = .(chrom)]$V1 == 1
  if (!all(max_unique)) {
    stop("Not all samples have same end positions")
  }

  ## check no samples missing
  n_counts_unique <- dt_tracts[, length(unique(.SD$sample)), by = .(chrom)]$V1 |>
    unique() |>
    length()
  if (!all(n_counts_unique == 1)) {
    stop("Not all samples exist for each chromosome")
  }
}

#' Convert an "RFMix" or "FLARE" file ancestry file to ".lanc" format
#'
#'
#' `convert_to_lanc` processes a local ancestry input file
#' from RFMIx or FLARE and splits it into ancestry tracts per-sample,
#' generating a file with the ".lanc" format defined by `Admix-kit`.
#'
#' @details
#' This function currently supports only RFMix or FLARE input.
#' It is assumed that the input ancestry file is sorted by position
#' and corresponds to a single chromosome. Samples and positions
#' are matched according to the plink2 pvar and psam files
#' defined by `plink_prefix`. For RFMix, it is assumed that CRF points
#' are contiguous. For FLARE input, which may not be not contiguous,
#' ancestry is imputed to the midpoint between variant positions.
#' After the ancestry input files are processed into tracts,
#' they are written in ".lanc" format to `output`,
#'
#' @param file A string with the file path for the ancestry file to be converted
#' @param file_fmt A string with the file format to be converted
#' (either "FLARE" or "RFMix")
#' @param plink_prefix A string with the plink2 prefix file
#' @param output A string with the file path to save the new ".lanc" file
#' @export
convert_to_lanc <- function(
    file, file_fmt,
    plink_prefix, output) {
  ## Read input to data frame
  if (file_fmt == "FLARE") {
    dt_tracts <- rcpp_read_flare(file)
  } else if (file_fmt == "RFMix") {
    dt_tracts <- rcpp_read_rfmix(file)
  } else {
    stop("Please specify either `FLARE` or `RFMix` input")
  }

  ## Read plink files
  pvar <- data.table::fread(paste0(plink_prefix, ".pvar"), skip = "#CHROM")
  psam <- data.table::fread(paste0(plink_prefix, ".psam"), skip = "#IID")

  ## Convert to data.table and sort/subset samples to pgen
  dt_tracts <- dt_tracts |> data.table::data.table()
  ## Check samples
  if (!all(psam$`#IID` %in% dt_tracts$sample)) {
    stop("Not all pgen samples exist in local ancestry input")
  }
  dt_tracts <- dt_tracts[sample %in% psam$`#IID`, ]
  dt_tracts$sample <- factor(dt_tracts$sample,
    levels = psam$`#IID`
  )
  data.table::setorder(dt_tracts, sample, chrom, spos)

  ## Check that tracts are valid
  validate_tracts(dt_tracts)

  ## Exclude tracts starting after or ending before pgen
  min_pvar <- min(pvar$POS)
  max_pvar <- max(pvar$POS)
  dt_tracts <- dt_tracts[(spos < max_pvar) & (epos > min_pvar), ]

  ## Clip tracts positions to pgen start, end
  dt_tracts$epos[dt_tracts$epos > max_pvar] <- max_pvar
  dt_tracts$spos[dt_tracts$spos < min_pvar] <- min_pvar

  ## Get index of first pvar pos >= tract epos
  dt_tracts$idx <- findInterval(dt_tracts$epos, pvar$POS, left.open = TRUE) + 1

  ## If multiple tracts have same idx, pick last one
  dt_tracts <- dt_tracts[, .SD[.N], by = .(sample, chrom, idx)]

  ## Get .lanc file lines
  dt_tracts[, switch := paste0(idx, ":", anc0, anc1)]
  lines <- dt_tracts[, paste(switch, collapse = " "), by = .(sample, chrom)]$V1

  ## Write output
  header <- paste(nrow(pvar), nrow(psam), sep = " ")
  write(header, output, append = FALSE)
  con <- file(output, "a")
  writeLines(lines, con = con)
  close(con)
}

#' Read .lanc file into list of tracts
#'
#' @param lanc_file A string with the file path to a .lanc file from `Admix-kit`
#' @return A list of data.frames, where each element corresponds to a sample,
#' and has columns `index` (representing indices
#' in a plink2 pvar file) and "anc0" and "anc1", corresponding to ancestries for
#' each haplotype.
#' @importFrom Rcpp sourceCpp
#' @export
read_lanc <- function(lanc_file) {
  rcpp_parse_lanc(readLines(lanc_file)[-1])
}

#' Query ancestry tracts for haplotype ancestry matrices
#'
#' @param query An integer vector of indices corresponding to a
#' .lanc and .pvar file
#' @param tracts A list of ancestry tracts (as returned by `read_lanc`)
#' @return A list with two elements "hap0" and "hap1", containing the
#' haplotype-level ancestry matrices for the query.
#' @noRd
query_tracts <- function(query, tracts) {
  rcpp_query_tracts(query, tracts)
}
