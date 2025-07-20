#' Read a ".lanc" format ancestry file
#'
#' This utility function read ".lanc" files into a data.table
#' where each row is a single per-sample/haplotype ancestry tract.
#'
#' @param lanc_file A string with the ".lanc" file path
#' @param plink_prefix A string with the prefix of the plink2 files
#' @return A data.table
read_lanc <- function(lanc_file, plink_prefix) {
  pvar <- data.table::fread(paste0(plink_prefix, ".pvar"), skip = "#CHROM")
  psam <- data.table::fread(paste0(plink_prefix, ".psam"), skip = "#IID")
  lines <- readLines(lanc_file)[-1]

  dt_tracts <- data.table::data.table()
  for (i in seq_len(length(lines))) {
    ## read line to positions, ancestry
    switches <- strsplit(lines[i], " ")[[1]]
    idx <- sapply(strsplit(switches, ":"), function(x) as.integer(x[[1]]))
    epos <- pvar$POS[idx]
    spos <- c(pvar$POS[1], epos[-length(epos)] + 1)
    anc <- sapply(strsplit(switches, ":"), function(x) x[[2]])
    anc0 <- substr(anc, start = 1, stop = 1) |> as.integer()
    anc1 <- substr(anc, start = 2, stop = 2) |> as.integer()

    ## identify where haplotyeps switch ancestry
    keep0 <- anc0[-length(anc0)] != anc0[-1]
    keep1 <- anc1[-length(anc1)] != anc1[-1]

    ## record ancestry tracts
    dt0 <- data.table::data.table(
      sample = psam$`#IID`[i],
      hap = 0,
      chrom = pvar$`#CHROM`[1],
      spos = spos[c(TRUE, keep0)],
      epos = epos[c(keep0, TRUE)],
      ancestry = anc0[c(TRUE, keep0)]
    )
    dt1 <- data.table::data.table(
      sample = psam$`#IID`[i],
      hap = 1,
      chrom = pvar$`#CHROM`[1],
      spos = spos[c(TRUE, keep1)],
      epos = epos[c(keep1, TRUE)],
      ancestry = anc1[c(TRUE, keep1)]
    )
    dt_tracts <- rbind(dt_tracts, dt0, dt1)
  }
  dt_tracts$chrom <- as.character(dt_tracts$chrom)
  return(dt_tracts)
}

#' Check ancestry tracts are valid
#'
#' @param dt_tracts A data.table with ancestry tracts
validate_tracts <- function(dt_tracts) {
  ## check contiguity
  is_discontig <- dt_tracts[,
    {
      data.table::setorder(.SD, spos)
      !all(.SD$spos[-1] == .SD$epos[-nrow(.SD)] + 1)
    },
    by = .(sample, hap, chrom)
  ]$V1
  if (any(is_discontig)) {
    stop("Ancestry tracts are not contiguous")
  }

  ## check min/max positions match
  min_unique <- dt_tracts[, .SD$spos[1],
    by = .(sample, hap, chrom)
  ][, length(unique(V1)), by = .(chrom)]$V1 == 1
  if (!all(min_unique)) {
    stop("Not all samples/haplotypes have same start positions")
  }

  max_unique <- dt_tracts[, .SD$epos[.N],
    by = .(sample, hap, chrom)
  ][, length(unique(V1)), by = .(chrom)]$V1 == 1
  if (!all(max_unique)) {
    stop("Not all samples/haplotypes have same start positions")
  }

  ## check no samples missing
  n_counts_unique <- dt_tracts[, length(unique(.SD$sample)), by = .(hap, chrom)]$V1 |>
    unique() |>
    length()
  if (!all(n_counts_unique == 1)) {
    stop("Not all samples exist for each chromosome/haplotype")
  }
}


#' Read ancestry tracts to a data.table
#'
#' @description
#' `read_ancestry_tracts` processes a local ancestry input file
#' and splits it into ancestry tracts per-sample and per-haplotype,
#' returning a data.table with one row for each tract.
#'
#' @details
#' This function currently supports RFMix, FLARE, or the .lanc input
#' format defined by `Admix-kit`. Multi-chromosome input files are
#' allowed for RFMix and FLARE but not .lanc. It is assumed that
#' input files are sorted by chromosome and position. For RFMix,
#' it is assumed that CRF points are contiguous. When the ancestry
#' of a haplotype switches (or when a new chromosome begins),
#' the haplotype is split into ancestry tracts.
#' For FLARE input, which is not contiguous,
#' ancestry is imputed to the midpoint between variants.
#' For .lanc input, tracts are contiguous by definition.
#'
#' @param file Path for the local ancestry input file
#' @param file_fmt Either "FLARE", "RFMix", or "lanc", specifying
#' the format of the local ancestry file.
#' @param extend_tracts Boolean indicating whether to
#' impute ancestry at before/after the minimum/maximum
#' position on a chromosome. This ensures that any
#' position on the chromosome may be queried for ancestry.
#' that all
#'
#' @return A data.table with one row per-ancestry tract, containing
#' columns "sample", "hap" (haplotype), "chrom" (chromosome), "spos"
#'  (start position), "epos" (end position) and "ancestry" (an integer).
#' @author Franklin Ockerman
#'
#' @importFrom data.table data.table
#' @useDynLib HAUDI
#' @export
read_ancestry_tracts <- function(
    file, file_fmt, extend_tracts = FALSE,
    plink_prefix = NULL) {
  ## Read inpput to data frame
  if (file_fmt == "FLARE") {
    dt_tracts <- rcpp_read_flare(file)
  } else if (file_fmt == "RFMix") {
    dt_tracts <- rcpp_read_rfmix(file)
  } else if (file_fmt == "lanc") {
    if (is.null(plink_prefix)) {
      stop("Plink prefix must be provided for `.lanc` format")
    }
    dt_tracts <- read_lanc(file, plink_prefix)
  } else {
    stop("Please specify either `FLARE`, `RFMix`, or `lanc` input")
  }

  ## Convert to data.table and sort
  dt_tracts <- dt_tracts |> data.table::data.table()
  data.table::setorder(dt_tracts, sample, hap, chrom, spos)

  ## Extend tracts if requested
  if (extend_tracts) {
    dt_tracts[, spos := ifelse(spos == min(spos), 0, spos), by = .(sample, hap, chrom)]
    dt_tracts[, epos := ifelse(epos == max(epos), 250e6, epos), by = .(sample, hap, chrom)]
  }

  ## Check that tracts are valid
  validate_tracts(dt_tracts)

  return(dt_tracts)
}


#' Convert an "RFMix" or "FLARE" file to ".lanc" format
#'
#' @param file A string with the file path for the ancestry file to be converted
#' @param file_fmt A string with the file format to be converted
#' (either "FLARE" or "RFMix")
convert_to_lanc <- function(
    file, file_fmt,
    plink_prefix, output) {
  ## Read input to data frame
  print(sprintf("Reading ancestry tracts for %s", file))
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

  ## Convert to data.table and sort
  dt_tracts <- dt_tracts |> data.table::data.table()
  data.table::setorder(dt_tracts, sample, hap, chrom, spos)

  ## Check that tracts are valid
  validate_tracts(dt_tracts)

  ## Extend end of tract and match with pgen
  dt_tracts[, epos := ifelse(epos == max(epos), max(pvar$POS), epos), by = .(sample, hap, chrom)]

  ## Check samples
  if (!all(psam$`#IID` %in% dt_tracts$sample)) {
    stop("Not all pgen samples exist in local ancestry input")
  }
  dt_tracts <- dt_tracts[sample %in% psam$`#IID`, ]

  print("Matching tracts to pvar")
  dt_tracts$idx <- sapply(seq_len(nrow(dt_tracts)), function(i) {
    which(pvar$POS >= dt_tracts$epos[i])[1]
  })


  ## Write output
  print("Writing output")
  header <- paste(nrow(pvar), nrow(psam), sep = " ")
  write(header, output, append = FALSE)
  lines <- sapply(psam$`#IID`, function(samp) {
    dt_samp <- dt_tracts[sample == samp, ]
    samp_idx <- unique(dt_samp)$idx |> sort()
    anc_mat <- matrix(NA, 2, length(samp_idx))
    for (i in seq_along(samp_idx)) {
      anc_mat[1, i:length(samp_idx)] <- dt_samp[hap == 0 & idx >= samp_idx[i], .SD[1]]$ancestry
      anc_mat[2, i:length(samp_idx)] <- dt_samp[hap == 1 & idx >= samp_idx[i], .SD[1]]$ancestry
    }
    anc <- apply(anc_mat, 2, function(x) paste(x, collapse = ""))
    switches <- paste(samp_idx, anc, sep = ":")
    paste(switches, collapse = " ")
  })
  con <- file(output, "a")
  writeLines(lines, con = con)
  close(con)
}
