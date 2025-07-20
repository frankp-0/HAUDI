##' Fit HAUDI model
##'
##' @title HAUDI
##' @param fbm_obj object of FBM class
##' @param fbm_info data frame containing information
##' for FBM (chrom/pos/samples/etc.)
##' @param y Vector of responses
##' @param gamma Integer specifying the multiplicative penalty
##' to apply to differences in ancestry-specific effects
##' @param ind_train Vector of indices specifying the
##' rows to use for training the model
##' @param family Either "gaussian" or "binomial"
##' @param snps Vector of SNPs to include in model
##' @param ... additional arguments to pass to big_spLinReg / big_spLogReg
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Frank Ockerman
##' @import bigstatsr
##' @export
haudi <- function(fbm_obj, fbm_info, y, gamma,
                  ind_train = NULL, family, snps = NULL, ...) {
  ## specify training data
  if (is.null(ind_train)) {
    ind_train <- seq_len(nrow(fbm_obj))
  }

  ## specify SNPs to retain
  if (is.null(snps)) {
    col_keep <- rep(TRUE, ncol(fbm_obj))
  } else {
    col_keep <- fbm_info$id %in% snps
  }
  ## remove reference-ancestry columns
  col_keep <- col_keep & (!fbm_info$anc_ref)
  ind_col <- which(col_keep) # subset SNPs

  ## specify multiplicative penalty
  pf_x <- rep(1, ncol(fbm_obj))
  pf_x[fbm_info$anc == "all"] <- gamma
  pf_x <- pf_x[ind_col]

  ## fit HAUDI model
  if (family == "gaussian") {
    mod <- bigstatsr::big_spLinReg(
      X = fbm_obj,
      y.train = y[ind_train],
      ind.train = ind_train,
      pf.X = pf_x,
      ind.col = ind_col,
      ...
    )
  } else if (family == "binomial") {
    mod <- bigstatsr::big_spLogReg(
      X = fbm_obj,
      y01.train = y[ind_train],
      ind.train = ind_train,
      pf.X = pf_x,
      ind.col = ind_col,
      ...
    )
  }
  return(mod)
}

##' @title LASSO
##' @param fbm_obj object of FBM class
##' @param fbm_info data frame containing information
##' for FBM (chrom/pos/samples/etc.)
##' @param y Vector of responses
##' @param ind_train Vector of indices specifying the
##' rows to use for training the model
##' @param family Either "gaussian" or "binomial"
##' @param snps Vector of SNPs to include in model
##' @param ... additional arguments to pass to big_spLinReg / big_spLogReg
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Frank Ockerman
##' @importFrom bigstatsr big_spLinReg big_spLogReg
##' @importClassesFrom bigstatsr FBM.code256
##' @export
lasso <- function(fbm_obj, fbm_info, y, ind_train = NULL,
                  family, snps = NULL, ...) {
  ## specify SNPs to retain
  if (is.null(snps)) {
    col_keep <- rep(TRUE, ncol(fbm_obj))
  } else {
    col_keep <- fbm_info$id %in% snps
  }
  ## remove ancestry-specific columns
  col_keep <- (col_keep) & (fbm_info$anc == "all")
  ind_col <- which(col_keep) # subset SNPs

  ## fit Lasso model
  if (family == "gaussian") {
    mod <- bigstatsr::big_spLinReg(
      X = fbm_obj,
      y.train = y[ind_train],
      ind.train = ind_train,
      ind.col = ind_col,
      ...
    )
  } else if (family == "binomial") {
    mod <- bigstatsr::big_spLogReg(
      X = fbm_obj,
      y01.train = y[ind_train],
      ind.train = ind_train,
      ind.col = ind_col,
      ...
    )
  }
  return(mod)
}



##' Extract population-specific SNP coefficients from a HAUDI model
##'
##' @title get_beta_haudi
##' @param fbm_info data frame containing information
##' for FBM (chrom/pos/samples/etc.)
##' @param haudi_model an object of class big_sp_list,
##' returned by `haudi` or `lasso`
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Frank Ockerman
##' @import bigstatsr
##' @import data.table
##' @export
get_beta_haudi <- function(fbm_info, haudi_model) {
  anc <- snp <- beta_all <- `:=` <- NULL # due to R CMD check

  dt_snp <- data.table::data.table(
    snp = fbm_info$id[attr(haudi_model, "ind.col")],
    beta = summary(haudi_model)$beta[[1]],
    anc = fbm_info$anc[attr(haudi_model, "ind.col")]
  )

  ancestries <- unique(fbm_info$anc)
  ancestries <- ancestries[ancestries != "all"]

  dt_ref <- data.table::data.table(dt_snp[dt_snp$anc == "all", ])
  dt_ref$anc <- NULL
  dt_ref$beta_all <- dt_ref$beta

  for (ancestry in ancestries) {
    x <- dt_snp[dt_snp$anc == ancestry, ]
    idx <- match(dt_ref$snp, x$snp)
    beta_diff <- rep(0, length(idx))
    beta_diff[!is.na(idx)] <- x[idx[!is.na(idx)], ]$beta
    new_col <- paste0("beta_", ancestry)
    dt_ref <- data.table::set(dt_ref,
      j = new_col,
      value = dt_ref$beta_all + beta_diff
    )
  }

  dt_ref$beta_all <- dt_ref$beta <- NULL
  return(dt_ref)
}
