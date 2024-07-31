##' Fit HAUDI model
##'
##' @title HAUDI
##' @param anc_FBM_obj List returned by the function `make_rf_FBM`
##' @param y Vector of responses
##' @param gamma Integer specifying the multiplicative penalty to apply to differences in ancestry-specific effects
##' @param ind_train Vector of indices specifying the rows to use for training the model
##' @param family Either "gaussian" or "binomial"
##' @param snps Vector of SNPs to include in model
##' @param ... additional arguments to pass to big_spLinReg / big_spLogReg
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Frank Ockerman
##' @importFrom bigstatsr big_spLinReg big_spLogReg
##' @importClassesFrom bigstatsr FBM.code256
##' @export
HAUDI <- function(anc_FBM_obj, y, gamma, ind_train=NULL, family, snps=NULL, ...){
  ## specify SNPs to retain
  if(is.null(snps)) {
    col_keep <- rep(TRUE, ncol(anc_FBM_obj$geno))
  } else {
    col_keep <- anc_FBM_obj$rsid %in% snps
  }
  ## remove reference-ancestry columns
  col_keep <- col_keep & (!anc_FBM_obj$anc_ref)
  ind_col <- which(col_keep) # subset SNPs

  ## specify multiplicative penalty
  pf_X <- rep(gamma, ncol(anc_FBM_obj$geno))
  pf_X[anc_FBM_obj$anc == "all"] <- 1
  pf_X <- pf_X[ind_col]

  ## fit HAUDI model
  if(family == "gaussian"){
    mod <- bigstatsr::big_spLinReg(
                        X=anc_FBM_obj$geno,
                        y.train=y[ind_train],
                        ind.train=ind_train,
                        pf.X=pf_X,
                        ind.col=ind_col,
                        ...)
  } else if (family == "binomial"){
    mod <- bigstatsr::big_spLogReg(
                        X=anc_FBM_obj$geno,
                        y01.train=y[ind_train],
                        ind.train=ind_train,
                        pf.X=pf_X,
                        ind.col=ind_col,
                        ...)
  }
  return(mod)
}

##' @export
LASSO <- function(anc_FBM_obj, y, gamma, ind_train=NULL, family, snps=NULL, ...){
  ## specify SNPs to retain
  if(is.null(snps)) {
    col_keep <- rep(TRUE, ncol(anc_FBM_obj$geno))
  } else {
    col_keep <- anc_FBM_obj$rsid %in% snps
  }
  ## remove ancestry-specific columns
  col_keep <- (col_keep) & (anc_FBM_obj$anc == "all")
  ind_col <- which(col_keep) # subset SNPs

  ## fit Lasso model
  if(family == "gaussian"){
    mod <- bigstatsr::big_spLinReg(
                        X=anc_FBM_obj$geno,
                        y.train=y[ind_train],
                        ind.train=ind_train,
                        ind.col=ind_col,
                        ...)
  } else if (family == "binomial"){
    mod <- bigstatsr::big_spLogReg(
                        X=anc_FBM_obj$geno,
                        y01.train=y[ind_train],
                        ind.train=ind_train,
                        ind.col=ind_col,
                        ...)
  }
  return(mod)
}



##' @export
get_beta_HAUDI <- function(anc_FBM_obj, HAUDI_model){
  dt_snp <- data.table(
    snp = anc_FBM_obj$rsid[attr(HAUDI_model, "ind.col")],
    beta = summary(HAUDI_model)$beta[[1]],
    anc = anc_FBM_obj$anc[attr(HAUDI_model, "ind.col")]
  )

  ancestries <- unique(anc_FBM_obj$anc)
  ancestries <- ancestries[ancestries != "all"]

  dt_ref <- dt_snp[anc == "all",-'anc']
  dt_ref$beta_all <- dt_ref$beta

  for (ancestry in ancestries){
    x <- dt_snp[anc == ancestry,]
    beta_diff <- dt_ref[, x[match(dt_ref$snp, snp), ]$beta]
    beta_diff[is.na(beta_diff)] <- 0
    dt_ref[, (paste0("beta_", ancestry)) := beta_all + beta_diff]
  }

  dt_ref$beta_all <- dt_ref$beta <- NULL
  return(dt_ref)
}

