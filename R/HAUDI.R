##' Fit HAUDI model
##'
##' @title HAUDI
##' @param FBM_obj object of FBM class
##' @param FBM_info data frame containing information for FBM (chrom/pos/samples/etc.)
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
HAUDI <- function(FBM_obj, FBM_info, y, gamma, ind_train=NULL, family, snps=NULL, ...){
    ## specify training data
    if(is.null(ind_train)){ind_train = 1:nrow(FBM)}
    
    ## specify SNPs to retain
    if(is.null(snps)) {
        col_keep <- rep(TRUE, ncol(FBM_obj))
    } else {
        col_keep <- FBM_info$rsid %in% snps
    }
    ## remove reference-ancestry columns
    col_keep <- col_keep & (!FBM_info$anc_ref)
    ind_col <- which(col_keep) # subset SNPs

    ## specify multiplicative penalty
    pf_X <- rep(gamma, ncol(FBM))
    pf_X[FBM_info$anc == "all"] <- 1
    pf_X <- pf_X[ind_col]

    ## fit HAUDI model
    if(family == "gaussian"){
        mod <- bigstatsr::big_spLinReg(
                              X=FBM_obj,
                              y.train=y[ind_train],
                              ind.train=ind_train,
                              pf.X=pf_X,
                              ind.col=ind_col,
                              ...)
    } else if (family == "binomial"){
        mod <- bigstatsr::big_spLogReg(
                              X=FBM_obm,
                              y01.train=y[ind_train],
                              ind.train=ind_train,
                              pf.X=pf_X,
                              ind.col=ind_col,
                              ...)
    }
    return(mod)
}

##' @title LASSO
##' @param FBM_obj object of FBM class
##' @param FBM_info data frame containing information for FBM (chrom/pos/samples/etc.)
##' @param y Vector of responses
##' @param ind_train Vector of indices specifying the rows to use for training the model
##' @param family Either "gaussian" or "binomial"
##' @param snps Vector of SNPs to include in model
##' @param ... additional arguments to pass to big_spLinReg / big_spLogReg
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Frank Ockerman
##' @importFrom bigstatsr big_spLinReg big_spLogReg
##' @importClassesFrom bigstatsr FBM.code256
##' @export
LASSO <- function(FBM_obj, FBM_info, y, gamma, ind_train=NULL, family, snps=NULL, ...){
    ## specify SNPs to retain
    if(is.null(snps)) {
        col_keep <- rep(TRUE, ncol(FBM_obj))
    } else {
        col_keep <- FBM_info$rsid %in% snps
    }
    ## remove ancestry-specific columns
    col_keep <- (col_keep) & (FBM_info$anc == "all")
    ind_col <- which(col_keep) # subset SNPs

    ## fit Lasso model
    if(family == "gaussian"){
        mod <- bigstatsr::big_spLinReg(
                              X=FBM_obj,
                              y.train=y[ind_train],
                              ind.train=ind_train,
                              ind.col=ind_col,
                              ...)
    } else if (family == "binomial"){
        mod <- bigstatsr::big_spLogReg(
                              X=FBM_obj,
                              y01.train=y[ind_train],
                              ind.train=ind_train,
                              ind.col=ind_col,
                              ...)
    }
    return(mod)
}



##' @export
get_beta_HAUDI <- function(FBM_obj, FBM_info, HAUDI_model){
    dt_snp <- data.table(
        snp = FBM_info$rsid[attr(HAUDI_model, "ind.col")],
        beta = summary(HAUDI_model)$beta[[1]],
        anc = FBM_info$anc[attr(HAUDI_model, "ind.col")]
    )

    ancestries <- unique(FBM_info$anc)
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
