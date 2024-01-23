##' Fit HAUDI model
##'
##' @title HAUDI
##' @param anc_FBM_obj List returned by the function `make_rf_FBM`
##' @param y Vector of responses
##' @param gamma Integer specifying the multiplicative penalty to apply to differences in ancestry-specific effects
##' @param K The number of sets to use for CMSA procedure (analagous to cross-validation)
##' @param ind_train Vector of indices specifying the rows to use for training the model
##' @param family Either "gaussian" or "binomial"
##' @param snps Vector of SNPs to include in model
##' @return An object of class `big_sp_list` from the `bigstatsr` package
##' @author Frank Ockerman
##' @export
HAUDI <- function(anc_FBM_obj, y, gamma, K=10, ind_train=NULL, family, snps=NULL){
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
            K=K,
            ind.col=ind_col
        )
    } else if (family == "binomial"){
        mod <- bigstatsr::big_spLogReg(
            X=anc_FBM_obj$geno,
            y.train=y[ind_train],
            ind.train=ind_train,
            pf.X=pf_X,
            K=K,
            ind.col=ind_col
        )        
    }
    return(mod)
}
