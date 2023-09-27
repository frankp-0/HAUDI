##' Train HAUDI model
##'
##' @param anc_FBM_obj a list containing ancestry-genotype information 
##' @param y a vector of phenotypes
##' @param gamma an integer specifying the multiplicative penalty applied to ancestry-specific differences in effect
##' @param ind_train a vector of indices for training data
##' @param family either "gaussian" or "binomial"
##' @return returns the model, an object of class `bigstatsr::big_sp_list`
##' @author Frank Ockerman
##' @importFrom bigstatsr big_spLinReg big_spLogReg

HAUDI <- function(anc_FBM_obj, y, gamma, ind_train, family){
    ind.col <- which(!anc_FBM_obj$anc_ref)    
    pf.X <- rep(gamma, ncol(anc_FBM_obj$geno))
    pf.X[anc_FBM_obj$anc == "all"] <- 1
    pf.X <- pf.X[ind.col]
    if(family == "gaussian"){
        mod <- big_spLinReg(
            X=anc_FBM_obj$geno,
            y.train=y[ind_train],
            ind.train=ind_train,
            pf.X=pf.X,
            ind.col=ind.col
        )
    } else if (family == "binomial"){
        mod <- big_spLogReg(
            X=anc_FBM_obj$geno,
            y.train=y[ind_train],
            ind.train=ind_train,
            pf.X=pf.X,
            ind.col=ind.col
        )        
    }
    return(mod)
}
