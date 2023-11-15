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

#library(data.table)
#pheno <- fread("/work/users/f/r/frocko/HAUDI/realDataAnalysis/derivedData/PAGE/HCHS_SOL_pheno.txt")
#anc_FBM_obj <- readRDS("../HAUDI_PAGE/realDataAnalysis/fbm.rds")
#fold1 <- fread("/work/users/f/r/frocko/HAUDI/realDataAnalysis/derivedData/PAGE/HCHS_SOL_samples_fold1.txt")[[1]]
#pheno <- pheno[match(anc_FBM_obj$samples, page_subject_id),]
#pheno_NA <- data.table(is.na(pheno))
#names(pheno_NA) <- names(pheno)
#pheno[is.na(pheno)] <- 0
#y <- pheno$height
#ind_train <- which(!(pheno_NA$height) & (pheno$page_subject_id %in% fold1))
#family <- "gaussian"
#gamma_vec <- c(1,1.5,2,2.5,3)
cv_HAUDI <- function(anc_FBM_obj, y, gamma_vec, ind_train, family){
    ind_train_shuffed <- sample(ind_train, size=length(ind_train), replace=FALSE)
    folds_shuffed <- rep(1:length(gamma_vec), length.out = length(ind_train))
    folds <- folds_shuffed[match(ind_train, ind_train_shuffed)]
    models <- list()
    for (fold in 1:length(gamma_vec)){
        ind_train_cv <- ind_train[folds != fold]
        mod <- HAUDI(anc_FBM_obj, y, gamma_vec[fold], ind_train_cv, family)
        }

    }
