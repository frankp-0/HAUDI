##' Fit GAUDI model
##'
##' @title GAUDI
##' @param FBM_obj object of FBM class
##' @param FBM_info data frame containing information for FBM (chrom/pos/samples/etc.)
##' @param y vector of responses
##' @param K number of folds to use in cross-validation
##' @param ind_train integer vector of samples to use in training model
##' @param snps vector of snp names to use in model
##' @param verbose boolean indicating verbosity
##' @param minlam minimum lambda value for genlasso
##' @param maxsteps_init max steps for genlasso in initial model
##' @param maxsteps_cv max steps for genlasso in cross-validation models
##' @param gamma_vec vector of gamma values
##' @return list containing information for best cross-validated fit 
##' @author Bryce Rowland, Frank Ockerman
##' @import genlasso
##' @importFrom splitTools create_folds
##' @importFrom stringr str_remove
##' @importFrom Matrix sparseMatrix
##' @importClassesFrom Matrix dgCMatrix
##' @export
GAUDI <- function(FBM_obj, FBM_info, y, gamma_vec, K=10, ind_train=NULL, snps=NULL,
                  verbose = F, minlam = 0, maxsteps_init = 2000,
                  maxsteps_cv = 2000 * 10){
    X <- construct_GAUDI(FBM_obj, FBM_info, snps)
    if(is.null(ind_train)) {ind_train <- 1:nrow(X)}
    X <- X[ind_train,]
    y <- y[ind_train]
    mod <- cv_fused_lasso(X=X, y=y, n_folds=K, verbose=verbose,
                          minlam=minlam, maxsteps_init=maxsteps_init,
                          maxsteps_cv=maxsteps_cv, gamma_vec=gamma_vec)
    return(mod)
}

cv_fused_lasso <- function(X, y, n_folds, 
                           verbose = F, 
                           minlam = 0,
                           maxsteps_init = 2000,
                           maxsteps_cv = 2000 * 10,
                           gamma_vec
                           ){
    
    gamma_scores <- list(gamma = gamma_vec, 
                         init_fit = vector("list", length(gamma_vec)), 
                         D = vector("list", length(gamma_vec)),
                         best_lambda = rep(0, length(gamma_vec)), 
                         best_r2 = rep(0, length(gamma_vec)), 
                         cv_r2_full = vector("list", length(gamma_vec)))
    
    i <- 1
    ## for each value of gamma, loop through this procedure. 
    for(gamma in gamma_vec){
        
        
        penalty_matrix <- get_penalty_matrix(X, gamma = gamma)
        
        if(is.null(penalty_matrix)){
            return(NA)
        }
        
        print(sprintf("Getting lambdas from initial genlasso call for gamma value %s", gamma))
        init_fit <- genlasso::genlasso(y = y, 
                             X = X, 
                             verbose = verbose, 
                             D = penalty_matrix,
                             maxsteps = maxsteps_init)
        print("Done!")
        
        splits <- create_folds(y = 1:nrow(X), 
                               k = n_folds, 
                               type = "basic")
        ## Set minlam for each fold to the min lam from the initial fit.
        ## Then, below we specify to run the cv fold models to the minimum lambda, and specify a large number
        ## of itterations to make sure it gets there. But stop once it does!
        minlam_fold <- min(init_fit$lambda)
        
        lambdas <- init_fit$lambda
        fold_scores <- vector("list", length = n_folds)
        
        for(k in 1:n_folds){
            
            print(sprintf("Fitting fold %s", k))
            
            fold_fit <- genlasso(y = y[splits[[k]]], 
                                 X = X[splits[[k]],], 
                                 D = penalty_matrix, 
                                 minlam = minlam_fold, 
                                 maxsteps = maxsteps_cv, 
                                 verbose = verbose)
            fold_pred <- predict.genlasso(fold_fit,
                                          lambda = lambdas,
                                          Xnew = X[-splits[[k]],])$fit
            fold_scores[[k]] <- cor(fold_pred, y[-splits[[k]]]) ^ 2
        }
        
        cv_r2_full <- do.call(cbind, fold_scores)
        cv_r2 <- rowMeans(do.call(cbind, fold_scores))
        best_lambda <- lambdas[which(cv_r2 == max(cv_r2, na.rm = T))][1]
        best_r2 <- cv_r2[which(cv_r2 == max(cv_r2, na.rm = T))][1]
        
        gamma_scores$init_fit[[i]] <- init_fit
        gamma_scores$D[[i]] <- penalty_matrix
        gamma_scores$cv_r2_full[[i]] <- cv_r2_full
        gamma_scores$best_lambda[i] <- best_lambda
        gamma_scores$best_r2[i] <- best_r2
        
        i <- i + 1  

    }
    
    best_index <- which(gamma_scores$best_r2 == max(gamma_scores$best_r2))
    
    out_list <- list("fit" = gamma_scores$init_fit[[best_index]], 
                     "D" = gamma_scores$D[[best_index]], 
                     "snps" = colnames(X),
                     "best_gamma" = gamma_scores$gamma[best_index],
                     "best_lambda" = gamma_scores$best_lambda[best_index], 
                     "cv_r2" = gamma_scores$best_r2[best_index], 
                     "full_cv_r2" = gamma_scores$cv_r2_full) 
    
    return(out_list)
}

get_penalty_matrix <- function(X, gamma){
    if(gamma > 0){
        rbind(get_upper_fusion_matrix(X), diag(ncol(X)) * gamma)
    } else{
        get_upper_fusion_matrix(X)
    }
}

get_cv_fl_best_betas <- function(.fl_fit){
    as.numeric(coef(.fl_fit$fit, lambda = .fl_fit$best_lambda)$beta)
}

get_pair_fusion_matrix <- function(p){
    genlasso::getD1dSparse(p)[c(1:p)[which(c(1:p) %% 2 != 0)],]
}

get_upper_fusion_matrix <- function(X){
    no_suffix <- stringr::str_remove(colnames(X), ".anc.*")
    
    p_pairs <- length(unique(no_suffix[duplicated(no_suffix)]))
    
    if(p_pairs <= 1){
        print("WARNING: Paired variants <= 1. Skipping GAUDI fit.")
        return(NULL)
    }
    
    fusion_D <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(p_pairs, ncol(X)))
    
    ## Get a logical vector for duplicate positions
    pair_col_indices <- (no_suffix %in% unique(no_suffix[duplicated(no_suffix)]))
    
    ## Get paired fused matrix for p pairs. 
    fm <- get_pair_fusion_matrix(sum(pair_col_indices))
    
    ## Set to the correct columns for our penalty. 0 for single snps indicates
    ##  no fusion. 
    fusion_D[,pair_col_indices] <- fm
    fusion_D
}

##' Constructs X matrix for GAUDI
##'
##' @title construct_GAUDI
##' @param anc_FBM_obj 
##' @param snps 
##' @return matrix for fitting GAUDI
##' @author Frank Ockerman
##' @export
construct_GAUDI <- function(FBM_obj, FBM_info, snps=NULL){
    if(is.null(snps)) {
        idx_col <- 1:ncol(FBM_obj)
        idx_col <- idx_col[FBM_info$anc[idx_col] != "all"]        
    } else {
        idx_col <- which(FBM_info$rsid %in% snps)
        idx_col <- idx_col[FBM_info$anc[idx_col] != "all"]        
    }
    rsid <- FBM_info$rsid[idx_col]
    idx_col <- idx_col[order(rsid)]
    anc <- FBM_info$anc[idx_col]
    rsid <- rsid[order(rsid)]
    X <- FBM_obj[,idx_col]
    colnames(X) <- paste(rsid, anc, sep=".anc.")
    return(X)
}
