#' Fit GAUDI model
#'
#' Calculates GAUDI admixture-aware polygenic scores
#' using a bigstatsr FBM.code256 object generated
#' by `make_fbm`.
#'
#' This function is a re-implementation of GAUDI
#' provided so that it may be easily incorporated
#' into HAUDI workflows.
#'
#' @inheritParams haudi
#' @param ind_train An optional vector of indices corresponding to rows
#' in the fbm matrix, used to subset samples for training. If provided,
#' y_train must be the same length and sorted to these indices.
#' @param verbose A boolean indicating verbosity
#' @param minlam A numeric with the minimum lambda value for genlasso
#' @param maxsteps_init An integer with the max steps for genlasso in
# the initial model
#' @param maxsteps_cv An integer with the max steps for genlasso in
# cross-validation models
#' @param gamma_vec A numeric vector with gamma values
#' @param approx A boolean indicating whether to use the approximate solution
#' @param splits A list of k vectors indicating the indices
#' (in the training set) for inner-CV folds
#' @importFrom stats cor
#' @importFrom genlasso fusedlasso getD1dSparse predict.genlasso
#' @importFrom splitTools create_folds
#' @importFrom stringr str_remove
#' @importFrom Matrix sparseMatrix
#' @importFrom stats na.omit
#' @export
gaudi <- function(fbm, fbm_info, y_train, ind_train = NULL,
                  gamma_vec, k = 10, variants = NULL,
                  verbose = FALSE, minlam = 0,
                  maxsteps_init = 2000, maxsteps_cv = 2000 * 10,
                  approx = FALSE, splits = NULL) {
  x_mat <- construct_gaudi(fbm, fbm_info, variants)
  if (is.null(ind_train)) {
    ind_train <- seq_len(nrow(x_mat))
  }
  x_mat <- x_mat[ind_train, ]

  if (is.null(splits)) {
    ind_sets <- sample(ind_train)
    ind_sets <- as.numeric(cut(seq_along(ind_sets), k))[order(ind_sets)]
    splits <- lapply(1:k, function(cv_fold) which(ind_sets != cv_fold))
  }

  model <- cv_fused_lasso(
    x = x_mat, y = y_train, n_folds = k, verbose = verbose,
    minlam = minlam, maxsteps_init = maxsteps_init,
    maxsteps_cv = maxsteps_cv, gamma_vec = gamma_vec,
    approx = approx, splits = splits
  )
  return(model)
}

cv_fused_lasso <- function(x, y, n_folds,
                           verbose = FALSE,
                           minlam = 0,
                           maxsteps_init = 2000,
                           maxsteps_cv = 2000 * 10,
                           gamma_vec, approx = FALSE,
                           splits) {
  gamma_scores <- list(
    gamma = gamma_vec,
    init_fit = vector("list", length(gamma_vec)),
    D = vector("list", length(gamma_vec)),
    best_lambda = rep(0, length(gamma_vec)),
    best_r2 = rep(0, length(gamma_vec)),
    cv_r2_full = vector("list", length(gamma_vec))
  )

  i <- 1
  ## for each value of gamma, loop through this procedure.
  for (gamma in gamma_vec) {
    penalty_matrix <- get_upper_fusion_matrix(x)

    if (is.null(penalty_matrix)) {
      return(NA)
    }

    print(sprintf(
      "Getting lambdas from initial genlasso call for gamma value %s",
      gamma
    ))
    init_fit <- genlasso::fusedlasso(
      y = y, X = x, D = penalty_matrix, gamma = gamma, maxsteps = maxsteps_init,
      verbose = verbose, approx = approx
    )
    print("Done!")

    ## Set minlam for each fold to the min lam from the initial fit.
    ## Then, below we specify to run the cv fold models to the
    ## minimum lambda, and specify a large number of iterations to
    ## make sure it gets there. But stop once it does!
    lambdas <- init_fit$lambda
    fold_scores <- vector("list", length = n_folds)
    n_ratio <- (n_folds - 1) / n_folds

    for (fold in 1:n_folds) {
      print(sprintf("Fitting fold %s", fold))
      fold_fit <- genlasso::fusedlasso(
        y = y[splits[[fold]]],
        X = x[splits[[fold]], ],
        D = penalty_matrix,
        gamma = gamma,
        minlam = n_ratio * min(init_fit$lambda),
        maxsteps = maxsteps_cv,
        verbose = verbose,
        approx = approx
      )
      fold_pred <- predict.genlasso(fold_fit,
        lambda = lambdas * n_ratio,
        Xnew = x[-splits[[fold]], ]
      )$fit
      fold_scores[[fold]] <- stats::cor(fold_pred, y[-splits[[fold]]])^2
    }

    cv_r2_full <- do.call(cbind, fold_scores)
    cv_r2 <- rowMeans(do.call(cbind, fold_scores))
    best_lambda <- lambdas[which(cv_r2 == max(cv_r2, na.rm = TRUE))][1]
    best_r2 <- cv_r2[which(cv_r2 == max(cv_r2, na.rm = TRUE))][1]

    gamma_scores$init_fit[[i]] <- init_fit
    gamma_scores$D[[i]] <- penalty_matrix
    gamma_scores$cv_r2_full[[i]] <- cv_r2_full
    gamma_scores$best_lambda[i] <- best_lambda
    gamma_scores$best_r2[i] <- best_r2
    i <- i + 1
  }

  best_index <- which(gamma_scores$best_r2 == max(gamma_scores$best_r2))

  out_list <- list(
    "fit" = gamma_scores$init_fit[[best_index]],
    "D" = gamma_scores$D[[best_index]],
    "snps" = colnames(x),
    "best_gamma" = gamma_scores$gamma[best_index],
    "best_lambda" = gamma_scores$best_lambda[best_index],
    "cv_r2" = gamma_scores$best_r2[best_index],
    "full_cv_r2" = gamma_scores$cv_r2_full
  )

  return(out_list)
}

get_pair_fusion_matrix <- function(p) {
  genlasso::getD1dSparse(p)[c(1:p)[which(c(1:p) %% 2 != 0)], ]
}

get_upper_fusion_matrix <- function(x) {
  no_suffix <- stringr::str_remove(colnames(x), ".anc.*")

  p_pairs <- length(unique(no_suffix[duplicated(no_suffix)]))

  if (p_pairs <= 1) {
    print("WARNING: Paired variants <= 1. Skipping gaudi fit.")
    return(NULL)
  }

  fusion_d <- Matrix::sparseMatrix(
    i = 1, j = 1, x = 0,
    dims = c(p_pairs, ncol(x))
  )

  ## Get a logical vector for duplicate positions
  pair_col_indices <- (no_suffix %in% unique(no_suffix[duplicated(no_suffix)]))

  ## Get paired fused matrix for p pairs.
  fm <- get_pair_fusion_matrix(sum(pair_col_indices))

  ## Set to the correct columns for our penalty. 0 for single snps indicates
  ##  no fusion.
  fusion_d[, pair_col_indices] <- fm
  fusion_d
}

#' Constructs X matrix for GAUDI
#'
#' @inheritParams haudi
#' @return The GAUDI model matrix
#' @export
construct_gaudi <- function(fbm, fbm_info, variants = NULL) {
  if (is.null(variants)) {
    ind_col <- seq_len(ncol(fbm))
    ind_col <- ind_col[fbm_info$anc[ind_col] != "all"]
  } else {
    ind_col <- which(fbm_info$id %in% variants)
    ind_col <- ind_col[fbm_info$anc[ind_col] != "all"]
  }
  id <- fbm_info$id[ind_col]
  ind_col <- ind_col[order(id)]
  anc <- fbm_info$anc[ind_col]
  id <- id[order(id)]
  x <- cbind(1, fbm[, ind_col])
  colnames(x) <- c("Intercept", paste(id, anc, sep = ".anc."))
  return(x)
}
