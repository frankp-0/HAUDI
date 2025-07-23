# TODO: Documentation

##' @export
haudi <- function(fbm, fbm_info, y_train, gamma_vec,
                  family, k = 10, variants = NULL, ...) {
  ## Check for disallowed arguments
  dots <- list(...)
  bad_args <- intersect(
    names(dots),
    c("X", "pf.X", "K", "y.train", "ind.col")
  )
  if (length(bad_args) > 0) {
    stop("Please do not specify: ", paste(bad_args, collapse = ", "))
  }



  ## Specify variants to retain
  if (is.null(variants)) {
    col_keep <- rep(TRUE, ncol(fbm))
  } else {
    col_keep <- fbm_info$id %in% variants
  }
  col_keep[fbm_info$anc_ref == TRUE] <- FALSE
  ind_col <- which(col_keep)

  cv_results <- data.frame(gamma = gamma_vec, validation_loss = NA_real_)

  for (gamma in gamma_vec) {
    ## specify multiplicative penalty
    pf_x <- rep(1, ncol(fbm))
    pf_x[fbm_info$anc == "all"] <- gamma
    pf_x <- pf_x[ind_col]

    ## fit HAUDI model
    if (family == "gaussian") {
      model <- bigstatsr::big_spLinReg(
        X = fbm,
        y.train = y_train,
        pf.X = pf_x,
        ind.col = ind_col,
        K = k,
        ...
      )
    } else if (family == "binomial") {
      model <- bigstatsr::big_spLogReg(
        X = fbm,
        y01.train = y_train,
        pf.X = pf_x,
        ind.col = ind_col,
        K = k,
        ...
      )
    }
    i <- which(cv_results$gamma == gamma)
    loss_val <- summary(model)$validation_loss
    cv_results$validation_loss[i] <- loss_val
    if (loss_val == max(cv_results$validation_loss, na.rm = TRUE)) {
      final_model <- model
    }
  }
  return(list(model = final_model, cv_results = cv_results))
}


##' @export
lasso <- function(fbm, fbm_info, y_train,
                  family, k = 10, variants = NULL, ...) {
  ## Check for disallowed arguments
  dots <- list(...)
  bad_args <- intersect(
    names(dots),
    c("X", "K", "y.train", "ind.col")
  )
  if (length(bad_args) > 0) {
    stop("Please do not specify: ", paste(bad_args, collapse = ", "))
  }

  ## specify SNPs to retain
  if (is.null(variants)) {
    col_keep <- rep(TRUE, ncol(fbm))
  } else {
    col_keep <- fbm_info$id %in% variants
  }

  ## remove ancestry-specific columns
  col_keep <- (col_keep) & (fbm_info$anc == "all")
  ind_col <- which(col_keep) # subset SNPs

  ## fit Lasso modelel
  if (family == "gaussian") {
    model <- bigstatsr::big_spLinReg(
      X = fbm,
      y.train = y_train,
      ind.col = ind_col,
      ...
    )
  } else if (family == "binomial") {
    model <- bigstatsr::big_spLogReg(
      X = fbm,
      y01.train = y_train,
      ind.col = ind_col,
      ...
    )
  }
  return(model)
}


##' @export
get_beta_haudi <- function(fbm_info, haudi_model) {
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
