#' Compute Approximate Complete Data Log-Likelihood
#'
#' Computes the approximate complete data log-likelihood for differential expression 
#' analysis under a mixture model framework. Calculates weighted log-likelihoods for 
#' each gene-cell type combination, assuming normal distributions with observed 
#' weights accounting for cell type-viral type assignments.
#'
#' @param Y Numeric matrix of observed gene expression data with dimensions 
#'   (cells x genes), where rows represent cells and columns represent genes.
#' @param M Numeric 3-dimensional array of estimated mean parameters with dimensions 
#'   (genes x cell types x viral types), containing the mean gene expression 
#'   for each combination of gene, cell type, and viral type.
#' @param W Numeric 3-dimensional array of cell-type and viral-type weights with 
#'   dimensions (cells x cell types x viral types), representing the soft assignment 
#'   or responsibility of each cell to each cell type-viral type combination.
#' @param sigma2 Numeric 3-dimensional array of estimated variance parameters with 
#'   dimensions (genes x cell types x viral types), containing the variance of gene 
#'   expression for each combination of gene, cell type, and viral type.
#'
#' @return A numeric matrix with dimensions (genes x cell types) containing the 
#'   approximate log-likelihood contribution for each gene in each cell type, 
#'   aggregated across viral types and weighted by the assignment weights.
#'
#' @keywords internal
#' @export
approx_complete_data_loglik_fast <- function(Y, M, W, sigma2) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  resout <- approx_complete_data_loglik_rcpp(Y, M, W, sigma2)
  
  rownames(resout) <- dimnames(M)[[1]]
  colnames(resout) <- dimnames(M)[[2]]
  return(resout)
}

approx_complete_data_loglik_pair_fast <- function(Y, M.alt, sigma2.alt, M.null, sigma2.null, W) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  approx_complete_data_loglik_pair_rcpp(Y, M.alt, sigma2.alt, M.null, sigma2.null, W)
}

#' Estimate Mean Expression per Gene, Cell Type, and Condition (DE)
#'
#' Computes weighted mean expression matrices for the alternative hypothesis
#' in differential expression testing, where specified conditions show expression
#' differences while others share effects.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)},
#'   where \eqn{C} is the number of cell types (or classes),
#'   and \eqn{V} is the number of conditions (or contrasts).
#'   Each slice \code{W[i,,]} contains the posterior weights for
#'   sample \eqn{i}.
#' @param compStatus Integer or character vector specifying which
#'   conditions should have separate mean estimates. Other conditions
#'   will have their means estimated from the combined data.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)}
#'   containing the weighted mean expression values. \code{M.out[j, c, v]}
#'   is the estimated mean of gene \eqn{j} in cell type \eqn{c} under
#'   condition \eqn{v}.
#'
#' @details
#' For conditions specified in \code{compStatus}, mean estimates are
#' computed separately using only samples from those conditions.
#' For other conditions, means are estimated as the average across all samples,
#' using mean weights to account for condition effects.
#'
#' keywords internal
resolve_comp_status_indices <- function(W, compStatus) {
  if (is.character(compStatus)) {
    status_names <- dimnames(W)[[3]]
    if (is.null(status_names)) {
      stop("Character compStatus requires named condition levels in W.")
    }
    comp_idx <- match(compStatus, status_names)
    if (anyNA(comp_idx)) {
      stop("All compStatus values must match condition names in W.")
    }
    return(comp_idx)
  }

  comp_idx <- as.integer(compStatus)
  if (anyNA(comp_idx) || any(comp_idx < 1L | comp_idx > dim(W)[3])) {
    stop("compStatus indices must be between 1 and the number of conditions.")
  }

  comp_idx
}

weighted_means_from_weights <- function(Y, weights) {
  weight_sums <- colSums(weights)
  weighted_totals <- t(crossprod(weights, Y))
  sweep(weighted_totals, 2, weight_sums, "/")
}

DE_mu <- function(Y, W, compStatus) {
  comp_idx <- resolve_comp_status_indices(W, compStatus)
  rest_idx <- setdiff(seq_len(dim(W)[3]), comp_idx)

  M_out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]),
                 dimnames = list(colnames(Y), dimnames(W)[[2]], dimnames(W)[[3]]))

  for (status_idx in comp_idx) {
    M_out[, , status_idx] <- weighted_means_from_weights(Y, W[, , status_idx, drop = FALSE][, , 1])
  }

  if (length(rest_idx) > 0) {
    W_rest <- apply(W[, , rest_idx, drop = FALSE], c(1, 2), mean)
    rest_means <- weighted_means_from_weights(Y, W_rest)
    for (status_idx in rest_idx) {
      M_out[, , status_idx] <- rest_means
    }
  }

  M_out
}

#' Estimate Mean Expression Under Null Hypothesis (DE)
#'
#' Computes weighted mean expression matrices for the null hypothesis
#' in differential expression testing, where specified conditions share
#' a common mean expression.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)},
#'   where \eqn{C} is the number of cell types (or classes),
#'   and \eqn{V} is the number of conditions (or contrasts).
#' @param compStatus Integer or character vector specifying which
#'   conditions are being compared. These conditions will share a
#'   common mean estimate under the null hypothesis.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)}
#'   containing estimated mean expression values where means for
#'   conditions in \code{compStatus} are identical, and other
#'   conditions have separate estimates.
#'
#' @details
#' Under the null hypothesis, conditions in \code{compStatus} share the same
#' expression levels, so their means are estimated from combined data
#' to improve statistical power.
#'
#' @keywords internal
DE_mu_null <- function(Y, W, compStatus) {
  comp_idx <- resolve_comp_status_indices(W, compStatus)
  rest_idx <- setdiff(seq_len(dim(W)[3]), comp_idx)

  M_out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]),
                 dimnames = list(colnames(Y), dimnames(W)[[2]], dimnames(W)[[3]]))

  W_combined <- apply(W[, , comp_idx, drop = FALSE], c(1, 2), mean)
  shared_mean <- weighted_means_from_weights(Y, W_combined)
  for (status_idx in comp_idx) {
    M_out[, , status_idx] <- shared_mean
  }

  if (length(rest_idx) > 0) {
    W_rest <- apply(W[, , rest_idx, drop = FALSE], c(1, 2), mean)
    rest_means <- weighted_means_from_weights(Y, W_rest)
    for (status_idx in rest_idx) {
      M_out[, , status_idx] <- rest_means
    }
  }

  M_out
}

#' Estimate Variance per Gene, Cell Type, and Condition (DE)
#'
#' Computes weighted variance estimates for differential expression analysis.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)} containing
#'   posterior weights for each sample, cell type, and condition.
#' @param M A numeric array of dimension \eqn{(p \times C \times V)} containing
#'   mean estimates.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)} with
#'   weighted variance estimates. Minimum variance floor of 0.1 is applied.
#'
#' @keywords internal
DE_sigma2 <- function(Y, W, M) {
  sigma2 <- array(0, dim = dim(M))
  dimnames(sigma2) <- dimnames(M)

  Y_sq <- Y * Y

  for (status_idx in seq_len(dim(W)[3])) {
    weights <- W[, , status_idx, drop = FALSE][, , 1]
    weight_sums <- colSums(weights)

    first_moment <- t(crossprod(weights, Y))
    first_moment <- sweep(first_moment, 2, weight_sums, "/")

    second_moment <- t(crossprod(weights, Y_sq))
    second_moment <- sweep(second_moment, 2, weight_sums, "/")

    sigma2[, , status_idx] <- pmax(
      second_moment - 2 * M[, , status_idx] * first_moment + M[, , status_idx]^2,
      0.1
    )
  }
  return(sigma2)
}

#' Estimate Mixing Proportions (DE)
#'
#' Computes mixing proportions for differential expression analysis
#' as the average posterior weights.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)} containing
#'   posterior weights.
#'
#' @return A numeric matrix of dimension \eqn{(C \times V)} with mixing
#'   proportions for each cell type-condition pair.
#'
#' @keywords internal
DE_probs <- function(Y, W){
  out.probs <- apply(W, c(2,3), sum)/dim(Y)[1]
  dimnames(out.probs) <- dimnames(W)[2:3]
  return(out.probs)
}

de_lrt_exceeds_from_weights <- function(Y, W, compGroups, lrt) {
    M.alt <- DE_mu(Y = Y, W = W, compStatus = compGroups)
    sigma2.alt <- DE_sigma2(Y, W, M.alt)

    M.null <- DE_mu_null(Y, W, compGroups)
    sigma2.null <- DE_sigma2(Y, W, M.null)

    ll_pair <- approx_complete_data_loglik_pair_fast(Y, M.alt, sigma2.alt, M.null, sigma2.null, W)
    perm_lrt <- -2 * (ll_pair[["ll.null"]] - ll_pair[["ll.alternative"]])

    (perm_lrt >= lrt) * 1L
}



#' Compute Log-Fold Changes Between Conditions
#'
#' Computes log-fold change (logFC) contrasts between all pairs of
#' infection level stauses for each gene and partitioning variable leveles
#' (e.g., cell types), given estimated mean expression values.
#'
#' @param M A numeric array of dimension \eqn{(p \times C \times V)},
#'   where \eqn{p} is the number of genes, \eqn{C} is the number of
#'   partitioning variable levels, and \eqn{V} is the number of conditions. Entries
#'   represent estimated mean expression levels.
#' @param compStatus Character vector of the infection status levels to output 
#'   pairwise log fold changes. Default will return a list of all pairwise comparisons.
#'
#' @return A named list of data frames, one for each pairwise
#'   comparison. Each data frame has \eqn{p} rows
#'   (genes) and \eqn{C} columns (partitioning variable/cell types).
#'
#' @importFrom utils combn
#' @keywords internal
compute_logFC <- function(M, compStatus = NULL) {
    
    num_genes <- dim(M)[1]
    num_celltypes <- dim(M)[2]
    v_dim <- dimnames(M)[[3]]
    
    if (is.null(compStatus)) {
        combs <- combn(v_dim, 2)
    } else {
        combs <- combn(compStatus, 2)
    }
    
    logFC_list <- list()
    
    for (i in seq_len(ncol(combs))) {
        v1 <- combs[1, i]
        v2 <- combs[2, i]
        
        df <- M[, , v1] - M[, , v2]

        colnames(df) <- paste0("logFC_", colnames(df), "_", v1, "_vs_", v2)
        
        logFC_list[[paste0(v1, "_vs_", v2)]] <- as.data.frame(df)
    }
    
    if(length(logFC_list) == 1) logFC_list <- logFC_list[[1]]
    
    return(logFC_list)
}




