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
#' @param c_obs Integer vector of length matching the number of cells, giving the 
#'   observed cell type labels for each cell (may contain \code{NA}).
#' @param v_obs Integer vector of length matching the number of cells, giving the 
#'   observed viral type labels for each cell (may contain \code{NA}).
#'
#' @return A numeric matrix with dimensions (genes x cell types) containing the 
#'   approximate log-likelihood contribution for each gene in each cell type, 
#'   aggregated across viral types and weighted by the assignment weights.
#'
#' @importFrom stats dnorm
#' @keywords internal
approx_complete_data_loglik <- function (Y, M, W, sigma2, c_obs, v_obs){

    n_genes <- dim(Y)[2]
    c_dim <- dim(M)[2]
    v_dim <- dim(M)[3]
    
    if (length(unique(c_obs[!is.na(c_obs)])) != c_dim) {
        warning(paste0(
            "Number of unique c_obs labels does not match that in the mean matrix."
        ))
    }
    if (length(unique(v_obs[!is.na(v_obs)])) != v_dim) {
        warning(paste0(
            "Number of unique v_obs labels does not match that in the mean matrix."
        ))
    }
    
    ll_gene_celltype <- matrix(0, nrow = n_genes, ncol = c_dim)
    
    for (cc in seq_len(c_dim)) {
      tot_contrib_1 <- array(0, dim=c(n_genes, v_dim))
      log_likelihood_kc <- NULL
      for(kk in seq_len(n_genes)){
          for(vv in seq_len(v_dim)) {
              all_cells_ll <- dnorm(Y[, kk],
                                    mean = M[kk, cc, vv],
                                    sd   = sqrt(sigma2[kk, cc, vv]),
                                    log  = TRUE) * W[, cc, vv]
              all_cells_ll <- all_cells_ll[is.finite(all_cells_ll)]
              tot_contrib_1[kk,vv] <- sum(all_cells_ll)
           }
      }
      ll_gene_celltype[, cc] <- rowSums(tot_contrib_1)
    }
    rownames(ll_gene_celltype) <- dimnames(M)[[1]]
    colnames(ll_gene_celltype) <- dimnames(M)[[2]]
    return(ll_gene_celltype)
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
DE_mu <- function(Y, W, compStatus) {
 
    M_out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]),
                   dimnames = list(colnames(Y), dimnames(W)[[2]], dimnames(W)[[3]]))
    
    for (j in seq_len(dim(Y)[2])) {
        M_out[j, , compStatus] <- apply(W[, , compStatus] * Y[, j], c(2, 3), sum) /
            apply(W[, , compStatus], c(2, 3), sum)
    }
    
    rest_status <- setdiff(dimnames(W)[[3]], compStatus)
    if (length(rest_status) > 0) {
        W_rest <- apply(W[, , rest_status, drop = FALSE], c(1, 2), mean)
        W_rest <- array(rep(W_rest, length(rest_status)),
                        dim = c(nrow(W_rest), ncol(W_rest), length(rest_status)))
        
        denom_rest <- apply(W_rest, c(2, 3), sum)
        for (j in seq_len(ncol(Y))) {
            M_out[j, , rest_status] <- apply(W_rest * Y[, j], c(2, 3), sum) / denom_rest
        }
    }
    
    return(M_out)
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
  V_dim <- dim(W)[3]
  M_out <- array(0, c(dim(Y)[2], dim(W)[2], V_dim),
                   dimnames = list(colnames(Y), dimnames(W)[[2]], dimnames(W)[[3]]))

  W_combined <- apply(W[, , compStatus, drop = FALSE], c(1, 2), mean)

  W_combined_array <- array(rep(W_combined, length(compStatus)),
                            dim = c(nrow(W_combined), ncol(W_combined), length(compStatus)))

  denom_combined <- apply(W_combined_array, c(2, 3), sum)

  for (j in seq_len(ncol(Y))) {
    shared_mean <- apply(W_combined_array * Y[, j], c(2, 3), sum) / denom_combined
    M_out[j, , compStatus] <- shared_mean
  }

  rest_status <- setdiff(seq_len(V_dim), compStatus)
  if (length(rest_status) > 0) {
    W_rest <- apply(W[, , rest_status, drop = FALSE], c(1, 2), mean)
    W_rest <- array(rep(W_rest, length(rest_status)),
                    dim = c(nrow(W_rest), ncol(W_rest), length(rest_status)))

    denom_rest <- apply(W_rest, c(2, 3), sum)
    for (j in seq_len(ncol(Y))) {
      M_out[j, , rest_status] <- apply(W_rest * Y[, j], c(2, 3), sum) / denom_rest
    }
  }

  return(M_out)
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
DE_sigma2 <- function(Y, W, M){
  p <- dim(Y)[2]
  sigma2 <- array(0, dim=dim(M))
  dimnames(sigma2) <- dimnames(M)
  W.tot <- apply(W, c(2,3), sum)
  for(kk in seq_len(p)){
    for(c.ind in seq_len(dim(M)[2])){
      for(v.ind in seq_len(dim(M)[3])){
        sigma2[kk, c.ind, v.ind] <- max(sum(W[, c.ind, v.ind] * (Y[, kk] - M[kk, c.ind, v.ind])^2) / W.tot[c.ind, v.ind], .1)      
      }
    }
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




