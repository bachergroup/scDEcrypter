
#' Estimate Variance for Differential Expression Model
#'
#' Computes the weighted variance of gene expression across samples
#' for each gene, cell type, and condition, given posterior weights
#' and mean estimates.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)},
#'   where \eqn{C} is the number of cell types and \eqn{V} is the
#'   number of conditions. Each slice \code{W[i,,]} contains the
#'   posterior weights for sample \eqn{i}.
#' @param M A numeric array of dimension \eqn{(p \times C \times V)}
#'   representing the estimated mean expression values.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)},
#'   where each entry \code{sigma2[j, c, v]} is the weighted variance
#'   of gene \eqn{j} for cell type \eqn{c} under condition \eqn{v}.
#'
#' @export
DE.sigma2 <- function(Y, W, M){
  p <- dim(Y)[2]
  sigma2 <- array(0, dim=dim(M))
  W.tot <- apply(W, c(2,3), sum)
  for(kk in 1:p){
    for(c.ind in 1:dim(M)[2]){
      for(v.ind in 1:dim(M)[3]){
        sigma2[kk, c.ind, v.ind] <- max(sum(W[, c.ind, v.ind] * (Y[, kk] - M[kk, c.ind, v.ind])^2) / W.tot[c.ind, v.ind], 1e-10)
      }
    }
  }
  return(sigma2)
}
