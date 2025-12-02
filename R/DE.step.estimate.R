
#' Estimate Mean Expression per Gene, Cell Type, and Condition
#'
#' Computes the weighted mean expression matrix across genes,
#' cell types, and conditions given observed data and posterior
#' weights.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)},
#'   where \eqn{C} is the number of cell types (or classes),
#'   and \eqn{V} is the number of conditions (or contrasts).
#'   Each slice \code{W[i,,]} contains the posterior weights for
#'   sample \eqn{i}.
#' @param M A numeric array of dimension \eqn{(p \times C \times V)}
#'   representing the estimated mean expression values.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)}
#'   containing the weighted mean expression values.
#'   Specifically, \code{M.out[j, c, v]} is the estimated mean of
#'   gene \eqn{j} in cell type \eqn{c} under condition \eqn{v}.
#' 
#' @return A numeric array of dimension \eqn{(p \times C \times V)},
#'   where each entry \code{sigma2[j, c, v]} is the weighted variance
#'   of gene \eqn{j} for cell type \eqn{c} under condition \eqn{v}.
#'   
#' @return A numeric matrix of dimension \eqn{(C \times V)} containing
#'   the estimated mixture probabilities. Each entry \code{probs[c, v]}
#'   represents the average posterior probability of cell type \eqn{c}
#'   under condition \eqn{v}.
#'
#' @export
DE.mu <- function(Y, W) {
  M.out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]))
  for (j in 1:dim(Y)[2]) {
    M.out[j,,] <- apply(W * Y[, j], c(2, 3), sum) / apply(W, c(2, 3), sum)
  }
  return(M.out)
}

DE.mu.null <- function (Y, W)  {
  V_dim <- dim(W)[3]
  M_out <- array(0, c(dim(Y)[2], dim(W)[2], V_dim))
  
  W_combined <- apply(W, c(1, 2), mean)
  W_combined <- array(rep(W_combined, V_dim),
                      dim = c(nrow(W_combined),
                              ncol(W_combined),
                              V_dim))
  
  denominator1 <- apply(W_combined, c(2, 3), sum)
  for (j in 1:ncol(Y)) {
    M_out[j, , ] <- apply(W_combined * Y[, j], c(2, 3), sum) / denominator1
  }
  return(M_out)
}

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

DE.probs <- function(Y, W){
  apply(W, c(2,3), sum)/dim(Y)[1]
}




