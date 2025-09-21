
#' Update Variance Parameters (M-step)
#'
#' Computes updated variance estimates \eqn{\sigma^2} in the M-step
#' of the multiway mixture model EM algorithm.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param W Numeric 3D array of posterior weights (responsibilities)
#'   with dimensions \code{(n x C x V)}, where \code{C} is the number
#'   of conditions and \code{V} is the number of infection status
#' @param M Numeric array of mean estimates with dimensions
#'   \code{(p x C x V)}.
#'
#' @return A numeric array of updated variance estimates with dimensions
#'   \code{(p x C x V)}.
#'
#' @export
M.step.variance <- function(Y, W, M){
  p <- dim(Y)[2]
  sigma2 <- array(0, dim=dim(M))
  W.tot <- apply(W, c(2,3), sum)
  for(kk in 1:p){
    for(c.ind in 1:dim(M)[2]){
      for(v.ind in 1:dim(M)[3]){
        #weighted_mean <- sum(Y[,kk]*W[,c.ind,v.ind])/W.tot[c.ind, v.ind]
        sigma2[kk, c.ind, v.ind] <- max(sum(W[, c.ind, v.ind] * (Y[, kk] - M[kk, c.ind, v.ind])^2) / W.tot[c.ind, v.ind], 1e-10)
      }
    }
  }
  return(sigma2)
}
