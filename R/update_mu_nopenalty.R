
#' Update Mean Parameters Without Penalty
#'
#' Computes updated estimates of the mean parameter array \eqn{M} in the
#' multiway mixture model without applying a penalty.
#' The update is based on weighted averages using the responsibility tensor \eqn{W}.
#'
#' @param Y Numeric matrix of observed data, with rows representing samples (cells)
#'   and columns representing features (e.g., genes).
#' @param M Numeric array of current mean estimates with dimensions
#'   \code{(p x C x V)}, where \code{p} is the number of features,
#'   \code{C} the number of conditions, and \code{V} the number of infection status
#' @param sigma2 Numeric array of variance estimates (same dimensions as \code{M}).
#'   (Note: not used in the current function but kept for compatibility.)
#' @param W Numeric 3D array of posterior weights (responsibilities) with dimensions
#'   \code{(n x C x V)}, where \code{n} is the number of samples.
#'
#' @return A numeric array of updated mean estimates with dimensions
#'   \code{(p x C x V)}.
#'
#' @export

update_mu_nopenalty <- function(Y, M, sigma2, W){
  M.out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]))
  for(j in 1:dim(M)[1]){
    M.out[j,,] <- apply(W*Y[,j], c(2,3), sum)/apply(W, c(2,3), sum)
  }
  return(M.out)
}
