
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
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)}
#'   containing the weighted mean expression values.
#'   Specifically, \code{M.out[j, c, v]} is the estimated mean of
#'   gene \eqn{j} in cell type \eqn{c} under condition \eqn{v}.
#'
#' @export
DE.mu <- function(Y, W) {
  M.out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]))
  for (j in 1:dim(Y)[2]) {
    M.out[j,,] <- apply(W * Y[, j], c(2, 3), sum) / apply(W, c(2, 3), sum)
  }
  return(M.out)
}
