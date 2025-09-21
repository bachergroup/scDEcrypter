

#' Estimate Mixture Probabilities for Differential Expression Model
#'
#' Computes the estimated mixture probabilities (class/condition proportions)
#' by averaging posterior weights across all samples.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes). This argument is only used to
#'   extract the number of samples \eqn{n}.
#' @param W A numeric array of dimension \eqn{(n \times C \times V)},
#'   where \eqn{C} is the number of cell types and \eqn{V} is the
#'   number of conditions. Each slice \code{W[i,,]} contains the
#'   posterior weights for sample \eqn{i}.
#'
#' @return A numeric matrix of dimension \eqn{(C \times V)} containing
#'   the estimated mixture probabilities. Each entry \code{probs[c, v]}
#'   represents the average posterior probability of cell type \eqn{c}
#'   under condition \eqn{v}.
#'
#' @export
DE.probs <- function(Y, W){
  apply(W, c(2,3), sum)/dim(Y)[1]
}
