
#' Estimate Null Model Mean Expression
#'
#' Computes the mean expression array under a null model, where
#' condition/variant-specific weights are averaged across conditions.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)},
#'   where \eqn{C} is the number of cell types (or classes),
#'   and \eqn{V} is the number of conditions (or contrasts).
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)}
#'   containing the estimated mean expression values under the null model.
#'   Specifically, the condition weights are replaced by their average
#'   across \eqn{V}.
#'
#' @export
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
