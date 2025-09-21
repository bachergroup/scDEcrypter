
#' Update Mixing Proportions (M-step)
#'
#' Computes updated mixing proportions (class probabilities)
#' in the M-step of the multiway mixture model EM algorithm.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param W Numeric 3D array of posterior weights (responsibilities)
#'   with dimensions \code{(n x C x V)}, where \code{C} is the number
#'   of conditions and \code{V} is the number of infection status
#'
#' @return A numeric matrix of updated mixing proportions with dimensions
#'   \code{(C x V)}. Each entry corresponds to the expected proportion
#'   of samples assigned to condition \code{c} and variant \code{v}.
#'
#' @export
M.step.probs <- function(Y, W){
  apply(W, c(2,3), sum)/dim(Y)[1]
}
