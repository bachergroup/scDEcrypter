
#' Observed-Data Log-Likelihood
#'
#' Computes the observed-data log-likelihood for the multiway mixture model,
#' given the observed data, mean estimates, variance estimates, mixing proportions,
#' and partially observed labels.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param M Numeric array of mean estimates with dimensions
#'   \code{(p x C x V)}, where \code{C} is the number of conditions
#'   and \code{V} is the number of infection status
#' @param sigma2 Numeric array of variance estimates with dimensions
#'   \code{(p x C x V)}.
#' @param probs Numeric matrix of mixing proportions with dimensions
#'   \code{(C x V)}.
#' @param C.obs Integer vector of length \code{n}, giving the observed
#'   condition labels for each sample (use \code{NA} if unknown).
#' @param V.obs Integer vector of length \code{n}, giving the observed
#'   variant labels for each sample (use \code{NA} if unknown).
#'
#' @return A numeric scalar: the observed-data log-likelihood value.
#'
#' @details
#' The function partitions samples based on whether their condition (\code{C.obs})
#' and variant (\code{V.obs}) labels are observed or missing:
#' \itemize{
#'   \item Both known: exact log-likelihood contribution.
#'   \item Condition unknown: marginalizes over all conditions.
#'   \item Variant unknown: marginalizes over all variants.
#'   \item Both unknown: marginalizes over both dimensions.
#' }
#'  @export

observed.data.loglik <- function(Y, M, sigma2, probs, C.obs, V.obs){

  n <- dim(Y)[1]; p <- dim(Y)[2]
  C.dim <- dim(M)[2]; V.dim <- dim(M)[3]
  cond.probs.row <- probs/rowSums(probs)
  cond.probs.col <- sweep(probs, 2, colSums(probs), FUN = "/")
  tot.contrib.vec <- rep(0, n)
  
  tot.contrib <- NULL; 
  both_unknown_contrib = 0
  v_unknown_contrib = 0
  c_unknown_contrib = 0
  both_unknown_contrib = 0
  
  # V observed and C observed
  both_known <- which(!is.na(V.obs) & !is.na(C.obs)) 
  if (length(both_known) > 0) {
    both_known_contrib <- sapply(both_known, function(index) {
      sum(dnorm(Y[index, ], mean = M[, C.obs[index], V.obs[index]],
                sd = sqrt(sigma2[, C.obs[index], V.obs[index]]),
                log = TRUE)) + log(probs[C.obs[index], V.obs[index]])
    })
  }
  
  # C unknown
  c_unknown <- which(!is.na(V.obs) & is.na(C.obs))
  if (length(c_unknown) > 0) {
    combos <- expand.grid(i = c_unknown, cc = 1:C.dim)
    combos$c_unknown_contrib <- mapply(function(i, cc) {
      logSumExp(dnorm(Y[i, ], mean = M[, cc, V.obs[i]],
                      sd = sqrt(sigma2[, cc, V.obs[i]]), log = TRUE)) +
        log(cond.probs.col[cc, V.obs[i]])}, combos$i, combos$cc)
    c_unknown_contrib <- tapply(combos$c_unknown_contrib, combos$i, logSumExp)
  }
  
  # V unknown
  v_unknown <- which(is.na(V.obs) & !is.na(C.obs))
  if (length(v_unknown) > 0) {
    combos <- expand.grid(i = v_unknown, v = 1:V.dim)
    combos$v_unknown_contrib <- mapply(function(i, v) {
      sum(dnorm(Y[i, ], mean = M[, C.obs[i], v],
                sd = sqrt(sigma2[, C.obs[i], v]), log = TRUE)) +
        log(cond.probs.row[C.obs[i], v])}, combos$i, combos$v)
    v_unknown_contrib <- tapply(combos$v_unknown_contrib, combos$i, logSumExp)
  }
  
  # Both unknown
  both_unknown <- which(is.na(V.obs) & is.na(C.obs))
  if (length(both_unknown) > 0) {
    combos <- expand.grid(i = both_unknown, v = 1:V.dim, cc=1:C.dim)
    combos$both_unknown_contrib <- mapply(function(i, v, cc) {
      logSumExp(dnorm(Y[i, ], mean = M[, cc, v],
                      sd = sqrt(sigma2[, cc, v]), log = TRUE)) +
        log(probs[cc, v])}, combos$i, combos$v, combos$cc)
    both_unknown_contrib <- tapply(combos$both_unknown_contrib, combos$i, logSumExp)
  }
  
  # Sum together
  tot.contrib <- sum(c(both_unknown_contrib, v_unknown_contrib, 
                       c_unknown_contrib, both_known_contrib))
  
  return(tot.contrib)
}
