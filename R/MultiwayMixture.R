

#' Multiway Mixture Model with Penalty
#'
#' Fits a multiway mixture model using the EM algorithm with
#' an \eqn{\ell_2}-type penalty on mean parameters to encourage
#' similarity across variants.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param C.obs Integer vector of length \code{n}, giving the observed
#'   condition labels for each sample (use \code{NA} if unknown).
#' @param V.obs Integer vector of length \code{n}, giving the observed
#'   variant labels for each sample (use \code{NA} if unknown).
#' @param max.iter Integer. Maximum number of EM iterations for each
#'   value of \code{lambda}.
#' @param tol Numeric. Convergence tolerance for relative change in
#'   mean estimates.
#' @param C.star Integer. Number of condition clusters to fit.
#' @param V.star Integer. Number of variant clusters to fit.
#' @param seed Integer. Random seed for reproducibility.
#' @param lambda.vec Numeric vector of penalty parameters. A separate
#'   model is fitted for each \code{lambda}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{M_list}}{List of estimated means (\code{p x C x V})
#'   for each \code{lambda}.}
#'   \item{\code{sigma2_list}}{List of estimated variances
#'   (\code{p x C x V}) for each \code{lambda}.}
#'   \item{\code{probs_list}}{List of estimated mixing proportions
#'   (\code{C x V}) for each \code{lambda}.}
#'   \item{\code{weights_list}}{List of posterior responsibilities
#'   (\code{n x C x V}) for each \code{lambda}.}
#' }
#'
#' @export
MultiwayMixture <- function(Y, C.obs, V.obs, max.iter, tol, C.star, V.star, seed, lambda.vec){

  set.seed(seed)
  M_lambda_list <- list()
  sigma2_lambda_list <- list()
  probs_lambda_list <- list()
  weights_lambda_list <- list()
  tmp <- MultiwayMixture_nopenalty(Y, C.obs, V.obs, max.iter, tol, C.star, V.star, seed)
  
  for(lambda in lambda.vec) {
    
    M <- tmp$M
    sigma2 <- tmp$sigma2
    probs <- tmp$probs
    
    for(mm in 1:max.iter){
      print("E.step")
      W1 <- E.step1(Y, C.obs, V.obs, M, probs, sigma2)
      print("M.step.prob")
      probs.new <- M.step.probs(Y, W1)
      print("update.mu")
      M.new <- update_mu(Y, M, sigma2, W1, lambda = lambda)
      print("M.step.variance")
      sigma2.new <- M.step.variance(Y, W1, M.new)
      print(mm)
      # log_likelihood <- observed.data.loglik(Y, M.new, sigma2.new, probs.new, C.obs, V.obs) - lambda*penalty.func(M.new)
      log_likelihood <- observed.data.loglik(Y, M.new, sigma2.new, probs.new, C.obs, V.obs) 
      
      cat(log_likelihood, "\n")
      cat(sum((M - M.new)^2)/sum(M^2), "\n")
      if(sum((M - M.new)^2)/sum(M^2) < tol){
        break
      }
      M <- M.new
      sigma2 <- sigma2.new
      probs <- probs.new
    }
    
    weights <- E.step1(Y, C.obs, V.obs, M, probs, sigma2)
    
    M_lambda_list[[as.character(lambda)]] <- M.new
    sigma2_lambda_list[[as.character(lambda)]] <- sigma2.new
    probs_lambda_list[[as.character(lambda)]] <- probs.new
    weights_lambda_list[[as.character(lambda)]] <- weights
  }
  
  return(list(
    "M_list" = M_lambda_list,
    "sigma2_list" = sigma2_lambda_list, 
    "probs_list" = probs_lambda_list,
    "weights_list" = weights_lambda_list
  ))
}