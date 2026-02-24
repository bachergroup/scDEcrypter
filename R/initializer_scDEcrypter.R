

#' Multiway Mixture Model (No Penalty)
#'
#' Fits a multiway mixture model using the EM algorithm without
#' applying any penalty term on the mean estimates.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param c_obs Integer vector of length \code{n}, giving the observed
#'   condition labels for each sample (use \code{NA} if unknown).
#' @param v_obs Integer vector of length \code{n}, giving the observed
#'   infection status labels for each sample (use \code{NA} if unknown).
#' @param max.iter Integer. Maximum number of EM iterations.
#' @param tol Numeric. Convergence tolerance for relative change in
#'   mean estimates.
#' @param c_star Integer. Number of condition clusters to fit.
#' @param v_star Integer. Number of infection status clusters to fit.
#' @importFrom stats na.omit
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{M}}{Array of estimated means (\code{p x C x V}).}
#'   \item{\code{sigma2}}{Array of estimated variances (\code{p x C x V}).}
#'   \item{\code{probs}}{Matrix of estimated mixing proportions (\code{C x V}).}
#'   \item{\code{weights}}{Array of posterior responsibilities (\code{n x C x V}).}
#' }
#'
#' @keywords internal
initializer_scDEcrypter <- function(Y, c_obs, v_obs, max.iter, tol, c_star, v_star){

  M <- array(0, dim=c(dim(Y)[2], c_star, v_star))
  sigma2 <- array(1, dim=c(dim(Y)[2], c_star, v_star))
  probs <- matrix(1/(c_star*v_star), c_star, v_star)
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  c_levels <- 1:c_star #sort(unique(na.omit(c_obs)))
  idx_c <- which(!is.na(c_obs)) # cells with known c
  idx_no_c <- which(is.na(c_obs)) # cells with unknown c
  
  # mean by cell type when c_obs not NA
  mean_by_c <- function(rows) {
    keep <- intersect(rows, idx_c)
    f_obs_fac <- factor(c_obs, levels = c_levels)
    Mmat <- sapply(seq_len(p), function(j) {
      tapply(Y[keep, j],
             f_obs_fac[keep],
             mean, na.rm = TRUE)
    })
    if (is.vector(Mmat)) Mmat <- matrix(Mmat, nrow = length(c_levels))
    rownames(Mmat) <- c_levels
    colnames(Mmat) <- colnames(Y)
    as.matrix(Mmat)
  }
  
  # Overall mean by c_obs (ignoring v)
  cmeans_by_c <- mean_by_c(idx_c)        
  
  # Global mean (ignoring c_obs and v_obs)
  g_vec  <- colMeans(Y, na.rm = TRUE)   
  
  fill_na_hier <- function(A, overall_by_c, by_v_vec, global_vec) {
    # 1. Fill v1_by_c or v2_by_c if both v and c known
    # 2. Fill overall mean by c_obs when v unknown
    na_idx <- is.na(A)
    if (any(na_idx)) A[na_idx] <- overall_by_c[na_idx]
    # 3. Fill overall mean by v_obs when c unknown
    if (anyNA(A)) {
      na_rows <- which(apply(A, 1, function(r) any(is.na(r))))
      for (r in na_rows) {
        miss <- which(is.na(A[r, ]))
        if (length(miss) > 0) A[r, miss] <- by_v_vec[miss]
      }
    }
    # 4. Fill overall mean of Y when c and v unknown
    if (anyNA(A)) {
      na2 <- which(is.na(A), arr.ind = TRUE)
      A[na2] <- global_vec[na2[, 2]]
    }
    A
  }
  
  for (vii in seq_len(v_star)) {
    idx_v <- which(v_obs == vii)
    
    v_by_c <- mean_by_c(idx_v)
    v_vec <- colMeans(Y[idx_v, , drop=FALSE], na.rm=TRUE)
    v_by_c_filled <- fill_na_hier(v_by_c, cmeans_by_c, v_vec, g_vec)
    
    M[, , vii] <- t(v_by_c_filled)
    dimnames(M)[[1]] <- colnames(v_by_c_filled)
  }
  
  dimnames(sigma2)[[1]] <- dimnames(M)[[1]] 
  for(mm in seq_len(max.iter)){
    W1 <- E_step(Y, c_obs, v_obs, M, probs, sigma2)
    probs.new <- M_step_probs(Y, W1)
    M.new <- update_mu_nopenalty(Y, M, sigma2, W1)
    sigma2.new <- M_step_variance(Y, W1, M.new)
    if(sum((M - M.new)^2)/sum(M^2) < tol){
      break
    }
    M <- M.new
    sigma2 <- sigma2.new
    probs <- probs.new
  }
  
  rtrn.weights <- E_step(Y, c_obs, v_obs, M, probs, sigma2)
  
  return(list(
    "M" = M.new,
    "sigma2" = sigma2.new, 
    "probs" = probs.new,
    "weights" = rtrn.weights
  ))
}

