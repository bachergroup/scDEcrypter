

#' Multiway Mixture Model (No Penalty)
#'
#' Fits a multiway mixture model using the EM algorithm without
#' applying any penalty term on the mean estimates.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param C.obs Integer vector of length \code{n}, giving the observed
#'   condition labels for each sample (use \code{NA} if unknown).
#' @param V.obs Integer vector of length \code{n}, giving the observed
#'   infection status labels for each sample (use \code{NA} if unknown).
#' @param max.iter Integer. Maximum number of EM iterations.
#' @param tol Numeric. Convergence tolerance for relative change in
#'   mean estimates.
#' @param C.star Integer. Number of condition clusters to fit.
#' @param V.star Integer. Number of infection status clusters to fit.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{M}}{Array of estimated means (\code{p x C x V}).}
#'   \item{\code{sigma2}}{Array of estimated variances (\code{p x C x V}).}
#'   \item{\code{probs}}{Matrix of estimated mixing proportions (\code{C x V}).}
#'   \item{\code{weights}}{Array of posterior responsibilities (\code{n x C x V}).}
#' }
#'
#'  @export
MultiwayMixture_nopenalty <- function(Y, C.obs, V.obs, max.iter, tol, C.star, V.star, seed){

  set.seed(seed)
  M <- array(0, dim=c(dim(Y)[2], C.star, V.star))
  sigma2 <- array(1, dim=c(dim(Y)[2], C.star, V.star))
  probs <- matrix(1/(C.star*V.star), C.star, V.star)
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  C.levels <- sort(unique(na.omit(C.obs)))
  idxC <- which(!is.na(C.obs)) # cells with known C
  idxNoC <- which(is.na(C.obs)) # cells with unknown C
  
  # mean by cell type when C.obs not NA
  mean_by_C <- function(rows) {
    keep <- intersect(rows, idxC)
    Mmat <- sapply(seq_len(p), function(j) {
      tapply(Y[keep, j],
             factor(C.obs[keep], levels = C.levels),
             mean, na.rm = TRUE)
    })
    if (is.vector(Mmat)) Mmat <- matrix(Mmat, nrow = length(C.levels))
    rownames(Mmat) <- C.levels
    colnames(Mmat) <- colnames(Y)
    as.matrix(Mmat)
  }
  
  # Overall mean by C.obs (ignoring V)
  cmeans_byC <- mean_by_C(idxC)           
  
  # Global mean (ignoring C.obs and V.obs)
  g_vec  <- colMeans(Y, na.rm = TRUE)   
  
  fill_na_hier <- function(A, overall_byC, byV_vec, global_vec) {
    # 1. Fill v1_byC or v2_byC if both V and C known
    # 2. Fill overall mean by C.obs when V unknown
    na_idx <- is.na(A)
    if (any(na_idx)) A[na_idx] <- overall_byC[na_idx]
    # 3. Fill overall mean by V.obs when C unknown
    if (anyNA(A)) {
      na_rows <- which(apply(A, 1, function(r) any(is.na(r))))
      for (r in na_rows) {
        miss <- which(is.na(A[r, ]))
        if (length(miss) > 0) A[r, miss] <- byV_vec[miss]
      }
    }
    # 4. Fill overall mean of Y when C and V unknown
    if (anyNA(A)) {
      na2 <- which(is.na(A), arr.ind = TRUE)
      A[na2] <- global_vec[na2[, 2]]
    }
    A
  }
  
  V.levels <- 1:V.star
  for (vi in seq_along(V.levels)) {
    v <- V.levels[vi]
    idxV <- which(V.obs == v)
    
    v_byC <- mean_by_C(idxV)
    v_vec <- colMeans(Y[idxV, , drop=FALSE], na.rm=TRUE)
    v_byC_filled <- fill_na_hier(v_byC, cmeans_byC, v_vec, g_vec)
    
    M[, , vi] <- t(v_byC_filled)
  }
  
  for(mm in 1:max.iter){
    W1 <- E.step1(Y, C.obs, V.obs, M, probs, sigma2)
    probs.new <- M.step.probs(Y, W1)
    M.new <- update_mu_nopenalty(Y, M, sigma2, W1)
    sigma2.new <- M.step.variance(Y, W1, M.new)
    log_likelihood00 <- observed.data.loglik(Y = Y, M = M.new, sigma2 = sigma2.new, probs  = probs.new, C.obs = C.obs, V.obs = V.obs)
    cat(log_likelihood00, "\n")
    cat(sum((M - M.new)^2)/sum(M^2), "\n")
    if(sum((M - M.new)^2)/sum(M^2) < tol){
      break
    }
    M <- M.new
    sigma2 <- sigma2.new
    probs <- probs.new
  }
  
  weights <- E.step1(Y, C.obs, V.obs, M, probs, sigma2)
  
  return(list(
    "M" = M.new,
    "sigma2" = sigma2.new, 
    "probs" = probs.new,
    "weights" = weights
  ))
}

