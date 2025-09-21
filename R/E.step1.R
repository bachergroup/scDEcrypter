#' E-step of Multiway Mixture Model
#'
#' @param Y Expression matrix (cell x gene).
#' @param C.obs Observed condition labels (vector).
#' @param V.obs Observed infection status labels (vector).
#' @param M Mean array (gene x condition x infection status).
#' @param probs Mixing proportions (condition x infection status matrix).
#' @param sigma2 Variance array (gene x condition x infection status).
#' @return Weight tensor (cell x condition x infection status).
#' @export


E.step1 <- function(Y, C.obs, V.obs, M, probs, sigma2){
  C.dim <- dim(M)[2]
  V.dim <- dim(M)[3]
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  W <- array(0, dim=c(n, C.dim, V.dim))
  
  for (kk in 1:n) {
    log.numerator <- array(0, dim=c(C.dim, V.dim))
    for (vv in 1:C.dim) {
      for (v in 1:V.dim) {
        mu <- M[, vv, v]
        sigma2.tmp <- sigma2[, vv, v]
        # sigma2.tmp[sigma2.tmp == 0] <- 1e-10
        log_likelihood <- -(1/2) * sum((Y[kk, ] - mu)^2 / sigma2.tmp) - (1/2) * sum(log(2 * pi * sigma2.tmp))
        log.numerator[vv, v] <- log_likelihood + log(probs[vv,v])
      }
    }
    if(!is.na(C.obs[kk]) & is.na(V.obs[kk])){
      W[kk,C.obs[kk],] <- exp(log.numerator[C.obs[kk],] - logSumExp(log.numerator[C.obs[kk],]))
    }
    if(!is.na(V.obs[kk]) & is.na(C.obs[kk])){
      W[kk,,V.obs[kk]] <- exp(log.numerator[,V.obs[kk]] - logSumExp(log.numerator[,V.obs[kk]]))
    }
    if(!is.na(V.obs[kk]) & !is.na(C.obs[kk])){
      W[kk,C.obs[kk],V.obs[kk]] <- 1
    }
    if(is.na(C.obs[kk]) & is.na(V.obs)[kk]){
      W[kk,,] <- exp(log.numerator - logSumExp(log.numerator))
    }
  }
  return(W)
}