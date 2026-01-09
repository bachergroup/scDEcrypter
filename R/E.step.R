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
  n_num_cells <- dim(Y)[1]
  p_num_genes <- dim(Y)[2]
  W <- array(0, dim=c(n_num_cells, C.dim, V.dim))

  for (kk in 1:n_num_cells) {
    log.numerator <- array(0, dim=c(C.dim, V.dim))
    for (cc in 1:C.dim) {
        for (vv in 1:V.dim) {
            log.numerator[cc, vv] <- sum(dnorm(Y[kk, ],
                                  mean = M[, cc, vv],
                                  sd   = sqrt(sigma2[, cc, vv]),
                                  log  = T)) + log(probs[cc, vv])
        }
    }
    if(!is.na(C.obs[kk]) & is.na(V.obs[kk])){
        W[kk,C.obs[kk],] <- exp(log.numerator[C.obs[kk],] - logSumExp(log.numerator[C.obs[kk],],na.rm=T))
    }
    if(!is.na(V.obs[kk]) & is.na(C.obs[kk])){
        W[kk,,V.obs[kk]] <- exp(log.numerator[,V.obs[kk]] - logSumExp(log.numerator[,V.obs[kk]],na.rm=T))
    }
    if(!is.na(V.obs[kk]) & !is.na(C.obs[kk])){
        W[kk,C.obs[kk],V.obs[kk]] <- 1
    }
    if(is.na(C.obs[kk]) & is.na(V.obs)[kk]){
        W[kk,,] <- exp(log.numerator - logSumExp(log.numerator,na.rm=T))
    }
  }
  return(W)
}