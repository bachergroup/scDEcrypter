#' E-step of Multiway Mixture Model
#'
#' @param Y Expression matrix (cell x gene).
#' @param c_obs Observed condition labels (vector).
#' @param v_obs Observed infection status labels (vector).
#' @param M Mean array (gene x condition x infection status).
#' @param probs Mixing proportions (condition x infection status matrix).
#' @param sigma2 Variance array (gene x condition x infection status).
#' @return Weight tensor (cell x condition x infection status).
#' @importFrom stats dnorm
#' @importFrom matrixStats logSumExp
#' @keywords internal
E_step <- function(Y, c_obs, v_obs, M, probs, sigma2){
  c_dim <- dim(M)[2]
  v_dim <- dim(M)[3]
  n_num_cells <- dim(Y)[1]
  p_num_genes <- dim(Y)[2]
  W <- array(0, dim=c(n_num_cells, c_dim, v_dim),
                dimnames = list(rownames(Y), NULL, NULL))
  dimnames(W)[2:3] <- dimnames(M)[2:3]
  genes <- dimnames(M)[[1]]
  
  for (kk in 1:n_num_cells) {
    log.numerator <- array(0, dim=c(c_dim, v_dim))
    for (cc in 1:c_dim) {
        for (vv in 1:v_dim) {
            log.numerator[cc, vv] <- sum(dnorm(Y[kk, genes],
                                  mean = M[, cc, vv],
                                  sd   = sqrt(sigma2[, cc, vv]),
                                  log  = T)) + log(probs[cc, vv])
        }
    }
    if(!is.na(c_obs[kk]) & is.na(v_obs[kk])){
        W[kk,c_obs[kk],] <- exp(log.numerator[c_obs[kk],] - logSumExp(log.numerator[c_obs[kk],],na.rm=T))
    }
    if(!is.na(v_obs[kk]) & is.na(c_obs[kk])){
        W[kk,,v_obs[kk]] <- exp(log.numerator[,v_obs[kk]] - logSumExp(log.numerator[,v_obs[kk]],na.rm=T))
    }
    if(!is.na(v_obs[kk]) & !is.na(c_obs[kk])){
        W[kk,c_obs[kk],v_obs[kk]] <- 1
    }
    if(is.na(c_obs[kk]) & is.na(v_obs)[kk]){
        W[kk,,] <- exp(log.numerator - logSumExp(log.numerator,na.rm=T))
    }
  }
  return(W)
}