
#' Estimate Mean Expression per Gene, Cell Type, and Condition
#'
#' Computes the weighted mean expression matrix across genes,
#' cell types, and conditions given observed data and posterior
#' weights.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)},
#'   where \eqn{C} is the number of cell types (or classes),
#'   and \eqn{V} is the number of conditions (or contrasts).
#'   Each slice \code{W[i,,]} contains the posterior weights for
#'   sample \eqn{i}.
#' @param M A numeric array of dimension \eqn{(p \times C \times V)}
#'   representing the estimated mean expression values.
#' @param target_status list of groups.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)}
#'   containing the weighted mean expression values.
#'   Specifically, \code{M.out[j, c, v]} is the estimated mean of
#'   gene \eqn{j} in cell type \eqn{c} under condition \eqn{v}.
#' 
#' @return A numeric array of dimension \eqn{(p \times C \times V)},
#'   where each entry \code{sigma2[j, c, v]} is the weighted variance
#'   of gene \eqn{j} for cell type \eqn{c} under condition \eqn{v}.
#'   
#' @return A numeric matrix of dimension \eqn{(C \times V)} containing
#'   the estimated mixture probabilities. Each entry \code{probs[c, v]}
#'   represents the average posterior probability of cell type \eqn{c}
#'   under condition \eqn{v}.
#'
#' @export
DE.mu <- function(Y, W, target_status) {
  V_dim <- dim(W)[3]
  M_out <- array(0, c(dim(Y)[2], dim(W)[2], V_dim))
  
  for (j in 1:dim(Y)[2]) {
    M_out[j, , target_status] <- apply(W[, , target_status] * Y[, j], c(2, 3), sum) /
      apply(W[, , target_status], c(2, 3), sum)
  }
  
  rest_status <- setdiff(1:V_dim, target_status)
  if (length(rest_status) > 0) {
    W_rest <- apply(W[, , rest_status, drop = FALSE], c(1, 2), mean)
    W_rest <- array(rep(W_rest, length(rest_status)),
                    dim = c(nrow(W_rest), ncol(W_rest), length(rest_status)))
    
    denom_rest <- apply(W_rest, c(2, 3), sum)
    for (j in 1:ncol(Y)) {
      M_out[j, , rest_status] <- apply(W_rest * Y[, j], c(2, 3), sum) / denom_rest
    }
  }
  
  return(M_out)
}
DE.mu.null <- function(Y, W, target_status) {
  V_dim <- dim(W)[3]
  M_out <- array(0, c(dim(Y)[2], dim(W)[2], V_dim))

  W_combined <- apply(W[, , target_status, drop = FALSE], c(1, 2), mean)

  W_combined_array <- array(rep(W_combined, length(target_status)),
                            dim = c(nrow(W_combined), ncol(W_combined), length(target_status)))

  denom_combined <- apply(W_combined_array, c(2, 3), sum)

  for (j in 1:ncol(Y)) {
    shared_mean <- apply(W_combined_array * Y[, j], c(2, 3), sum) / denom_combined
    M_out[j, , target_status] <- shared_mean
  }

  rest_status <- setdiff(1:V_dim, target_status)
  if (length(rest_status) > 0) {
    W_rest <- apply(W[, , rest_status, drop = FALSE], c(1, 2), mean)
    W_rest <- array(rep(W_rest, length(rest_status)),
                    dim = c(nrow(W_rest), ncol(W_rest), length(rest_status)))

    denom_rest <- apply(W_rest, c(2, 3), sum)
    for (j in 1:ncol(Y)) {
      M_out[j, , rest_status] <- apply(W_rest * Y[, j], c(2, 3), sum) / denom_rest
    }
  }

  return(M_out)
}

DE.sigma2 <- function(Y, W, M){
  p <- dim(Y)[2]
  sigma2 <- array(0, dim=dim(M))
  W.tot <- apply(W, c(2,3), sum)
  for(kk in 1:p){
    for(c.ind in 1:dim(M)[2]){
      for(v.ind in 1:dim(M)[3]){
        sigma2[kk, c.ind, v.ind] <- max(sum(W[, c.ind, v.ind] * (Y[, kk] - M[kk, c.ind, v.ind])^2) / W.tot[c.ind, v.ind], .1)      
      }
    }
  }
  return(sigma2)
}

DE.probs <- function(Y, W){
  apply(W, c(2,3), sum)/dim(Y)[1]
}




