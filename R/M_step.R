
#' Update Mean Parameters With Penalty
#'
#' Updates the mean parameter array \eqn{M} in the multiway mixture model
#' using penalized accelerated proximal gradient descent (APGD).
#' This version applies a row-wise penalty controlled by \code{lambda}.
#'
#' @param Y Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param M Numeric array of current mean estimates with dimensions
#'   \code{(p x C x V)}, where \code{C} is the number of conditions
#'   and \code{V} is the number of infection status
#' @param sigma2 Numeric array of variance estimates with dimensions
#'   \code{(p x C x V)}.
#' @param W Numeric 3D array of posterior weights (responsibilities)
#'   with dimensions \code{(n x C x V)}.
#' @param lambda Numeric scalar; penalty parameter applied in the update.
#'
#' @return A numeric array of updated estimates.
#'
#' @export
update_mu <- function(Y, M, sigma2, W, lambda) {

  C_dim <- dim(M)[2]
  V_dim <- dim(M)[3]
  n <- dim(Y)[1]
  p <- dim(Y)[2]

  mu_j_mat <- array(0, dim = c(p, C_dim, V_dim))
  dimnames(mu_j_mat) <- dimnames(M)
  
  y_bar_jcv <- array(0, dim = c(p, C_dim, V_dim))
  epsilon_cv <- array(0, dim = c(p, C_dim, V_dim))
  for (j in seq_len(p)) { # for each gene
    for (c in seq_len(C_dim)) {
      for (v in seq_len(V_dim)) {
        y_j <- Y[, j] # per gene, all cells
        epsilon_icv <- 1 / W[, c, v] # weight per cell

        y_bar_jcv[j,c,v] <- sum(y_j / epsilon_icv) # store for gene
        epsilon_cv[j, c, v] <- sum(1 / epsilon_icv) # store for gene
      }
    }
    for (c in seq_len(C_dim)) {
      # print(c)
     if(all(epsilon_cv[j,c,] > 0)) {
      y_tilde_jc <- y_bar_jcv[j,c,] / epsilon_cv[j,c,]
      a <- as.vector(t(y_tilde_jc))
      S <- as.vector(t(sigma2[j,c,] / epsilon_cv[j,c,])) 
      # Store the updated mean values
      mu_j_updated <- AccPGD.Dm(a = a, S = S, lambda = lambda,
                                M.init = as.vector(t(M[j,c,])),
                                tol=1e-12)
     } else {
         # print("AHHH")
        mu_j_updated <- as.vector(t(M[j,c,]))
     }
      mu_j_mat[j,c , ] <- mu_j_updated
    }
  }

  return(mu_j_mat)
}

#' Update Mean Parameters (M-step, No Penalty)
#'
#' Computes updated mean estimates for a multiway mixture model
#' without penalty regularization.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param M A numeric array of dimension \eqn{(p \times C \times V)} containing
#'   current estimates of mean expression.
#' @param sigma2 A numeric array of dimension \eqn{(p \times C \times V)} containing
#'   variance estimates (not used in calculation but kept for function signature).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)} containing
#'   posterior weights for each sample, cell type, and condition.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)} with updated
#'   mean estimates computed as weighted averages.
#'
#' @details
#' Updates mean estimates as:
#' \eqn{M[j, c, v] = \frac{\sum_i W[i, c, v] * Y[i, j]}{\sum_i W[i, c, v]}}
#'
#' @keywords internal
update_mu_nopenalty <- function(Y, M, sigma2, W){
  M.out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]))
  dimnames(M.out) <- dimnames(M)
  W.tot <- apply(W, c(2,3), sum)
  W.tot[W.tot ==0] <- .1
  for(j in seq_len(dim(M)[1])){
    M.out[j,,] <- apply(W*Y[,j], c(2,3), sum)/W.tot
  }
  return(M.out)
}

#' Estimate Variance Parameters (M-step)
#'
#' Computes variance estimates for a multiway mixture model using
#' weighted squared deviations from current mean estimates.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)} containing
#'   posterior weights for each sample, cell type, and condition.
#' @param M A numeric array of dimension \eqn{(p \times C \times V)} containing
#'   current mean estimates.
#'
#' @return A numeric array of dimension \eqn{(p \times C \times V)} with updated
#'   variance estimates. Each entry is the weighted variance with a minimum floor of 0.1.
#'
#' @details
#' Computes variance as:
#' \eqn{\sigma^2[j, c, v] = \max\left(\frac{\sum_i W[i, c, v] * (Y[i, j] - M[j, c, v])^2}{\sum_i W[i, c, v]}, 0.1\right)}
#' 
#' The floor of 0.1 prevents numerical issues with very small variances.
#'
#' @keywords internal
M_step_variance <- function(Y, W, M){
	p <- dim(Y)[2]
	sigma2 <- array(0, dim=dim(M))
	dimnames(sigma2) <- dimnames(M)
	W.tot <- apply(W, c(2,3), sum)
	W.tot[W.tot ==0] <- .1
	for(kk in seq_len(p)){
	  for(c.ind in 1:dim(M)[2]){
	    for(v.ind in 1:dim(M)[3]){
	      sigma2[kk, c.ind, v.ind] <- max(sum(W[, c.ind, v.ind] * (Y[, kk] - M[kk, c.ind, v.ind])^2) / W.tot[c.ind, v.ind], .1)
	    }
	  }
	}
	return(sigma2)
	}

#' Estimate Mixing Proportions (M-step)
#'
#' Computes updated mixing proportions for cell type and condition
#' clusters as the average posterior weights.
#'
#' @param Y A numeric matrix of dimension \eqn{(n \times p)}, where
#'   \eqn{n} is the number of samples (cells) and \eqn{p} is the
#'   number of features (genes).
#' @param W A numeric array of dimension \eqn{(n \times C \times V)} containing
#'   posterior weights for each sample, cell type, and condition.
#'
#' @return A numeric matrix of dimension \eqn{(C \times V)} containing
#'   the mixing proportions (average posterior probabilities).
#'
#' @details
#' Computes the proportion for each cell type-condition pair as:
#' \eqn{\pi[c, v] = \frac{1}{n}\sum_i W[i, c, v]}
#'
#' @keywords internal
M_step_probs <- function(Y, W){
  apply(W, c(2,3), sum)/dim(Y)[1]
}

