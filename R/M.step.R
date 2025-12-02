
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

  y_bar_jcv <- array(0, dim = c(p, C_dim, V_dim))
  epsilon_cv <- array(0, dim = c(p, C_dim, V_dim))
  for (j in 1:p) { # for each gene
    for (c in 1:C_dim) {
      for (v in 1:V_dim) {
        y_j <- Y[, j] # per gene, all cells
        epsilon_icv <- 1 / W[, c, v] # weight per cell

        y_bar_jcv[j,c,v] <- sum(y_j / epsilon_icv) # store for gene
        epsilon_cv[j, c, v] <- sum(1 / epsilon_icv) # store for gene
      }
    }
    for (c in 1:C_dim) {
      # print(c)
      y_tilde_jc <- y_bar_jcv[j,c,] / epsilon_cv[j,c,]
      a <- as.vector(t(y_tilde_jc))
      S <- as.vector(t(sigma2[j,c,] / epsilon_cv[j,c,])) ## why inverted
      # Store the updated mean values
      mu_j_updated <- AccPGD.Dm(a = a, S = S, lambda = lambda,
                                M.init = as.vector(t(M[j,c,])),
                                tol=1e-12)
      mu_j_mat[j,c , ] <- mu_j_updated
    }
    # print(j)
  }

  M_out <- array(0, dim = c(p, C_dim, V_dim))
  for (j in 1:p) {
    for (c in 1:C_dim) {
      M_out[j, c, ] <- mu_j_mat[j, c, ]
    }
  }
  return(M_out)
}

update_mu_nopenalty <- function(Y, M, sigma2, W){
  M.out <- array(0, c(dim(Y)[2], dim(W)[2], dim(W)[3]))
  for(j in 1:dim(M)[1]){
    M.out[j,,] <- apply(W*Y[,j], c(2,3), sum)/apply(W, c(2,3), sum)
  }
  return(M.out)
}

M.step.variance <- function(Y, W, M){
  p <- dim(Y)[2]
  sigma2 <- array(0, dim=dim(M))
  W.tot <- apply(W, c(2,3), sum)
  for(kk in 1:p){
    for(c.ind in 1:dim(M)[2]){
      for(v.ind in 1:dim(M)[3]){
        #weighted_mean <- sum(Y[,kk]*W[,c.ind,v.ind])/W.tot[c.ind, v.ind]
        sigma2[kk, c.ind, v.ind] <- max(sum(W[, c.ind, v.ind] * (Y[, kk] - M[kk, c.ind, v.ind])^2) / W.tot[c.ind, v.ind], 1e-10)
      }
    }
  }
  return(sigma2)
}

M.step.probs <- function(Y, W){
  apply(W, c(2,3), sum)/dim(Y)[1]
}

