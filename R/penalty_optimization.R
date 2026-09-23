
#' Internal: Proximal operator for L2 penalty
#'
#' Proximal operator used in proximal gradient descent.
#' This is an internal function and not intended to be called directly by users.
#'
#' @param input Numeric vector.
#' @param lambda Penalty parameter (numeric scalar).
#' @param D Difference matrix.
#' @param Us Matrix from SVD decomposition.
#'
#' @return Numeric vector after applying proximal operator.
#' @keywords internal
proxOp1 <- function(input, lambda, D, Us) {
  centered <- as.vector(D %*% input)
  size <- sqrt(sum(centered^2))
  if (size <= lambda) return(rep(mean(input), length(input)))
  input - (lambda / size) * centered
}


#' Internal: Objective function
#'
#' Computes the objective value for the proximal gradient descent algorithm.
#' This is an internal helper function and not intended to be called directly.
#'
#' @param M Numeric vector of current estimates.
#' @param a Numeric vector of observed values.
#' @param S Numeric weight vector.
#' @param D Difference matrix.
#' @param lambda Penalty parameter (numeric scalar).
#'
#' @return Numeric value of the objective function.
#' @keywords internal
obj.func <- function(M, a, S, D, lambda) {
  crossprod(M-a, S*(M-a))/2 + lambda*(sqrt(sum(crossprod(D, M)^2)))
}

#' Proximal Gradient Descent
#'
#' Performs proximal gradient descent with row-wise penalties.
#' Used internally for updating mean parameters in the multiway mixture model.
#'
#' @param a Numeric vector of observed values.
#' @param S Numeric weight vector (same length as \code{a}).
#' @param lambda Penalty parameter (numeric scalar).
#' @param M.init Optional initial vector for iteration. Defaults to zero vector.
#' @param max.iter Maximum number of iterations (default: 100).
#' @param tol Convergence tolerance (default: 1e-10).
#'
#' @return Numeric vector of updated parameter estimates.
#' @export
accpgd_basis_cache <- new.env(parent = emptyenv())

get_accpgd_basis <- function(m.len) {
  key <- as.character(m.len)
  basis <- accpgd_basis_cache[[key]]
  if (is.null(basis)) {
    D <- diag(rep(1, m.len)) - matrix(1, m.len, m.len)/m.len
    basis <- list(D = D, Us = svd(D)$u)
    accpgd_basis_cache[[key]] <- basis
  }
  basis
}

AccPGD.Dm <- function(a, S, lambda, M.init = NULL, max.iter = 100, tol = 1e-10)  {
  # ---------------------------------
  # preliminaries
  # --------------------------------
  m.len <- length(a)
  basis <- get_accpgd_basis(m.len)
  D <- basis$D
  Us <- basis$Us
  
  comp.grad <- function(a, M, S) {
    S*(M-a)
  }
  
  if(is.null(M.init)) {
    M <- rep(0, m.len) 
  } else {
    M <- M.init
  }
  
  if (any(!is.finite(c(a, S, lambda, M))) || any(S < 0) || lambda < 0) {
    stop("Mean update requires finite inputs, nonnegative weights, and nonnegative lambda.")
  }
  if (max(S) == 0) return(rep(mean(M), m.len))
  alpha <- 1 / max(S)
  # A proximal-gradient step majorizes the quadratic when alpha <= 1/max(S).
  # Every step decreases the M-step objective in exact arithmetic.
  for (kk in seq_len(max.iter)) {
    M.new <- proxOp1(M - alpha * comp.grad(a, M, S), alpha * lambda, D, Us)
    change <- sum((M.new - M)^2)
    M <- M.new
    if (change <= tol^2 * max(1, sum(M^2))) break
  }

  M
}
