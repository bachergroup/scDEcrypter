
#' Internal: Proximal operator for L2 penalty
#'
#' Proximal operator used in accelerated proximal gradient descent.
#' This is an internal function and not intended to be called directly by users.
#'
#' @param input Numeric vector.
#' @param lambda Penalty parameter (numeric scalar).
#' @param D Difference matrix.
#' @param Us Matrix from SVD decomposition.
#'
#' @return Numeric vector after applying proximal operator.
#' @keywords internal
#' 
proxOp1 <- function(input, lambda, D, Us) {
  s <- length(input)
  v <- crossprod(Us, input)[-s]
  t1 <- crossprod(D, input)
  t0 <- sum(t1^2)
  if (t0 <= lambda^2) {
    result <- mean(input)*rep(1, s)
  } else {
    result <- input - crossprod(D, input)/(sqrt(sum(v^2))/lambda)
  }
  return(result)
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

#' Accelerated Proximal Gradient Descent
#'
#' Performs accelerated proximal gradient descent (APGD) with row-wise penalties.
#' Used internally for updating mean parameters in the multiway mixture model.
#'
#' @param a Numeric vector of observed values.
#' @param S Numeric weight vector (same length as \code{a}).
#' @param lambda Penalty parameter (numeric scalar).
#' @param M.init Optional initial vector for iteration. Defaults to zero vector.
#' @param max.iter Maximum number of iterations (default: 1e4).
#' @param tol Convergence tolerance (default: 1e-10).
#'
#' @return Numeric vector of updated parameter estimates.
#' @export
AccPGD.Dm <- function(a, S, lambda, M.init = NULL, max.iter = 1e4, tol = 1e-10)  {
  # ---------------------------------
  # preliminaries
  # --------------------------------
  m.len <- length(a)
  D <- diag(rep(1, m.len)) - matrix(1, m.len, m.len)/m.len
  Us <- svd(D)$u
  
  comp.grad <- function(a, M, S) {
    S*(M-a)
  }
  
  if(is.null(M.init)) {
    M <- rep(0, m.len) 
  } else {
    M <- M.init
  }
  
  alpha <- 1/(max(S))
  M.prev <- M
  converged <- FALSE
  obj.current <- obj.func(M, a, S, D, lambda)[1,1]
  kk <- 1
  
  
  while (!converged) {
    # extrapolate
    M.extra <- M + ((kk - 1)/(kk + 2))*(M - M.prev)
    
    # prox operator of ||DM||_2
    M.new <- proxOp1(input = M.extra - alpha*comp.grad(a, M.extra, S), lambda = alpha*lambda, D = D, Us = Us)
    # epsilon <- 1e-10
    # M.new <- ifelse(is.na(M.new) | M.new == 0, epsilon, M.new)
    obj.new <- obj.func(M.new, a, S, D, lambda)[1,1]
    
    M.prev <- M
    M <- M.new
    
    if(abs(obj.current - obj.new) <= tol * abs(obj.current)) {
      # if(kk > 416) {
      converged <- TRUE
    } else {
      obj.current <- obj.new
    }
    kk <- kk + 1
  }
  
  return(M)
}

