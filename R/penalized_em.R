#' Penalized observed-data objective
#'
#' The fitted penalty is n * lambda times the Euclidean norm of each
#' condition-specific mean vector after centering across viral states.
#' @keywords internal
penalized_em_objective <- function(Y, c_obs, v_obs, M, sigma2, probs, lambda) {
  centered <- sweep(M, c(1, 2), apply(M, c(1, 2), mean))
  penalty <- sum(sqrt(apply(centered^2, c(1, 2), sum)))
  -observed.data.loglik(Y, M, sigma2, probs, c_obs, v_obs) +
    nrow(Y) * lambda * penalty
}

#' Run penalized generalized EM and retain every iterate
#' @keywords internal
run_penalized_em <- function(Y, c_obs, v_obs, M, sigma2, probs, lambda,
                             max.iter, tol) {
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite, nonnegative number.")
  }
  if (length(max.iter) != 1L || !is.finite(max.iter) || max.iter < 1 ||
      max.iter != as.integer(max.iter)) {
    stop("max.iter must be a positive integer.")
  }
  objective_trace <- numeric(max.iter + 1L)
  objective_trace[1L] <- penalized_em_objective(
    Y, c_obs, v_obs, M, sigma2, probs, lambda)

  for (mm in seq_len(max.iter)) {
    W <- E_step(Y, c_obs, v_obs, M, probs, sigma2)
    probs.new <- M_step_probs(Y, W)
    M.new <- update_mu(Y, M, sigma2, W, lambda = lambda * nrow(Y))
    sigma2.new <- M_step_variance(Y, W, M.new)
    rel_change <- sum((M - M.new)^2) / max(sum(M^2), .Machine$double.eps)

    # Every computed iterate is retained, including the converged iterate.
    M <- M.new
    sigma2 <- sigma2.new
    probs <- probs.new
    objective_trace[mm + 1L] <- penalized_em_objective(
      Y, c_obs, v_obs, M, sigma2, probs, lambda)

    if (mm == 1L || mm %% 10L == 0L) {
      message(sprintf("  iteration %d/%d (rel_change=%.3e)", mm, max.iter, rel_change))
    }
    if (rel_change < tol) {
      message(sprintf("  converged at iteration %d (rel_change=%.3e)", mm, rel_change))
      break
    }
  }

  list(M = M, sigma2 = sigma2, probs = probs,
       weights = E_step(Y, c_obs, v_obs, M, probs, sigma2),
       objective_trace = objective_trace[seq_len(mm + 1L)])
}
