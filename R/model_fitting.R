

#' Multiway Mixture Model with Penalty
#'
#' Fits a multiway mixture model using the EM algorithm with
#' an \eqn{\ell_2}-type penalty on mean parameters to encourage
#' similarity across variants.
#'
#' @param Data Numeric matrix of observed data with dimensions
#'   \code{(n x p)}, where \code{n} is the number of samples (cells)
#'   and \code{p} is the number of features (e.g., genes).
#' @param c_obs Integer vector of length \code{n}, giving the pre-labeled
#'   partitioning variable status for each cell (may contain \code{NA}).
#' @param v_obs Integer vector of length \code{n}, giving the pre-labeled
#'   infection status labels for each cell (may contain \code{NA}).
#' @param max.iter Integer. Maximum number of EM iterations for each
#'   value of \code{lambda}.
#' @param tol Numeric. Convergence tolerance for relative change in
#'   mean estimates.
#' @param c_star Integer. Number of condition clusters to fit.
#' @param v_star Integer. Number of variant clusters to fit.
#' @param lambda.vec Numeric vector of per-cell penalty parameters. For
#'   \eqn{n} cells, each fit maximizes \eqn{\ell_n - n\lambda P(\mu)}.
#'   A separate model is fitted for each value.
#' @param infectionLabels Character vector of length equal to the number of viral types,
#'   specifying biological names for each infection/variant status. 
#'   Example: \code{c("Uninfected", "Infected", "Bystander", etc)}
#' @param partitionLabels Character vector of length equal to the number of cell types,
#'   specifying biological names for each cell type cluster.
#'   Example: \code{c("T_cells", "B_cells", "Monocytes", etc)}
#' @return A list containing:
#' \describe{
#'   \item{\code{M_generation}}{Estimated means (\code{p x C x V})
#'   for each \code{lambda}.}
#'   \item{\code{sigma2_generation}}{Estimated variances
#'   (\code{p x C x V}) for each \code{lambda}.}
#'   \item{\code{probs_generation}}{Estimated mixing proportions
#'   (\code{C x V}) for each \code{lambda}.}
#'   \item{\code{weights_generation}}{Posterior responsibilities
#'   (\code{n x C x V}) for each \code{lambda}.}
#'   \item{\code{objective_generation}}{Penalized observed-data objective
#'   at initialization and after every retained EM iteration.}
#' }
#'
#' @export
fit_scDEcrypter <- function(Data, c_obs=NULL, v_obs=NULL, infectionLabels=NULL, partitionLabels=NULL,
                      max.iter=250, tol=1e-8, c_star, v_star, lambda.vec){


  if (!"RNA" %in% names(Data@assays)) {
      stop("Seurat object must contain an RNA assay.")
  }
 
  if(is.null(c_obs)) c_obs <- Data$C.preLabel
  if(is.null(v_obs)) v_obs <- Data$V.preLabel

  # Extract and rotate data (as.matrix: layer may be sparse; EM needs dense)
  Y <- t(as.matrix(Data[["RNA"]]$data.Generation))

  M_lambda_list <- list()
  sigma2_lambda_list <- list()
  probs_lambda_list <- list()
  weights_lambda_list <- list()
  objective_lambda_list <- list()
  message("Initializing...")
  init_t0 <- proc.time()[3]
  tmp <- initializer_scDEcrypter(Y, c_obs, v_obs, max.iter, tol, c_star, v_star)
  message(sprintf("Initialization done (%.1fs)", proc.time()[3] - init_t0))
  
  
  n_lambda <- length(lambda.vec)
  for(lambda_idx in seq_along(lambda.vec)) {
    lambda <- lambda.vec[[lambda_idx]]
    lambda_t0 <- proc.time()[3]
    message(sprintf("[lambda %d/%d] Running lambda=%s", lambda_idx, n_lambda, format(lambda, scientific = TRUE)))
    
    message("EM in progress...")
    fit <- run_penalized_em(Y, c_obs, v_obs, tmp$M, tmp$sigma2,
                            tmp$probs, lambda, max.iter, tol)
    message("Preparing outputs...")
    M.new <- fit$M
    sigma2.new <- fit$sigma2
    probs.new <- fit$probs
    out.weights <- fit$weights
    
    dimnames(M.new) <- list(dimnames(M.new)[[1]], 
                               paste0("partitionStatus_", seq_len(c_star)),
                               paste0("infectionStatus_", seq_len(v_star)))
    dimnames(sigma2.new) <- dimnames(M.new)
    dimnames(out.weights) <- list(dimnames(out.weights)[[1]],
                                 paste0("partitionStatus_", seq_len(c_star)),
                                 paste0("infectionStatus_", seq_len(v_star)))
    dimnames(probs.new) <- list(paste0("partitionStatus_", seq_len(c_star)),
                                 paste0("infectionStatus_", seq_len(v_star)))
                                 
    M_lambda_list[[as.character(lambda)]] <- M.new
    sigma2_lambda_list[[as.character(lambda)]] <- sigma2.new
    probs_lambda_list[[as.character(lambda)]] <- probs.new
    weights_lambda_list[[as.character(lambda)]] <- out.weights
    objective_lambda_list[[as.character(lambda)]] <- fit$objective_trace
    message(sprintf("[lambda %d/%d] Complete (%.1fs)", lambda_idx, n_lambda, proc.time()[3] - lambda_t0))
  }
  
  if (length(lambda.vec) == 1) {
    out.list <- list(
      "M_generation" = M_lambda_list[[1]],
      "sigma2_generation" = sigma2_lambda_list[[1]], 
      "probs_generation" = probs_lambda_list[[1]],
      "weights_generation" = weights_lambda_list[[1]],
      "objective_generation" = objective_lambda_list[[1]]
    )
    if(!is.null(infectionLabels) & !is.null(partitionLabels)) {
     out.list <- switchNames(out.list, infectionLabels, partitionLabels)
    }
    return(out.list)
  } else {                             
    return(list(
      "M_generation" = M_lambda_list,
      "sigma2_generation" = sigma2_lambda_list, 
      "probs_generation" = probs_lambda_list,
      "weights_generation" = weights_lambda_list,
      "objective_generation" = objective_lambda_list
    ))
  }

}
